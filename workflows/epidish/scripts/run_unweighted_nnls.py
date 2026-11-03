#!/usr/bin/env python3

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import nnls
from recovery_contract import PREPARED_FILES, atomic_json, check_receipt, check_recovery_outputs, write_receipt, utc_now


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run unweighted NNLS on prepared EpiDISH comparison matrices."
    )
    parser.add_argument("work_dir", nargs="?", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test == (args.work_dir is not None):
        parser.error("provide either WORK_DIR or --self-test")
    return args


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def solve_unweighted(reference, observed):
    coefficients, _ = nnls(reference, observed)
    total = coefficients.sum()
    if not np.isfinite(total) or total <= 0:
        raise ValueError("unweighted NNLS returned a zero or non-finite composition")
    return coefficients / total


def validate_fractions(values, label):
    if not np.isfinite(values).all() or (values < -1e-10).any():
        raise ValueError(f"invalid fractions for {label}")
    if not np.allclose(values.sum(axis=0), 1.0, atol=1e-4, rtol=0):
        raise ValueError(f"fractions do not sum to one for {label}")


def atomic_tsv(frame, path):
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(temporary, sep="\t", index=False)
    os.replace(temporary, path)


def run_workflow(work_dir):
    started_utc = utc_now()
    work_dir = work_dir.resolve(strict=True)
    check_recovery_outputs(work_dir)
    required = ['provenance/comparison_jobs.tsv', 'provenance/pat_manifest.tsv',
                'provenance/sample_metadata.tsv', 'provenance/input_manifest.json']
    if not check_receipt(work_dir, 'preparation', required):
        raise ValueError('Missing preparation receipt')
    timing_path = work_dir / 'provenance/nnls_execution_attempts.json'
    timing = json.loads(timing_path.read_text()) if timing_path.exists() else {
        'started_utc': started_utc, 'attempts': []}
    attempt = {'attempt_id': len(timing['attempts']) + 1, 'started_utc': started_utc, 'status': 'running'}
    timing['attempts'].append(attempt)
    atomic_json(timing_path, timing)
    try:
        _run_workflow(work_dir, started_utc)
    except BaseException:
        attempt['status'] = 'failed'
        attempt['ended_utc'] = utc_now()
        atomic_json(timing_path, timing)
        raise
    attempt.update(status='completed', completed_utc=utc_now())
    attempt['ended_utc'] = attempt['completed_utc']
    atomic_json(timing_path, timing)


def _run_workflow(work_dir, started_utc):
    work_dir = work_dir.resolve(strict=True)
    jobs_path = work_dir / "provenance" / "comparison_jobs.tsv"
    pats_path = work_dir / "provenance" / "pat_manifest.tsv"
    jobs = pd.read_csv(jobs_path, sep="\t")
    samples = pd.read_csv(pats_path, sep="\t")["sample"].tolist()
    if len(jobs) != 6 or len(samples) == 0 or len(set(samples)) != len(samples):
        raise ValueError("prepared jobs or sample manifest is invalid")

    outputs = []
    for job in jobs.itertuples(index=False):
        comparison_dir = (work_dir / job.comparison_dir).resolve(strict=True)
        relative_output = f"{job.comparison_dir}/nnls_unweighted_fractions.tsv"
        stage = f"nnls/{job.profile}/{job.threshold_label}"
        if not check_receipt(work_dir, f"prepared/{job.profile}/{job.threshold_label}",
                             [f"{job.comparison_dir}/{name}" for name in PREPARED_FILES]):
            raise ValueError("Missing prepared comparison receipt")
        cached = check_receipt(work_dir, stage, [relative_output])
        if cached:
            outputs.append({"profile": job.profile, "coverage_threshold": int(job.coverage_threshold),
                            "path": relative_output, "sha256": sha256_file(work_dir / relative_output)})
            continue
        reference_frame = pd.read_csv(comparison_dir / "reference_matrix.tsv", sep="\t")
        sample_frame = pd.read_csv(
            comparison_dir / "sample_regional_values.tsv.gz", sep="\t"
        )
        if reference_frame["region_id"].tolist() != sample_frame["region_id"].tolist():
            raise ValueError(f"reference/sample row mismatch: {comparison_dir}")
        cell_types = reference_frame.columns[1:].tolist()
        if sample_frame.columns.tolist() != ["region_id", *samples]:
            raise ValueError(f"sample order mismatch: {comparison_dir}")
        reference = reference_frame[cell_types].to_numpy(dtype=float)
        observed = sample_frame[samples].to_numpy(dtype=float)
        if not np.isfinite(reference).all() or not np.isfinite(observed).all():
            raise ValueError(f"non-finite NNLS input: {comparison_dir}")
        estimates = np.column_stack(
            [solve_unweighted(reference, observed[:, index]) for index in range(len(samples))]
        )
        validate_fractions(estimates, str(comparison_dir))
        output = pd.DataFrame({"CellType": cell_types})
        for index, sample in enumerate(samples):
            output[sample] = estimates[:, index]
        output_path = comparison_dir / "nnls_unweighted_fractions.tsv"
        if not cached:
            atomic_tsv(output, output_path)
            write_receipt(work_dir, stage, [relative_output])
        outputs.append(
            {
                "profile": job.profile,
                "coverage_threshold": int(job.coverage_threshold),
                "path": relative_output,
                "sha256": sha256_file(output_path),
            }
        )

    provenance = {
        "started_utc": started_utc,
        "completed_utc": utc_now(),
        "method": "scipy.optimize.nnls on unweighted marker rows; normalized to sum one",
        "scipy_version": scipy.__version__,
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "outputs": outputs,
    }
    provenance_path = work_dir / "provenance" / "unweighted_nnls.json"
    if provenance_path.exists():
        saved = json.loads(provenance_path.read_text())
        if saved["outputs"] != outputs:
            raise ValueError("Recorded unweighted NNLS outputs differ")
        print(f"Reused unweighted NNLS for {len(outputs)} comparisons")
        return
    temporary = provenance_path.with_name(provenance_path.name + ".tmp")
    temporary.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, provenance_path)
    write_receipt(work_dir, "nnls_provenance", ["provenance/unweighted_nnls.json"])
    print(f"Completed unweighted NNLS for {len(outputs)} comparisons")


def self_test():
    rng = np.random.default_rng(20260907)
    reference = rng.uniform(0.05, 0.95, size=(200, 5))
    expected = np.array([0.05, 0.15, 0.30, 0.20, 0.30])
    observed = reference @ expected
    estimated = solve_unweighted(reference, observed)
    if np.max(np.abs(estimated - expected)) >= 0.02:
        raise AssertionError("unweighted NNLS failed noiseless synthetic recovery")
    print("Unweighted NNLS synthetic test passed")


if __name__ == "__main__":
    parsed = parse_args()
    if parsed.self_test:
        self_test()
    else:
        run_workflow(parsed.work_dir)
