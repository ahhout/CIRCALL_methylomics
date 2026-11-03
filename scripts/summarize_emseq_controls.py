#!/usr/bin/env python3
"""Summarize EM-seq control methylation from Bismark coverage files."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / 'scripts/lib'))
import output_guard

CONTROLS = ("phage_lambda", "plasmid_puc19c")


@dataclass
class Counts:
    covered_cpgs: int = 0
    methylated: int = 0
    unmethylated: int = 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a per-sample QC table for lambda and pUC19 controls."
    )
    parser.add_argument(
        "--coverage-dir",
        type=Path,
        required=True,
        help="Directory containing one Bismark .cov.gz file per sample.",
    )
    parser.add_argument(
        "--samplesheet",
        type=Path,
        required=True,
        help="nf-core/methylseq samplesheet containing a sample column.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output TSV path.",
    )
    parser.add_argument('--owned-output-root', type=Path, required=True,
                        help='Exact owned incomplete methylseq result bundle.')
    parser.add_argument('--attempt-id', required=True,
                        help='Active attempt ID recorded by the project launcher.')
    return parser.parse_args()


def validate_owner(output: Path, root: Path, attempt_id: str) -> dict:
    """Guard standalone invocation before reading payloads or replacing any file."""
    record = output_guard.read_record(root)
    if root != Path(record['output']):
        raise ValueError('Owned output root must be its canonical absolute path')
    if record['status'] != 'running':
        raise ValueError('Control QC requires an owned running attempt; completed/historical results are immutable')
    identity = record['identity']
    if identity.get('project_source') != str(PROJECT):
        raise ValueError('Output belongs to a different project source')
    for repository, key in ((PROJECT, 'project_commit'),
                            (Path(identity['pipeline_source']), 'pipeline_commit')):
        prefix = ['git', '--no-optional-locks', '-C', str(repository)]
        current = subprocess.check_output(prefix + ['rev-parse', 'HEAD'], text=True).strip()
        dirty = subprocess.check_output(prefix + ['status', '--porcelain', '--untracked-files=all'], text=True)
        if current != identity.get(key) or dirty:
            raise ValueError('Source identity changed before control QC')
    if not re.fullmatch(r'[0-9]+-[0-9a-f]{12}', attempt_id):
        raise ValueError('Invalid attempt identity')
    if not record['attempts'] or record['attempts'][-1].get('attempt_id') != attempt_id:
        raise ValueError('Control QC attempt differs from the current owned attempt')
    expected = root / 'control_qc/control_qc.tsv'
    if output != expected or output.is_symlink() or output.parent.is_symlink():
        raise ValueError('Control QC must publish to the owned bundle without symlinks')
    return record


def read_samples(samplesheet: Path) -> list[str]:
    if not samplesheet.is_file():
        raise ValueError(f"samplesheet does not exist: {samplesheet}")
    with samplesheet.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "sample" not in reader.fieldnames:
            raise ValueError(f"samplesheet has no 'sample' column: {samplesheet}")
        samples = [row["sample"].strip() for row in reader]
    if not samples or any(not sample for sample in samples):
        raise ValueError(f"samplesheet contains no samples or a blank sample: {samplesheet}")
    if len(samples) != len(set(samples)):
        raise ValueError(f"samplesheet contains duplicate sample identifiers: {samplesheet}")
    return samples


def assign_coverage_files(coverage_dir: Path, samples: list[str]) -> dict[str, Path]:
    if not coverage_dir.is_dir():
        raise ValueError(f"coverage directory does not exist: {coverage_dir}")
    files = sorted(coverage_dir.glob("*.cov.gz"))
    assignments: dict[str, Path] = {}
    for path in files:
        matches = [
            sample
            for sample in samples
            if path.name.startswith(f"{sample}_") or path.name.startswith(f"{sample}.")
        ]
        if not matches:
            continue
        sample = max(matches, key=len)
        if sample in assignments:
            raise ValueError(
                f"multiple coverage files match sample {sample}: "
                f"{assignments[sample]} and {path}"
            )
        assignments[sample] = path
    missing = [sample for sample in samples if sample not in assignments]
    if missing:
        raise ValueError("missing coverage files for samples: " + ", ".join(missing))
    return assignments


def summarize_file(path: Path) -> dict[str, Counts]:
    counts = {control: Counts() for control in CONTROLS}
    try:
        with gzip.open(path, mode="rt", encoding="utf-8", newline="") as handle:
            for line_number, line in enumerate(handle, start=1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    raise ValueError(f"{path}:{line_number}: expected at least 6 columns")
                control = fields[0]
                if control not in counts:
                    continue
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                    float(fields[3])
                    methylated = int(fields[4])
                    unmethylated = int(fields[5])
                except ValueError as exc:
                    raise ValueError(f"{path}:{line_number}: malformed numeric value") from exc
                if start < 1 or end < start or methylated < 0 or unmethylated < 0:
                    raise ValueError(f"{path}:{line_number}: invalid coordinate or count")
                target = counts[control]
                total = methylated + unmethylated
                if total > 0:
                    target.covered_cpgs += 1
                target.methylated += methylated
                target.unmethylated += unmethylated
    except (OSError, EOFError) as exc:
        raise ValueError(f"cannot read gzip coverage file {path}: {exc}") from exc
    return counts


def metrics(control: str, counts: Counts) -> tuple[str, str, str]:
    total = counts.methylated + counts.unmethylated
    if total == 0:
        return "NA", "NA", "NO_DATA"
    methylation = 100.0 * counts.methylated / total
    if control == "phage_lambda":
        conversion = 100.0 - methylation
        status = "PASS" if methylation <= 0.5 else "WARN"
        return f"{methylation:.6f}", f"{conversion:.6f}", status
    status = "PASS" if methylation >= 95.0 else "WARN"
    return f"{methylation:.6f}", "NA", status


def write_output(
    output: Path,
    samples: list[str],
    summaries: dict[str, dict[str, Counts]],
    *, owned_output_root: Path, attempt_id: str,
) -> None:
    validate_owner(output, owned_output_root, attempt_id)
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".tmp", dir=output.parent, text=True
    )
    try:
        with os.fdopen(fd, mode="w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "sample",
                    "control",
                    "covered_cpgs",
                    "methylated_calls",
                    "unmethylated_calls",
                    "weighted_cpg_methylation_percent",
                    "lambda_conversion_efficiency_percent",
                    "status",
                ]
            )
            for sample in samples:
                for control in CONTROLS:
                    counts = summaries[sample][control]
                    methylation, conversion, status = metrics(control, counts)
                    writer.writerow(
                        [
                            sample,
                            control,
                            counts.covered_cpgs,
                            counts.methylated,
                            counts.unmethylated,
                            methylation,
                            conversion,
                            status,
                        ]
                    )
        validate_owner(output, owned_output_root, attempt_id)
        os.replace(temporary_name, output)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def main() -> int:
    args = parse_args()
    provenance = None
    provenance_created = False
    started = datetime.now(timezone.utc).isoformat()
    try:
        owner = validate_owner(args.output, args.owned_output_root, args.attempt_id)
        submitted = args.owned_output_root / 'pipeline_info/submitted'
        expected_coverage = args.owned_output_root / 'bismark/methylation_calls/methylation_coverage'
        if args.coverage_dir != expected_coverage or args.samplesheet != submitted / 'samplesheet.csv':
            raise ValueError('Control-QC inputs must belong to the same owned output bundle')
        if args.samplesheet.is_symlink() or hashlib.sha256(args.samplesheet.read_bytes()).hexdigest() != owner['identity']['inputs']['samplesheet.csv']:
            raise ValueError('Submitted samplesheet identity changed')
        attempt = submitted / 'attempts' / args.attempt_id
        if not attempt.is_dir() or any(path.is_symlink() for path in (submitted.parent, submitted, attempt.parent, attempt)):
            raise ValueError('Missing or aliased attempt directory')
        provenance = attempt / 'control-qc.json'
        with provenance.open('x') as handle:
            json.dump({'started_utc': started, 'status': 'running', 'attempt_id': args.attempt_id}, handle)
        provenance_created = True
        samples = read_samples(args.samplesheet)
        assignments = assign_coverage_files(args.coverage_dir, samples)
        summaries = {sample: summarize_file(assignments[sample]) for sample in samples}
        write_output(args.output, samples, summaries,
                     owned_output_root=args.owned_output_root, attempt_id=args.attempt_id)
        output_guard_record = {'started_utc': started, 'completed_utc': datetime.now(timezone.utc).isoformat(),
                              'status': 'completed', 'attempt_id': args.attempt_id}
        temporary = provenance.with_suffix('.json.tmp')
        with temporary.open('x') as handle:
            json.dump(output_guard_record, handle)
        os.replace(temporary, provenance)
    except (ValueError, OSError, KeyError, subprocess.CalledProcessError) as exc:
        if provenance_created:
            descriptor, temporary_name = tempfile.mkstemp(prefix='.control-qc-failure-', dir=provenance.parent)
            try:
                with os.fdopen(descriptor, 'w') as handle:
                    json.dump({'started_utc': started, 'ended_utc': datetime.now(timezone.utc).isoformat(),
                               'status': 'failed', 'attempt_id': args.attempt_id}, handle)
                os.replace(temporary_name, provenance)
            finally:
                if os.path.exists(temporary_name):
                    os.unlink(temporary_name)
        print(
            f"Control-QC reporting failed: {exc}. "
            "Any completed methylseq outputs remain intact.",
            file=sys.stderr,
        )
        return 1
    print(f"Control-QC report written: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
