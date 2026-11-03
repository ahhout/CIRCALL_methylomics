#!/usr/bin/env python3
"""Check regional PAT counts against hand calculations and native wgbstools.

This is a small synthetic validation, not a patient/reference pipeline run.
The fixture deliberately includes a fragment spanning more than 100 CpGs,
adjacent and separated regions, missing calls, repeated rows and short reads.
"""

import argparse
import ast
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import output_guard


def utc():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def require_run_name(path):
    match = re.fullmatch(r"(\d{4}-\d{2}-\d{2}_\d{2}-\d{2})-.+-v\d+", path.name)
    if not match:
        raise ValueError("Output name must be YYYY-MM-DD_HH-MM-purpose-vN (UTC)")
    datetime.strptime(match[1], "%Y-%m-%d_%H-%M")


def fixture():
    regions = [
        ("chr1", 1000, 1050, 110, 116),
        ("chr1", 1100, 1150, 116, 122),
        ("chr1", 2000, 2100, 250, 256),
        ("chr2", 3000, 3050, 1000, 1006),
        ("chr2", 4000, 4050, 1200, 1206),
    ]
    long_pattern = "." * 109 + "CTC.TC" + "TTCTCC" + "." * 128 + "CCCTTT"
    records = [
        ("chr1", 1, long_pattern, 2),
        ("chr1", 108, "CCCTTTCC", 3),
        ("chr1", 110, "T", 7),
        ("chr1", 112, "....", 5),
        ("chr1", 114, "CTCTCTCT", 4),
        ("chr1", 116, "CCCT", 5),
        ("chr1", 118, "CTTT", 3),
        ("chr1", 118, "CTTT", 3),
        ("chr1", 121, "TCC", 11),
        ("chr1", 249, "TCCCTTTT", 6),
        ("chr2", 1000, "CCCCCC", 13),
        ("chr2", 1002, "T..T", 17),
        ("chr2", 1206, "CCC", 100),
    ]
    # Every value here is independently hand-calculated from the fixture.
    expected = [
        {"methylated": 19, "total": 43, "U": 0, "X": 5, "M": 0},
        {"methylated": 39, "total": 91, "U": 6, "X": 6, "M": 5},
        {"methylated": 24, "total": 48, "U": 0, "X": 8, "M": 0},
        {"methylated": 78, "total": 112, "U": 0, "X": 0, "M": 13},
        {"methylated": 0, "total": 0, "U": 0, "X": 0, "M": 0},
    ]
    return regions, records, expected


def independent_counts(regions, records):
    """Deliberately simple per-character oracle; no production helper code."""
    result = []
    for chromosome, _, _, first, last in regions:
        methylated = total = 0
        for chrom, start, pattern, multiplicity in records:
            if chrom != chromosome:
                continue
            for offset, call in enumerate(pattern):
                if first <= start + offset < last and call in "CT":
                    total += multiplicity
                    methylated += multiplicity * (call == "C")
        result.append({"methylated": methylated, "total": total})
    return result


def native_counts(root, regions, text, tools):
    result = []
    counter = tools / "pat2beta/stdin2beta"
    for index, (_, _, _, first, last) in enumerate(regions):
        command = [str(counter), str(first), str(last)]
        completed = subprocess.run(command, input=text, text=True,
                                   capture_output=True, check=True)
        values = [int(value) for value in completed.stdout.split()]
        if len(values) != (last - first) * 2:
            raise AssertionError("Native pat2beta output has incorrect dimensions")
        result.append({"methylated": sum(values[::2]), "total": sum(values[1::2])})
        (root / f"native-pat2beta-{index}.txt").write_text(completed.stdout)
    return result


def check_celfie(source):
    """Run native functions without its CLI or unavailable optional imports."""
    import numpy as np

    tree = ast.parse(source.read_text())
    functions = ast.Module(body=[node for node in tree.body
                                if isinstance(node, ast.FunctionDef)], type_ignores=[])
    namespace = {"np": np}
    exec(compile(functions, str(source), "exec"), namespace)
    x = np.array([[5., 8., 1.]])
    depth = np.array([[10., 10., 10.]])
    reference = np.array([[0., 1., 0.], [1., 0., 1.]])
    reference_depth = np.array([[0., 2., 1.], [2., 2., 1.]])
    with np.errstate(divide="ignore", invalid="ignore"):
        alpha, gamma, likelihood = namespace["em"](
            x, depth, reference, reference_depth, 2, .0001, 42)
    assert reference[0, 0] == 1 and reference_depth[0, 0] == 2
    assert np.isfinite(alpha).all() and np.isfinite(gamma).all()
    assert np.isfinite(likelihood)
    return {"source": str(source), "sha256": sha(source),
            "zero_depth_reference_before": [0, 0],
            "zero_depth_reference_after_native_initialization": [1, 2],
            "finite_after_two_native_em_iterations": True,
            "preprocessing_imputation": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tools-dir", type=Path, default=Path(
        "/path/to/hpc-home/.conda/envs/nextflow"))
    parser.add_argument("--celfie-source", type=Path, default=Path(__file__).resolve().parents[1]
                        / "software/DecoNFlow-4195f65ad986ea5a68129d2da5b2b1effed25bc8/bin/celfie.py")
    args = parser.parse_args()
    require_run_name(args.output)
    revision = output_guard.require_clean_source(Path(__file__))
    started = utc()
    output_guard.reserve_path(args.output, {'project_revision': revision, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    manifest = {"started_utc": started, "status": "running",
                "scope": "synthetic count and native pseudocount validation",
                "script": str(Path(__file__).resolve()), "script_sha256": sha(__file__),
                "argv": sys.argv, "tools_dir": str(args.tools_dir), "checks": {}}
    write_json(args.output / "execution.json", manifest)
    try:
        from pat_counts import count_pat

        regions, records, expected = fixture()
        text = "".join("\t".join(map(str, row)) + "\n" for row in records)
        pat = args.output / "fixture.pat.gz"
        with gzip.open(pat, "wt") as stream:
            stream.write(text)
        atlas = args.output / "fixture.atlas.tsv"
        with atlas.open("w") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(["chr", "start", "end", "startCpG", "endCpG"])
            writer.writerows(regions)
        oracle = independent_counts(regions, records)
        native_source = args.tools_dir / "share/wgbs_tools/src"
        native = native_counts(args.output, regions, text, native_source)
        assert oracle == [{key: row[key] for key in ("methylated", "total")}
                          for row in expected]
        assert native == oracle
        output = args.output / "candidate"
        count_pat(pat, atlas, output, args.tools_dir)
        count_path = output / "counts.tsv"
        with count_path.open() as stream:
            measured = [{key: int(row[key]) for key in expected[0]}
                        for row in csv.DictReader(stream, delimiter="\t")]
        assert measured == expected, (measured, expected)
        first_execution = json.loads((output / "execution.json").read_text())
        first_counts_sha = sha(count_path)
        count_pat(pat, atlas, output, args.tools_dir)
        assert json.loads((output / "execution.json").read_text()) == first_execution
        assert sha(count_path) == first_counts_sha
        corrupt = args.output / "truncated.pat.gz"
        corrupt.write_bytes(pat.read_bytes()[:-8])
        corrupt_output = args.output / "truncated-input"
        try:
            count_pat(corrupt, atlas, corrupt_output, args.tools_dir)
        except ValueError as error:
            assert "PAT counting pipeline failed" in str(error)
        else:
            raise AssertionError("Truncated gzip was accepted")
        assert json.loads((corrupt_output / "execution.json").read_text())["status"] == "failed"
        assert not (corrupt_output / "counts.tsv").exists()
        manifest["checks"] = {
            "hand_calculations_equal_independent_character_counter": True,
            "independent_counter_equal_native_pat2beta": True,
            "candidate_equal_hand_calculations_including_native_UXM": True,
            "completed_cache_reuses_unchanged_counts_and_execution_manifest": True,
            "truncated_gzip_is_rejected_without_completed_count_table": True,
            "synthetic_regions": len(regions), "synthetic_PAT_records": len(records),
            "longest_fragment_span_CpGs": max(len(row[2]) for row in records),
            "expected": expected, "measured": measured,
            "native_celfie": check_celfie(args.celfie_source),
        }
        manifest["candidate_helper"] = first_execution["settings"]
        manifest["native_sources"] = {
            name: {"path": str(native_source / name), "sha256": sha(native_source / name)}
            for name in ("pat2beta/stdin2beta", "pat2beta/stdin2beta.cpp", "homog/homog", "homog/homog.cpp")}
        manifest["status"] = "completed"
        manifest["completed_utc"] = utc()
        write_json(args.output / "execution.json", manifest)
        output_guard.finish(args.output, manifest["status"] == "completed")
        print(json.dumps({"status": manifest["status"], "output": str(args.output)}))
    except Exception as exc:
        manifest.update(status="failed", failed_utc=utc(), error=f"{type(exc).__name__}: {exc}")
        write_json(args.output / "execution.json", manifest)
        output_guard.finish(args.output, manifest["status"] == "completed")
        raise


if __name__ == "__main__":
    main()
