#!/usr/bin/env python3
"""Validate a portable PAT bundle against an explicit methylseq run and samplesheet."""
import argparse
import copy
import csv
import json
from pathlib import Path
from prepare_pat import (HG38_CPG_SHA256, LineageInputError, discover_bams,
                         read_bundle, reference_cpg_count, sha256, validate_reuse)


def validate(args):
    with args.samplesheet.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    samples = [row["sample"] for row in rows]
    if not samples or len(samples) != len(set(samples)):
        raise LineageInputError("Samplesheet must contain unique, nonempty sample IDs")
    bams = discover_bams(args.results_dir, rows)
    provenance, manifest = read_bundle(args.pat_dir)
    if list(manifest) != samples:
        raise LineageInputError("Bundle must match the complete samplesheet in order")
    reference_dir = args.reference_root / "references/hg38"
    reference = {"cpg_bed_sha256": sha256(reference_dir / "CpG.bed.gz"),
                 "cpg_count": reference_cpg_count(reference_dir)}
    if reference != {"cpg_bed_sha256": HG38_CPG_SHA256, "cpg_count": 29401795}:
        raise LineageInputError("Unexpected hg38 CpG dictionary")
    for sample, bam in bams.items():
        validate_reuse(sample, bam, args.pat_dir, provenance, manifest, reference,
                       {"wgbstools": "wgbstools version 0.3.0", "samtools": "samtools 1.22.1"},
                       args.tabix, hash_bam=args.hash_bams)
    result = {"status": "PASS", "sample_order": samples,
              "results_dir": str(args.results_dir.resolve()),
              "pat_dir": str(args.pat_dir.resolve()), "reference": reference,
              "provenance_sha256": sha256(args.pat_dir / "pat_provenance.json"),
              "manifest_sha256": sha256(args.pat_dir / "pat_manifest.tsv"),
              "bam_content_rehashed": args.hash_bams}
    if args.compare_with:
        other_args = copy.copy(args)
        other_args.pat_dir, other_args.compare_with = args.compare_with, None
        validate(other_args)
        _, other_manifest = read_bundle(args.compare_with)
        for sample in samples:
            for kind in ("pat", "csi", "lbeta"):
                first = (args.pat_dir / manifest[sample][kind]).resolve()
                second = (args.compare_with / other_manifest[sample][kind]).resolve()
                if first != second:
                    raise LineageInputError(f"Full/heme inputs do not share the same artifact: {sample} {kind}")
        result["same_artifacts_as"] = str(args.compare_with.resolve())
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--pat-dir", type=Path, required=True, help="Directory containing pat_manifest.tsv")
    parser.add_argument("--samplesheet", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path,
                        default=Path("/path/to/hpc-home/.conda/envs/nextflow/share/wgbs_tools"))
    parser.add_argument("--tabix", type=Path, default=Path("/path/to/hpc-home/.conda/envs/nextflow/bin/tabix"))
    parser.add_argument("--compare-with", type=Path, help="Require a second validated run to use the same PAT artifacts")
    parser.add_argument("--hash-bams", action="store_true", help="Rehash whole BAMs; use Slurm for large files")
    try:
        validate(parser.parse_args())
    except (LineageInputError, OSError, ValueError, KeyError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc


if __name__ == "__main__":
    main()
