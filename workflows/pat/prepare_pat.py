#!/usr/bin/env python3
"""Validate methylseq BAMs and prepare reusable UXM-compatible PAT inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import signal
import subprocess
from datetime import datetime, timezone
from pathlib import Path


HG38_CPG_SHA256 = "28180d18667dc0a6e66dd4250593e03a8bf10be50e9212aefd4e9282afea8590"
PRIMARY_CONTIGS = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])
BAM2PAT_DEFAULTS = {
    "mapq": 10,
    "exclude_flags": 1796,
    "paired_end_include_flags": 3,
    "min_cpg_per_read": 1,
    "mate_overlap_handling": "wgbstools default mate merge",
}


class LineageInputError(RuntimeError):
    pass


def run(command: list[str | Path], *, capture: bool = True) -> subprocess.CompletedProcess[str]:
    cmd = [str(value) for value in command]
    print("RUN", " ".join(cmd), flush=True)
    try:
        return subprocess.run(
            cmd,
            check=True,
            text=True,
            stdout=subprocess.PIPE if capture else None,
            stderr=subprocess.PIPE if capture else None,
        )
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or exc.stdout or "").strip()
        raise LineageInputError(
            f"Command failed ({exc.returncode}): {' '.join(cmd)}"
            + (f"\n{detail}" if detail else "")
        ) from exc


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_metadata(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise LineageInputError(f"Missing or empty metadata: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows or not {"sample", "patient", "compartment"}.issubset(reader.fieldnames or []):
        raise LineageInputError("Metadata lacks sample, patient, or compartment")
    samples = [row["sample"] for row in rows]
    if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", sample) for sample in samples) or len(samples) != len(set(samples)):
        raise LineageInputError("Metadata sample identifiers must be non-empty and unique")
    return rows


def discover_bams(results_dir: Path, metadata: list[dict[str, str]], *, allow_subset: bool = False) -> dict[str, Path]:
    bam_dir = results_dir / "bismark/deduplicated"
    if not bam_dir.is_dir():
        raise LineageInputError(f"Missing Bismark BAM directory: {bam_dir}")
    available = {
        path.name.removesuffix(".deduplicated.sorted.bam"): path
        for path in bam_dir.glob("*.deduplicated.sorted.bam")
    }
    metadata_order = [row["sample"] for row in metadata]
    selected = {sample: available[sample] for sample in metadata_order if sample in available}
    unknown = sorted(set(available) - set(metadata_order))
    if unknown and not allow_subset:
        raise LineageInputError(f"BAM samples absent from metadata: {', '.join(unknown)}")
    missing = set(metadata_order) - set(available)
    if missing:
        raise LineageInputError(f"Missing BAM samples: {sorted(missing)}")
    if not selected:
        raise LineageInputError("No metadata samples have BAMs in the selected methylseq result")
    for sample, bam in selected.items():
        bai = Path(f"{bam}.bai")
        if bam.stat().st_size == 0 or not bai.is_file() or bai.stat().st_size == 0:
            raise LineageInputError(f"Missing or empty BAM/BAI for {sample}: {bam}")
    return selected


def sample_xm_tags(samtools: Path, bam: Path, max_records: int = 1000) -> tuple[int, int]:
    process = subprocess.Popen(
        [str(samtools), "view", str(bam), "chr1:1-10000000"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert process.stdout is not None
    sampled = tagged = 0
    for line in process.stdout:
        sampled += 1
        tagged += "\tXM:Z:" in line
        if sampled >= max_records:
            process.terminate(); break
    process.stdout.close()
    _, stderr = process.communicate()
    if process.returncode not in (0, -15):
        raise LineageInputError(f"samtools view failed for {bam}: {stderr.strip()}")
    return sampled, tagged


def validate_bam(sample: str, bam: Path, samtools: Path) -> dict[str, object]:
    run([samtools, "quickcheck", "-v", bam])
    header = run([samtools, "view", "-H", bam]).stdout
    if not any(line.startswith("@HD") and "SO:coordinate" in line for line in header.splitlines()):
        raise LineageInputError(f"{sample} BAM is not declared coordinate-sorted")
    idxstats = run([samtools, "idxstats", bam]).stdout
    mapped = sum(
        int(fields[2])
        for line in idxstats.splitlines()
        if len(fields := line.split("\t")) == 4 and fields[0] in PRIMARY_CONTIGS
    )
    sampled, tagged = sample_xm_tags(samtools, bam)
    if mapped == 0 or sampled == 0 or tagged == 0:
        raise LineageInputError(f"{sample} lacks usable primary-contig XM-tagged alignments")
    return {
        "sample": sample,
        "bam": str(bam.resolve()),
        "bam_bytes": bam.stat().st_size,
        "bai": str(Path(f"{bam}.bai").resolve()),
        "primary_mapped_reads": mapped,
        "alignments_sampled": sampled,
        "xm_tagged_alignments": tagged,
    }


def reference_cpg_count(reference_dir: Path) -> int:
    size_table = reference_dir / "CpG.chrome.size"
    if not size_table.is_file():
        raise LineageInputError(f"Missing CpG size table: {size_table}")
    total = 0
    with size_table.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise LineageInputError(f"Malformed CpG size table: {size_table}")
            total += int(fields[1])
    if total <= 0:
        raise LineageInputError("CpG dictionary contains no sites")
    return total


def validate_pat(sample: str, base: Path, tabix: Path, expected_lbeta_bytes: int) -> dict[str, Path]:
    paths = {
        "pat": base / f"{sample}.pat.gz",
        "csi": base / f"{sample}.pat.gz.csi",
        "lbeta": base / f"{sample}.lbeta",
    }
    for label, path in paths.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise LineageInputError(f"Missing or empty {label} for {sample}: {path}")
    if paths["lbeta"].stat().st_size != expected_lbeta_bytes:
        raise LineageInputError(f"Unexpected LBETA byte count for {sample}")
    if not run([tabix, "-l", paths["pat"]]).stdout.strip():
        raise LineageInputError(f"PAT index contains no contigs for {sample}")
    return paths


def file_identity(path: Path, *, hash_content: bool = False) -> dict[str, object]:
    stat = path.stat()
    record = {"path": str(path.resolve()), "bytes": stat.st_size,
              "mtime_ns": stat.st_mtime_ns, "device": stat.st_dev, "inode": stat.st_ino}
    if hash_content:
        record["sha256"] = sha256(path)
        if file_identity(path) != {key: value for key, value in record.items() if key != "sha256"}:
            raise LineageInputError(f"File changed while hashing: {path}")
    return record


def check_identity(path: Path, record: dict, *, hash_content: bool = False, portable: bool = False) -> None:
    current = file_identity(path)
    # st_dev identifies a mount on one host, not a shared filesystem globally.
    # Keep other identity checks and require content proof across device IDs.
    if any(current[key] != record.get(key) for key in current
           if key != "device" and not (portable and key == "path")):
        raise LineageInputError(f"File identity changed: {path}")
    if hash_content or current["device"] != record.get("device"):
        checked = file_identity(path, hash_content=True)
        if {key: value for key, value in checked.items() if key != "sha256"} != current:
            raise LineageInputError(f"File changed while validating: {path}")
        if checked["sha256"] != record.get("sha256"):
            raise LineageInputError(f"Checksum changed: {path}")


def unique_records(rows: list[dict], label: str) -> dict[str, dict]:
    result = {row["sample"]: row for row in rows}
    if len(result) != len(rows) or not result:
        raise LineageInputError(f"Empty or duplicate samples in {label}")
    return result


def read_bundle(directory: Path, *, legacy: bool = False) -> tuple[dict, dict[str, dict]]:
    """Resolve published paths exclusively through the portable manifest."""
    try:
        provenance = json.loads((directory / "pat_provenance.json").read_text())
        with (directory / "pat_manifest.tsv").open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["sample", "pat", "csi", "lbeta"]:
                raise LineageInputError(f"Invalid PAT manifest columns: {directory}")
            rows = list(reader)
        manifest = unique_records(rows, "PAT manifest")
        inputs = unique_records(provenance["bam_inputs"], "BAM provenance")
        artifacts = unique_records(provenance["pat_artifacts"], "PAT provenance")
        if (list(manifest) != provenance["sample_order"] or
                set(manifest) != set(inputs) or set(manifest) != set(artifacts)):
            raise LineageInputError(f"Inconsistent sample sets: {directory}")
        if not legacy and not (directory / "PAT_COMPLETED").is_file():
            raise LineageInputError(f"Incomplete PAT bundle: {directory}")
        if not legacy and provenance.get("schema_version") != 2:
            raise LineageInputError(f"Legacy bundle needs explicit reconciliation: {directory}")
        for sample, row in manifest.items():
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", sample):
                raise LineageInputError(f"Unsafe sample identifier: {sample}")
            for kind, suffix in (("pat", "pat.gz"), ("csi", "pat.gz.csi"), ("lbeta", "lbeta")):
                if row[kind] != f"pat/{sample}.{suffix}":
                    raise LineageInputError(f"Unsafe or unexpected {kind} manifest path: {row[kind]}")
                if not (directory / row[kind]).is_file():
                    raise LineageInputError(f"Stale PAT manifest entry: {directory / row[kind]}")
                if not legacy and artifacts[sample][kind] != row[kind]:
                    raise LineageInputError(f"Manifest/provenance path mismatch: {sample}")
        return provenance, manifest
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise LineageInputError(f"Invalid PAT bundle {directory}: {exc}") from exc


def validate_reuse(sample: str, bam: Path, directory: Path, provenance: dict,
                   manifest: dict, reference: dict, versions: dict, tabix: Path,
                   *, legacy: bool = False, hash_bam: bool = False) -> dict[str, Path]:
    if (provenance.get("genome") != "hg38" or
            provenance.get("bam2pat_defaults") != BAM2PAT_DEFAULTS or
            any(provenance.get("reference", {}).get(key) != reference[key]
                for key in ("cpg_bed_sha256", "cpg_count")) or
            any(provenance.get("versions", {}).get(key) != versions[key]
                for key in ("wgbstools", "samtools"))):
        raise LineageInputError(f"Incompatible reference, conversion settings or versions: {directory}")
    if sample not in manifest:
        raise LineageInputError(f"Reuse bundle lacks required sample {sample}: {directory}")
    original = unique_records(provenance["bam_inputs"], "BAM inputs")[sample]
    bai = Path(f"{bam}.bai")
    if (Path(original["bam"]).resolve() != bam.resolve() or
            Path(original["bai"]).resolve() != bai.resolve() or
            original["bam_bytes"] != bam.stat().st_size):
        raise LineageInputError(f"Source BAM/BAI mismatch for {sample}: {directory}")
    if not legacy:
        expected_source = bam.parent.parent.parent.resolve()
        if Path(provenance["results_dir"]).resolve() != expected_source:
            raise LineageInputError(f"Source methylseq run mismatch for {sample}")
        check_identity(bam, original["bam_identity"], hash_content=hash_bam)
        check_identity(bai, original["bai_identity"], hash_content=True)
    elif provenance.get("conversion_mode") != "rebuild":
        raise LineageInputError("Legacy reconciliation requires original rebuilt PATs")
    paths = validate_pat(sample, directory / "pat", tabix, reference["cpg_count"] * 4)
    artifact = unique_records(provenance["pat_artifacts"], "PAT artifacts")[sample]
    for kind, path in paths.items():
        if kind in ("pat", "lbeta") and path.stat().st_size != artifact[f"{kind}_bytes"]:
            raise LineageInputError(f"Artifact size differs from original provenance: {path}")
        if not legacy:
            check_identity(path, artifact["identities"][kind], hash_content=True, portable=True)
        elif Path(artifact[kind]).name != path.name:
            raise LineageInputError(f"Legacy artifact name mismatch: {path}")
    return paths


def atomic_symlink(source: Path, destination: Path) -> None:
    # symlink creation itself is atomic and must never replace existing data.
    destination.symlink_to(source.resolve())


def prepare(args: argparse.Namespace) -> None:
    started_utc = datetime.now(timezone.utc).isoformat()
    # This is only set after exclusive output reservation; never touch another run.
    args._timing_reserved = False
    try:
        _prepare(args, started_utc)
    except BaseException:
        if args._timing_reserved:
            try:
                record = {"started_utc": started_utc, "ended_utc": datetime.now(timezone.utc).isoformat(),
                          "status": "failed_or_interrupted"}
                (args.output_dir / "pat_execution.json").write_text(json.dumps(record, indent=2) + "\n")
                provenance_path = args.output_dir / "pat_provenance.json"
                if provenance_path.is_file():
                    provenance = json.loads(provenance_path.read_text())
                    provenance.pop("completed_utc", None)
                    provenance.update(record)
                    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
            except (OSError, ValueError):
                pass  # Keep the original failure when timing cannot be saved.
        raise


def _prepare(args: argparse.Namespace, started_utc: str) -> None:
    for executable in (args.wgbstools, args.samtools, args.tabix):
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise LineageInputError(f"Missing executable: {executable}")
    metadata = load_metadata(args.metadata)
    bams = discover_bams(args.results_dir, metadata, allow_subset=args.allow_subset)
    reference_dir = args.reference_root / "references" / args.genome
    executable_root = args.wgbstools.resolve().parents[2]
    if (executable_root / "references" / args.genome).resolve() != reference_dir.resolve():
        raise LineageInputError("wgbstools executable and requested reference root disagree")
    reference = {"cpg_bed": str((reference_dir / "CpG.bed.gz").resolve()),
                 "cpg_bed_sha256": sha256(reference_dir / "CpG.bed.gz"),
                 "cpg_count": reference_cpg_count(reference_dir)}
    if args.genome != "hg38" or reference["cpg_bed_sha256"] != HG38_CPG_SHA256 or reference["cpg_count"] != 29401795:
        raise LineageInputError("PAT conversion requires the canonical wgbstools hg38 dictionary")
    if args.reference_manifest:
        reference["alignment_manifest"] = file_identity(args.reference_manifest, hash_content=True)
    versions = {name: run([executable, "--version"]).stdout.strip().splitlines()[0]
                for name, executable in (("wgbstools", args.wgbstools), ("samtools", args.samtools))}
    if versions != {"wgbstools": "wgbstools version 0.3.0", "samtools": "samtools 1.22.1"}:
        raise LineageInputError(f"Unreviewed conversion tool versions: {versions}")
    reuse_samples = list(bams) if args.conversion_mode == "reuse" else args.reuse_sample
    if len(reuse_samples) != len(set(reuse_samples)) or not set(reuse_samples).issubset(bams):
        raise LineageInputError("Reuse samples must be unique selected BAM sample IDs")
    if args.conversion_mode == "rebuild" and (args.reuse_dir or args.legacy_reuse_dir or reuse_samples):
        raise LineageInputError("Rebuild mode does not accept reuse inputs")
    bundles = [(directory, legacy, *read_bundle(directory, legacy=legacy))
               for directory, legacy in [(p, False) for p in args.reuse_dir] +
                                        [(p, True) for p in args.legacy_reuse_dir]]
    reuse = {}
    for sample in reuse_samples:
        matches = [(directory, legacy, prov, manifest) for directory, legacy, prov, manifest in bundles
                   if sample in manifest]
        if len(matches) != 1:
            raise LineageInputError(f"Expected exactly one explicit reuse source for {sample}; found {len(matches)}")
        directory, legacy, prov, manifest = matches[0]
        paths = validate_reuse(sample, bams[sample], directory, prov, manifest,
                               reference, versions, args.tabix, legacy=legacy)
        reuse[sample] = (paths, directory, legacy, prov)
    if args.preflight_only:
        for sample, bam in bams.items():
            record = validate_bam(sample, bam, args.samtools)
            if sample in reuse and reuse[sample][2]:
                original = unique_records(reuse[sample][3]["bam_inputs"], "legacy BAM inputs")[sample]
                if record["primary_mapped_reads"] != original["primary_mapped_reads"]:
                    raise LineageInputError(f"Legacy source alignment counts changed for {sample}")
        print(json.dumps({"status": "PREFLIGHT_PASS", "samples": list(bams),
                          "reuse": reuse_samples, "generate": [s for s in bams if s not in reuse],
                          "reference": reference, "versions": versions}, indent=2))
        return
    if args.output_dir.is_symlink() or args.output_dir.exists():
        raise LineageInputError(f"Output must be absent: {args.output_dir}")
    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir()  # Atomic exclusive reservation, including empty/dangling collisions.
    # Exclusive reservation: a failed attempt is retained and cannot be reused silently.
    with (args.output_dir / "PAT_STARTED").open("x") as handle:
        handle.write(f"pid={os.getpid()}\n")
    args._timing_reserved = True
    (args.output_dir / "pat_execution.json").write_text(json.dumps({
        "started_utc": started_utc, "status": "running"}, indent=2) + "\n")
    pat_dir = args.output_dir / "pat"
    pat_dir.mkdir()
    integrity = []
    for sample, bam in bams.items():
        record = validate_bam(sample, bam, args.samtools)
        if sample in reuse and reuse[sample][2]:
            original = unique_records(reuse[sample][3]["bam_inputs"], "legacy BAM inputs")[sample]
            if record["primary_mapped_reads"] != original["primary_mapped_reads"]:
                raise LineageInputError(f"Legacy source alignment counts changed for {sample}")
        print(f"Recording BAM/BAI identities: {sample}", flush=True)
        # Reuse validated hashes for schema-2 consumers; no repeated whole-BAM scan.
        if sample in reuse and not reuse[sample][2]:
            old = unique_records(reuse[sample][3]["bam_inputs"], "BAM inputs")[sample]
            record["bam_identity"] = old["bam_identity"]
        else:
            record["bam_identity"] = file_identity(bam, hash_content=True)
        record["bai_identity"] = file_identity(Path(f"{bam}.bai"), hash_content=True)
        integrity.append(record)
    artifacts, commands, reconciliation = [], [], []
    for sample, bam in bams.items():
        if sample in reuse:
            paths, directory, legacy, original = reuse[sample]
            mode = "reused_symlink"
            for path in paths.values():
                atomic_symlink(path, pat_dir / path.name)
            source = {"bundle": str(directory.resolve()),
                      "provenance_sha256": sha256(directory / "pat_provenance.json"),
                      "manifest_sha256": sha256(directory / "pat_manifest.tsv")}
            if legacy:
                prior = unique_records(original["pat_artifacts"], "legacy artifacts")[sample]
                reconciliation.append({"sample": sample, **source,
                    "original_results_dir": original["results_dir"],
                    "verified_results_dir": str(args.results_dir.resolve()),
                    "original_paths": {kind: prior[kind] for kind in paths},
                    "published_paths": {kind: str(path.resolve()) for kind, path in paths.items()},
                    "limitation": "Original provenance has no content hashes. Hashes here establish current identity, not a retrospective conversion-time checksum."})
        else:
            work_dir = args.output_dir / "intermediates" / sample
            work_dir.mkdir(parents=True)
            bam_link = work_dir / f"{sample}.bam"
            bam_link.symlink_to(bam.resolve())
            (work_dir / f"{sample}.bam.bai").symlink_to(Path(f"{bam}.bai").resolve())
            temp_dir = work_dir / "sorttmp"
            temp_dir.mkdir()
            command = [str(args.wgbstools), "bam2pat", str(bam_link), "--genome", args.genome,
                       "--out_dir", str(work_dir), "--threads", str(args.threads),
                       "--lbeta", "--temp_dir", str(temp_dir)]
            commands.append(command)
            run(command, capture=False)
            paths = validate_pat(sample, work_dir, args.tabix, reference["cpg_count"] * 4)
            # Hardlink then unlink only newly generated files on the same filesystem.
            # link() refuses replacement and avoids copying the payload.
            for path in paths.values():
                os.link(path, pat_dir / path.name)
                path.unlink()
            mode, source = "rebuilt", None
        published = validate_pat(sample, pat_dir, args.tabix, reference["cpg_count"] * 4)
        artifacts.append({"sample": sample, "mode": mode, "source": source,
            **{kind: f"pat/{path.name}" for kind, path in published.items()},
            "pat_bytes": published["pat"].stat().st_size,
            "lbeta_bytes": published["lbeta"].stat().st_size,
            "identities": {kind: file_identity(path, hash_content=True) for kind, path in published.items()}})
    # Moving a bundle changes absolute paths of newly built artifacts. Validate those
    # by relative path, size, inode, timestamp and hash; retain initial path as provenance.
    for record in integrity:
        check_identity(Path(record["bam"]), record["bam_identity"])
        check_identity(Path(record["bai"]), record["bai_identity"])
    provenance = {"schema_version": 2, "started_utc": started_utc,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "results_dir": str(args.results_dir.resolve()),
        "metadata": {"path": str(args.metadata.resolve()), "sha256": sha256(args.metadata)},
        "sample_order": list(bams), "conversion_mode": args.conversion_mode,
        "genome": args.genome, "bam2pat_defaults": BAM2PAT_DEFAULTS,
        "reference": reference, "bam_inputs": integrity, "pat_artifacts": artifacts,
        "commands": commands, "versions": versions,
        "helper_sha256": sha256(Path(__file__)),
        "legacy_reconciliation": reconciliation}
    with (args.output_dir / "pat_provenance.json").open("x") as handle:
        json.dump(provenance, handle, indent=2); handle.write("\n")
    with (args.output_dir / "pat_manifest.tsv").open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("sample", "pat", "csi", "lbeta"))
        for row in artifacts:
            writer.writerow([row[key] for key in ("sample", "pat", "csi", "lbeta")])
    (args.output_dir / "PAT_COMPLETED").write_text(provenance["completed_utc"] + "\n")
    (args.output_dir / "pat_execution.json").write_text(json.dumps({
        "started_utc": started_utc, "completed_utc": provenance["completed_utc"],
        "ended_utc": provenance["completed_utc"], "status": "completed"}, indent=2) + "\n")


def validate_cli_output(output, owner=None):
    """Standalone bundles are dated; internal stages require an incomplete owner."""
    if '..' in output.parts:
        raise LineageInputError('Output traversal is forbidden')
    output = output.absolute()
    if output.parent.resolve(strict=True) != output.parent:
        raise LineageInputError('Output parent must be canonical, not a symlink alias')
    for ancestor in output.parents:
        ownership = ancestor / '.run-ownership.json'
        if ownership.is_file() and json.loads(ownership.read_text()).get('status') == 'completed':
            raise LineageInputError('Completed output ancestors are immutable')
        if any((ancestor / marker).exists() for marker in ('PAT_COMPLETED', 'EPIDISHR_COMPLETED', 'COMPLETE', 'COHORT_COMPLETE', 'DISPLAY_COMPLETED')):
            raise LineageInputError('Completed output ancestors are immutable')
    if owner:
        owner = Path(owner).absolute()
        receipt = owner / '.run-ownership.json'
        if owner.is_symlink() or receipt.is_symlink():
            raise LineageInputError('Output owner cannot be an alias')
        record = json.loads(receipt.read_text())
        if record.get('output') != str(owner) or record.get('status') not in ('reserved', 'running', 'failed'):
            raise LineageInputError('Internal PAT stage requires an owned incomplete bundle')
        staged = any(p.name.startswith('.' + owner.name + '.running-') for p in output.parents)
        if not output.is_relative_to(owner) and not staged:
            raise LineageInputError('PAT stage does not belong to the selected output owner')
    else:
        if not re.fullmatch(r'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-[A-Za-z0-9][A-Za-z0-9._-]*', output.name):
            raise LineageInputError('Standalone PAT output requires YYYY-MM-DD_HH-MM-description (UTC)')
        datetime.strptime(output.name[:16], '%Y-%m-%d_%H-%M')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--owner-output", type=Path, default=os.environ.get("CIRCALL_OUTPUT_OWNER") or os.environ.get("EXPERIMENT_OWNER_OUTPUT"), help="Owned incomplete parent bundle for an internal prepared stage")
    parser.add_argument("--reuse-dir", type=Path, action="append", default=[])
    parser.add_argument("--conversion-mode", choices=("reuse", "rebuild", "complete"), default="reuse")
    parser.add_argument("--legacy-reuse-dir", type=Path, action="append", default=[],
                        help="Explicit original legacy bundle to reconcile; never a filename fallback")
    parser.add_argument("--reuse-sample", action="append", default=[])
    parser.add_argument("--allow-subset", action="store_true", help="Select metadata samples from a larger source run")
    parser.add_argument("--reference-manifest", type=Path)
    parser.add_argument("--preflight-only", action="store_true", help="Validate sources without hashing BAMs or creating outputs")
    parser.add_argument("--genome", default="hg38")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--wgbstools", type=Path, required=True)
    parser.add_argument("--samtools", type=Path, required=True)
    parser.add_argument("--tabix", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    args = parser.parse_args()
    if not args.preflight_only:
        # Submitted helper copies execute only while their committed source remains
        # identical. Historical snapshots retain their original entrypoint.
        import sys
        project = Path('/path/to/hpc-home/analysis/CIRCALL_methylomics')
        sys.path.insert(0, str(project / 'scripts/lib'))
        import output_guard
        output_guard.require_clean_source(project)
        if sha256(Path(__file__)) != sha256(project / 'workflows/pat/prepare_pat.py'):
            raise LineageInputError('Submitted PAT helper differs from committed project source')
        validate_cli_output(args.output_dir, args.owner_output)
    if args.threads < 1:
        parser.error("--threads must be positive")
    return args


if __name__ == "__main__":
    def interrupted(signum, frame):
        raise SystemExit(128 + signum)
    signal.signal(signal.SIGTERM, interrupted)
    try:
        prepare(parse_args())
    except (LineageInputError, OSError, ValueError, KeyError) as error:
        raise SystemExit(f"ERROR: {error}") from error
