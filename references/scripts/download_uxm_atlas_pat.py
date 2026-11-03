#!/usr/bin/env python3
"""Download and verify official UXM reference PAT files from NCBI GEO GSE186458.

Usage:
    python3 download_uxm_atlas_pat.py --output-dir /path/to/hpc-home/analysis/CIRCALL_methylomics/references/uxm_atlas_pat
"""

import argparse
import concurrent.futures
import csv
import json
import subprocess
import sys
import time
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Download and verify official UXM reference PAT files from GSE186458."
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("/path/to/hpc-home/analysis/CIRCALL_methylomics/references/uxm_atlas_pat/download_manifest.tsv"),
        help="Path to download manifest TSV file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/path/to/hpc-home/analysis/CIRCALL_methylomics/references/uxm_atlas_pat"),
        help="Target output directory for downloaded PAT and CSI files.",
    )
    parser.add_argument(
        "--concurrency",
        type=int,
        default=4,
        help="Number of concurrent download processes (default: 4).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report download actions without transferring files.",
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Verify size and presence of existing files without downloading.",
    )
    parser.add_argument(
        "--sample",
        type=str,
        default=None,
        help="Optional single sample GSM ID or name to filter.",
    )
    args = parser.parse_args()
    if args.concurrency < 1:
        parser.error("--concurrency must be positive")
    if args.dry_run and args.verify_only:
        parser.error("--dry-run and --verify-only are mutually exclusive")
    return args


def load_manifest(manifest_path, sample_filter=None):
    if not manifest_path.is_file():
        sys.exit(f"Error: Manifest file not found: {manifest_path}")

    records = []
    with manifest_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if sample_filter:
                if sample_filter not in (row["gsm_id"], row["sample_name"]):
                    continue
            records.append({
                "gsm_id": row["gsm_id"],
                "sample_name": row["sample_name"],
                "cell_type": row["cell_type"],
                "pat_filename": row["pat_filename"],
                "pat_url": row["pat_url"],
                "pat_size_bytes": int(row["pat_size_bytes"]),
                "csi_filename": row["csi_filename"],
                "csi_url": row["csi_url"],
                "csi_size_bytes": int(row["csi_size_bytes"]),
            })
    if not records:
        sys.exit("Error: No samples matched the manifest/filter")
    return records


def download_single_file(url, destination, expected_size, dry_run=False):
    """Download a file with curl, resume capability, and size verification."""
    if destination.is_file():
        actual_size = destination.stat().st_size
        if actual_size == expected_size:
            return {"status": "CACHED", "path": str(destination), "bytes": actual_size}
        return {
            "status": "SIZE_MISMATCH", "path": str(destination),
            "error": f"Existing file has {actual_size} bytes; expected {expected_size}. "
                     "Refusing to delete or overwrite; resolve explicitly before retrying.",
        }

    if dry_run:
        return {"status": "DRY_RUN", "path": str(destination), "bytes": expected_size}

    destination.parent.mkdir(parents=True, exist_ok=True)
    part_path = destination.with_suffix(destination.suffix + ".part")
    if part_path.exists() and part_path.stat().st_size > expected_size:
        return {
            "status": "SIZE_MISMATCH", "path": str(part_path),
            "error": "Partial file exceeds expected size; refusing to modify it.",
        }
    if part_path.exists() and part_path.stat().st_size == expected_size:
        part_path.rename(destination)
        return {"status": "DOWNLOADED", "path": str(destination), "bytes": expected_size}

    # Use curl with resume (--continue-at -) and retries
    curl_cmd = [
        "curl",
        "--fail",
        "--location",
        "--silent",
        "--show-error",
        "--retry", "5",
        "--retry-delay", "2",
        "--retry-max-time", "120",
        "--continue-at", "-",
        "--output", str(part_path),
        url,
    ]

    result = subprocess.run(curl_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return {
            "status": "ERROR",
            "path": str(destination),
            "error": f"curl failed ({result.returncode}): {result.stderr.strip()}",
        }

    actual_size = part_path.stat().st_size
    if actual_size != expected_size:
        return {
            "status": "SIZE_MISMATCH",
            "path": str(destination),
            "error": f"Expected {expected_size} bytes, got {actual_size} bytes",
        }

    part_path.rename(destination)
    return {"status": "DOWNLOADED", "path": str(destination), "bytes": actual_size}


def verify_file(destination, expected_size):
    if not destination.is_file():
        return {"status": "MISSING", "path": str(destination), "expected": expected_size}
    actual_size = destination.stat().st_size
    if actual_size != expected_size:
        return {
            "status": "SIZE_MISMATCH",
            "path": str(destination),
            "actual": actual_size,
            "expected": expected_size,
        }
    return {"status": "OK", "path": str(destination), "bytes": actual_size}


def main():
    args = parse_args()
    manifest_records = load_manifest(args.manifest, args.sample)
    output_dir = args.output_dir.resolve()
    if not args.dry_run and not args.verify_only:
        output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loaded {len(manifest_records)} sample records from manifest: {args.manifest}")
    print(f"Target directory: {output_dir}")

    total_pat_bytes = sum(r["pat_size_bytes"] for r in manifest_records)
    total_csi_bytes = sum(r["csi_size_bytes"] for r in manifest_records)
    total_bytes = total_pat_bytes + total_csi_bytes
    print(
        f"Total payload: {total_bytes / (1024**3):.2f} GiB "
        f"({total_bytes / (1000**3):.2f} GB) across {len(manifest_records) * 2} files"
    )

    if args.verify_only:
        print("\n--- Verifying Existing Files ---")
        missing = 0
        mismatched = 0
        verified = 0
        for r in manifest_records:
            pat_res = verify_file(output_dir / r["pat_filename"], r["pat_size_bytes"])
            csi_res = verify_file(output_dir / r["csi_filename"], r["csi_size_bytes"])
            for res in (pat_res, csi_res):
                if res["status"] == "OK":
                    verified += 1
                elif res["status"] == "MISSING":
                    missing += 1
                    print(f"MISSING: {res['path']}")
                elif res["status"] == "SIZE_MISMATCH":
                    mismatched += 1
                    print(f"SIZE MISMATCH: {res['path']} (got {res['actual']}, expected {res['expected']})")

        print(f"\nVerification summary: {verified} verified OK, {missing} missing, {mismatched} mismatched.")
        if missing > 0 or mismatched > 0:
            sys.exit(1)
        return

    # Build tasks list: each sample has 2 files: CSI (index) then PAT
    tasks = []
    for r in manifest_records:
        tasks.append((r["csi_url"], output_dir / r["csi_filename"], r["csi_size_bytes"]))
        tasks.append((r["pat_url"], output_dir / r["pat_filename"], r["pat_size_bytes"]))

    print(f"\nStarting transfer with concurrency {args.concurrency}...")
    start_time = time.time()
    results = {"DOWNLOADED": 0, "CACHED": 0, "DRY_RUN": 0, "ERROR": 0, "SIZE_MISMATCH": 0}
    errors = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.concurrency) as executor:
        futures = {
            executor.submit(
                download_single_file, url, dest, size, args.dry_run
            ): (dest.name, size)
            for url, dest, size in tasks
        }

        for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
            fname, size = futures[future]
            try:
                res = future.result()
                status = res["status"]
                results[status] = results.get(status, 0) + 1
                if status in ("ERROR", "SIZE_MISMATCH"):
                    print(f"[{i}/{len(tasks)}] FAILED: {fname} - {res.get('error')}", file=sys.stderr)
                    errors.append(res)
                elif status == "CACHED":
                    print(f"[{i}/{len(tasks)}] CACHED: {fname}")
                elif status == "DOWNLOADED":
                    print(f"[{i}/{len(tasks)}] DOWNLOADED: {fname} ({size / (1024**2):.1f} MB)")
                elif status == "DRY_RUN":
                    print(f"[{i}/{len(tasks)}] DRY_RUN: {fname} ({size / (1024**2):.1f} MB)")
            except Exception as exc:
                print(f"[{i}/{len(tasks)}] EXCEPTION: {fname} - {exc}", file=sys.stderr)
                results["ERROR"] += 1
                errors.append({"path": fname, "error": str(exc)})

    elapsed = time.time() - start_time
    print(f"\nCompleted in {elapsed:.1f} seconds.")
    print(f"Results: {results}")

    if args.dry_run:
        if errors:
            sys.exit(f"Dry-run found {len(errors)} errors; no files changed.")
        print("Dry-run complete; no files changed.")
        return

    # Write summary record
    summary_path = output_dir / "download_summary.json"
    summary_data = {
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "manifest": str(args.manifest),
        "output_dir": str(output_dir),
        "total_samples": len(manifest_records),
        "total_files": len(tasks),
        "total_bytes": total_bytes,
        "results": results,
        "errors": errors,
    }
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(summary_data, f, indent=2)

    if errors:
        sys.exit(f"Download completed with {len(errors)} errors. Check {summary_path}")
    print(f"All files match expected byte counts (content checksums not verified). Summary written to {summary_path}")


if __name__ == "__main__":
    main()
