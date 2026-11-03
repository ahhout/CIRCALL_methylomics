#!/bin/bash
#SBATCH --job-name=uxm_atlas_pat_dl
#SBATCH --output=/path/to/hpc-home/analysis/CIRCALL_methylomics/logs/uxm_atlas_pat_dl_%j.out
#SBATCH --error=/path/to/hpc-home/analysis/CIRCALL_methylomics/logs/uxm_atlas_pat_dl_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=mediumq
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4

set -euo pipefail

readonly PROJECT_DIR="/path/to/hpc-home/analysis/CIRCALL_methylomics"
readonly SCRIPT_DIR="${PROJECT_DIR}/references/scripts"
readonly OUTPUT_DIR="${PROJECT_DIR}/references/uxm_atlas_pat"
readonly MANIFEST="${OUTPUT_DIR}/download_manifest.tsv"
readonly DOWNLOADER="${SCRIPT_DIR}/download_uxm_atlas_pat.py"
readonly LOG_DIR="${PROJECT_DIR}/logs"

[[ -d "${PROJECT_DIR}" ]] || { echo "Missing project directory: ${PROJECT_DIR}" >&2; exit 1; }
[[ -d "${SCRIPT_DIR}" ]] || { echo "Missing script directory: ${SCRIPT_DIR}" >&2; exit 1; }
[[ -s "${MANIFEST}" ]] || { echo "Missing download manifest: ${MANIFEST}" >&2; exit 1; }
[[ -x "${DOWNLOADER}" ]] || { echo "Missing executable downloader: ${DOWNLOADER}" >&2; exit 1; }

command -v python3 >/dev/null 2>&1 || { echo "python3 is unavailable" >&2; exit 1; }
command -v curl >/dev/null 2>&1 || { echo "curl is unavailable" >&2; exit 1; }

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

echo "======================================================================"
echo "Starting UXM reference PAT download (GSE186458, hg38)"
echo "Job ID:            ${SLURM_JOB_ID:-manual}"
echo "Start time:        $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
echo "Manifest:          ${MANIFEST}"
echo "Output directory:  ${OUTPUT_DIR}"
echo "Concurrency:       4"
echo "======================================================================"

python3 -u "${DOWNLOADER}" \
    --manifest "${MANIFEST}" \
    --output-dir "${OUTPUT_DIR}" \
    --concurrency 4

echo "======================================================================"
echo "Running post-download file-size verification (no content checksums)..."
echo "======================================================================"

python3 "${DOWNLOADER}" \
    --manifest "${MANIFEST}" \
    --output-dir "${OUTPUT_DIR}" \
    --verify-only

echo "======================================================================"
echo "UXM reference PAT download completed successfully."
echo "Finish time: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
echo "======================================================================"
