#!/bin/bash
#SBATCH --job-name=emseq_pat
#SBATCH --output=/path/to/hpc-home/analysis/CIRCALL_methylomics/runs/pat/pat_%j.out
#SBATCH --error=/path/to/hpc-home/analysis/CIRCALL_methylomics/runs/pat/pat_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=mediumq
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
set -euo pipefail
export PAT_EXECUTION_STARTED_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
readonly SCRIPT_DIR="/path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/pat"
[[ -d "${SCRIPT_DIR}" && -s "${SCRIPT_DIR}/launch_pat.py" ]] || { echo 'Missing PAT launcher implementation' >&2; exit 1; }
if [[ "${1:-}" == --help ]]; then
    exec python3 "${SCRIPT_DIR}/launch_pat.py" --help
fi
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'Submit with sbatch' >&2; exit 1; }
command -v conda >/dev/null 2>&1 || { echo 'conda unavailable' >&2; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nextflow
exec python3 "${SCRIPT_DIR}/launch_pat.py" "$@"
