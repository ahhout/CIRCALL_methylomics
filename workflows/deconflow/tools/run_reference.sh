#!/usr/bin/env bash
#SBATCH --job-name=deconflow-reference
#SBATCH --partition=longq
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=72:00:00
set -euo pipefail
workflow="${DECONFLOW_WORKFLOW:-/path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/deconflow}"
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'Submit with sbatch' >&2; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nextflow
[[ "$(wgbstools --version)" == 'wgbstools version 0.3.0' ]]
python "$workflow/tools/reference.py" build "$@"
