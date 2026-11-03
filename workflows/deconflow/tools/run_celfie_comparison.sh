#!/usr/bin/env bash
#SBATCH --job-name=uxm-celfie-comparison
#SBATCH --partition=shortq
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
set -euo pipefail
workflow="${DECONFLOW_WORKFLOW:-/path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/deconflow}"
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'Submit with sbatch' >&2; exit 1; }
[[ -d "$workflow" && -f "$workflow/tools/compare_celfie_uxm.py" ]] || { echo 'Comparison helper missing' >&2; exit 1; }
command -v conda >/dev/null || { echo 'Required command missing: conda' >&2; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nextflow
for executable in python scontrol; do command -v "$executable" >/dev/null; done
export TZ=UTC
exec python "$workflow/tools/compare_celfie_uxm.py" "$@"
