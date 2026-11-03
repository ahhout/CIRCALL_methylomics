#!/usr/bin/env bash
#SBATCH --job-name=deconflow-fixtures
#SBATCH --partition=shortq
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail
workflow="${DECONFLOW_WORKFLOW:-/path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/deconflow}"
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'Submit with sbatch' >&2; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nextflow
module load containers/apptainer/1.3.6
[[ "$(apptainer --version)" == 'apptainer version 1.3.6' ]]
exec 8>"/path/to/hpc-scratch/tmp/em_seq/.methylseq-cache-maintenance.lock"
flock -n 8 || { echo "Shared Nextflow work cache is active" >&2; exit 1; }
python "$workflow/tools/validate_runtime.py" "$@"
