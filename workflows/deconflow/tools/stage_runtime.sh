#!/usr/bin/env bash
#SBATCH --job-name=deconflow-stage
#SBATCH --partition=shortq
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
set -euo pipefail
workflow="${DECONFLOW_WORKFLOW:-/path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/deconflow}"
[[ -d "$workflow" && -f "$workflow/tools/stage_runtime.py" ]]
root="${1:?Usage: sbatch stage_runtime.sh NEW_DATED_RUNTIME}"
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'Submit with sbatch' >&2; exit 1; }
source "$workflow/../../scripts/lib/result_names.sh"
require_date_first_result_name "$(basename "$root")"
[[ ! -e "$root" ]] || { echo "Runtime already exists: $root" >&2; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nextflow
module load containers/apptainer/1.3.6
[[ "$(apptainer --version)" == 'apptainer version 1.3.6' ]]
export APPTAINER_CACHEDIR=/path/to/hpc-scratch/nextflow/apptainer_cache
export APPTAINER_TMPDIR="${SLURM_TMPDIR:-/tmp}/deconflow-build-${SLURM_JOB_ID}"
mkdir -p "$APPTAINER_TMPDIR"
python "$workflow/tools/stage_runtime.py" "$root"
