#!/bin/bash
#SBATCH --job-name=nf_core_methylseq_integrated
#SBATCH --time=7-00:00:00
#SBATCH --partition=longq
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

# Project-owned entrypoint; --submit chooses scheduler logs before sbatch.
set -euo pipefail
PROJECT_ROOT="${CIRCALL_PROJECT_ROOT:-/path/to/hpc-home/analysis/CIRCALL_methylomics}"
[[ -d "$PROJECT_ROOT" ]] || { echo "Missing project directory" >&2; exit 1; }
HELPER="$PROJECT_ROOT/scripts/methylseq_launch.py"
[[ -f "$HELPER" ]] || { echo "Missing project launcher helper" >&2; exit 1; }
command -v python3 >/dev/null || { echo "python3 is unavailable" >&2; exit 1; }
case "${1:-}" in
    --preflight|--submit)
        ACTION="${1#--}"
        shift
        exec python3 "$HELPER" "$ACTION" --project-root "$PROJECT_ROOT" "$@"
        ;;
    "") ;;
    *) echo "Use --preflight or --submit with --run-name NAME [--resume-existing]" >&2; exit 1 ;;
esac
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo "Submit full executions with this launcher's --submit option" >&2; exit 1; }
cmp -s -- "${BASH_SOURCE[0]}" "$PROJECT_ROOT/run_methylseq_integrated.sh" || {
    echo "Submitted launcher differs from committed project source; resubmit" >&2; exit 1;
}
command -v conda >/dev/null || { echo "conda is unavailable" >&2; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nextflow
command -v module >/dev/null || { echo "module is unavailable" >&2; exit 1; }
module load containers/apptainer/1.3.6
command -v nextflow >/dev/null || { echo "nextflow is unavailable" >&2; exit 1; }
command -v apptainer >/dev/null || { echo "apptainer is unavailable" >&2; exit 1; }
exec python3 "$HELPER" run --project-root "$PROJECT_ROOT"
