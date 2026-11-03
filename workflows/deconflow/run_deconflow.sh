#!/usr/bin/env bash
#SBATCH --job-name=deconflow
#SBATCH --partition=longq
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=96:00:00
set -euo pipefail
workflow="${DECONFLOW_WORKFLOW:-/path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/deconflow}"
method=uxm
if [[ "${1:-}" == celfie ]]; then
    method=celfie
    shift
fi
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'Submit with sbatch' >&2; exit 1; }
[[ -d "$workflow" ]] || { echo "Workflow directory missing: $workflow" >&2; exit 1; }
for required in analysis.py common.py native.py fragments.py pat_inputs.py fragments.nf config/fragments.config config/flamingo.config config/params.uxm.json; do
    [[ -f "$workflow/$required" ]] || { echo "Required file missing: $workflow/$required" >&2; exit 1; }
done
if [[ "$method" == celfie ]]; then
    for required in tools/celfie_patients.py tools/celfie_patients.nf tools/celfie_diagnostics.py; do
        [[ -f "$workflow/$required" ]] || { echo "Required CELFIE file missing: $workflow/$required" >&2; exit 1; }
    done
fi
command -v conda >/dev/null || { echo 'Required command missing: conda' >&2; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nextflow
module load containers/apptainer/1.3.6
[[ "$(apptainer --version)" == 'apptainer version 1.3.6' ]]
for executable in python nextflow apptainer sbatch; do command -v "$executable" >/dev/null; done
export TZ=UTC
exec 8>"/path/to/hpc-scratch/tmp/em_seq/.methylseq-cache-maintenance.lock"
flock -n 8 || { echo "Shared Nextflow work cache is active" >&2; exit 1; }
if [[ "$method" == celfie ]]; then
    exec python "$workflow/tools/celfie_patients.py" "$@"
fi
exec python "$workflow/analysis.py" "$@"
