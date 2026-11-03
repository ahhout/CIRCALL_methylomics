#!/bin/bash
#SBATCH --job-name=emseq_uxm_heme
#SBATCH --output=/path/to/hpc-home/analysis/CIRCALL_methylomics/runs/uxm/slurm_hematopoietic_%j.out
#SBATCH --error=/path/to/hpc-home/analysis/CIRCALL_methylomics/runs/uxm/slurm_hematopoietic_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=mediumq
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail
export UXM_EXECUTION_STARTED_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

# Slurm executes a spool copy; BASH_SOURCE does not locate repository siblings.
readonly SCRIPT_DIR="/path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/uxm"
readonly SHARED_LAUNCHER="${SCRIPT_DIR}/run_lineage_deconvolution.sh"

[[ -x "${SHARED_LAUNCHER}" ]] || { echo "Missing executable: ${SHARED_LAUNCHER}" >&2; exit 1; }

export UXM_LAUNCHER_PROFILE="hematopoietic"
exec "${SHARED_LAUNCHER}"
