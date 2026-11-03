#!/bin/bash
#SBATCH --job-name=emseq_uxm_full
#SBATCH --output=/path/to/hpc-home/analysis/CIRCALL_methylomics/runs/uxm/slurm_full_%j.out
#SBATCH --error=/path/to/hpc-home/analysis/CIRCALL_methylomics/runs/uxm/slurm_full_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=mediumq
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail
readonly STARTED_UTC="${UXM_EXECUTION_STARTED_UTC:-$(date -u +%Y-%m-%dT%H:%M:%SZ)}"

readonly PROJECT_DIR="/path/to/hpc-home/analysis/CIRCALL_methylomics"
readonly OUTPUT_CONTRACT="${PROJECT_DIR}/workflows/output_contract.py"
readonly LINEAGE_DIR="${PROJECT_DIR}/workflows/uxm"
# Allow a run-specific cohort snapshot while preserving the default metadata.
readonly METADATA="${EMSEQ_METADATA:-${PROJECT_DIR}/metadata/metadata.tsv}"
readonly PREPARE_PAT="${PROJECT_DIR}/workflows/pat/prepare_pat.py"
readonly UXM_DIR="${LINEAGE_DIR}/software/UXM_deconv"
readonly ATLAS="${LINEAGE_DIR}/reference/Atlas.U25.l4.hg38.full.tsv"
readonly EXPECTED_UXM_COMMIT="8d0bb456742ece21802bc8e651e080869c67c944"
readonly EXPECTED_ATLAS_SHA256="e58fbec4fe6200a311540c6e1a16be901a2d1e9ed9d74b6a57888f2feb635f51"
: "${METHYLSEQ_RESULTS_DIR:?METHYLSEQ_RESULTS_DIR is required}"
: "${UXM_RUN_BASE_NAME:?UXM_RUN_BASE_NAME is required}"
readonly METHYLSEQ_RESULTS_DIR
readonly UXM_ANALYSIS_PROFILE="${UXM_LAUNCHER_PROFILE:-full-atlas}"
readonly UXM_RUN_BASE_NAME
readonly LINEAGE_RUN_NAME="${UXM_RUN_BASE_NAME}-uxm-${UXM_ANALYSIS_PROFILE}"
: "${PAT_INPUT_DIR:?PAT_INPUT_DIR must name an explicitly validated PAT bundle}"
readonly PAT_INPUT_DIR
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo "Submit with sbatch" >&2; exit 1; }
readonly WORK_BASE="/path/to/hpc-scratch/tmp/em_seq/downstream/lineage_deconvolution"
readonly RESULTS_BASE="/path/to/hpc-scratch/nextflow/nf-core-methylseq/results/uxm_deconvolution"
readonly REPOSITORY_RESULTS="${LINEAGE_DIR}/results"
readonly JOB_ID="${SLURM_JOB_ID:-manual}"
readonly FINAL_DIR="${RESULTS_BASE}/${LINEAGE_RUN_NAME}"
readonly WORK_DIR="${WORK_BASE}/.${LINEAGE_RUN_NAME}.running-${JOB_ID}"
readonly RESULT_DIR="${WORK_DIR}/results"
readonly INTERMEDIATE_DIR="${WORK_DIR}/intermediates"
readonly RESULT_CSV="${RESULT_DIR}/lineage_fractions.csv"
readonly RESULT_PDF="${RESULT_DIR}/lineage_fractions.pdf"
readonly RUN_LOG="${PROJECT_DIR}/runs/uxm/${LINEAGE_RUN_NAME}_${JOB_ID}.log"

case "${UXM_ANALYSIS_PROFILE}" in
    full-atlas)
        declare -ar UXM_DECONV_ARGS=()
        ;;
    hematopoietic)
        declare -ar UXM_DECONV_ARGS=(
            --include
            Blood-B
            Blood-T
            Blood-NK
            Blood-Mono+Macro
            Blood-Granul
            Eryth-prog
            Megakaryocytes
        )
        ;;
    *)
        echo "UXM_ANALYSIS_PROFILE must be full-atlas or hematopoietic" >&2
        exit 1
        ;;
esac

[[ "${UXM_RUN_BASE_NAME}" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "Unsafe UXM run base name" >&2; exit 1; }
source /path/to/hpc-home/analysis/CIRCALL_methylomics/scripts/lib/result_names.sh
require_date_first_result_name "${LINEAGE_RUN_NAME}"
[[ "${LINEAGE_RUN_NAME}" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "Unsafe lineage run name" >&2; exit 1; }
[[ "${CONVERSION_MODE:-reuse}" == "reuse" ]] || {
    echo "Build one shared PAT bundle with run_prepare_pat.sh before UXM" >&2; exit 1;
}
[[ -s "${PAT_INPUT_DIR}/pat_provenance.json" ]] || { echo "Missing explicit PAT provenance" >&2; exit 1; }
[[ -x "${PREPARE_PAT}" ]] || { echo "Missing executable: ${PREPARE_PAT}" >&2; exit 1; }
[[ -d "${UXM_DIR}/.git" ]] || { echo "Missing pinned UXM checkout: ${UXM_DIR}" >&2; exit 1; }
[[ -s "${ATLAS}" ]] || { echo "Missing UXM atlas: ${ATLAS}" >&2; exit 1; }
[[ -s "${METADATA}" ]] || { echo "Missing metadata: ${METADATA}" >&2; exit 1; }
[[ -d "${METHYLSEQ_RESULTS_DIR}" ]] || { echo "Missing methylseq result: ${METHYLSEQ_RESULTS_DIR}" >&2; exit 1; }
[[ -L "${REPOSITORY_RESULTS}" ]] || { echo "Missing repository results symlink: ${REPOSITORY_RESULTS}" >&2; exit 1; }
[[ "$(readlink -f "${REPOSITORY_RESULTS}")" == "${RESULTS_BASE}" ]] || {
    echo "Repository results symlink does not resolve to ${RESULTS_BASE}" >&2; exit 1;
}
[[ ! -e "${FINAL_DIR}" && ! -L "${FINAL_DIR}" ]] || { echo "Lineage output already exists: ${FINAL_DIR}" >&2; exit 1; }
[[ ! -e "${WORK_DIR}" && ! -L "${WORK_DIR}" ]] || { echo "Lineage work directory already exists: ${WORK_DIR}" >&2; exit 1; }
for executable in flock python3; do
    command -v "${executable}" >/dev/null 2>&1 || { echo "Missing executable: ${executable}" >&2; exit 1; }
done
mkdir -p "${RESULTS_BASE}"
# Persistent lock beside the final destination, shared by every job for this name.
exec 9>"${FINAL_DIR}.lock"
flock -n 9 || { echo "UXM destination is locked: ${FINAL_DIR}" >&2; exit 1; }
[[ ! -e "${FINAL_DIR}" && ! -L "${FINAL_DIR}" ]] || { echo "Lineage output already exists: ${FINAL_DIR}" >&2; exit 1; }
[[ ! -e "${WORK_DIR}" && ! -L "${WORK_DIR}" ]] || { echo "Lineage work already exists: ${WORK_DIR}" >&2; exit 1; }
python3 "${OUTPUT_CONTRACT}" reserve "${FINAL_DIR}" --workflow "${LINEAGE_DIR}" --input "${METADATA}" --input "${PAT_INPUT_DIR}/pat_provenance.json" --work-root "${WORK_DIR}"
mkdir -p "${WORK_BASE}" "${PROJECT_DIR}/runs/uxm"
mkdir "${WORK_DIR}"
mkdir "${RESULT_DIR}" "${INTERMEDIATE_DIR}"
execution_dir="${WORK_DIR}"
record_execution() {
    python3 - "${execution_dir}/execution.json" "${STARTED_UTC}" "$1" "${JOB_ID}" <<'PY_EXECUTION'
import json, sys
from datetime import datetime, timezone
from pathlib import Path
path, started, status, job = sys.argv[1:]
record = {'started_utc': started, 'status': status, 'slurm_job_id': job}
if status != 'running':
    record['ended_utc'] = datetime.now(timezone.utc).isoformat()
if status == 'completed':
    record['completed_utc'] = record['ended_utc']
Path(path).write_text(json.dumps(record, indent=2) + '\n')
PY_EXECUTION
}
record_execution running
record_uxm_exit() {
    local status=$?
    trap - EXIT
    if [[ "$status" -ne 0 ]]; then
        python3 "${OUTPUT_CONTRACT}" fail "${FINAL_DIR}" || true
        record_execution failed_or_interrupted || echo "Failed to record UXM failure timing" >&2
    fi
    exit "$status"
}
trap record_uxm_exit EXIT
trap 'exit 143' TERM
trap 'exit 130' INT

command -v conda >/dev/null 2>&1 || { echo "conda is unavailable" >&2; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nextflow
export PATH="${UXM_DIR}:${PATH}"
export MPLBACKEND="Agg"
export MPLCONFIGDIR="${INTERMEDIATE_DIR}/matplotlib"
mkdir -p "${MPLCONFIGDIR}"
for executable in python3 git sha256sum samtools tabix wgbstools uxm; do
    command -v "${executable}" >/dev/null 2>&1 || { echo "Missing executable: ${executable}" >&2; exit 1; }
done
python3 -c "import pandas, numpy, scipy, matplotlib" || { echo "UXM Python dependencies are unavailable" >&2; exit 1; }

actual_commit="$(git -C "${UXM_DIR}" rev-parse HEAD)"
[[ "${actual_commit}" == "${EXPECTED_UXM_COMMIT}" ]] || { echo "Unexpected UXM commit: ${actual_commit}" >&2; exit 1; }
echo "${EXPECTED_ATLAS_SHA256}  ${ATLAS}" | sha256sum --check -

exec > >(tee -a "${RUN_LOG}") 2>&1
echo "Lineage workflow started: $(date --iso-8601=seconds)"
echo "Analysis profile: ${UXM_ANALYSIS_PROFILE}"
echo "Run name: ${LINEAGE_RUN_NAME}"
python3 "${PREPARE_PAT}" \
    --results-dir "${METHYLSEQ_RESULTS_DIR}" \
    --metadata "${METADATA}" \
    --output-dir "${WORK_DIR}/prepared" \
    --owner-output "${FINAL_DIR}" \
    --reuse-dir "${PAT_INPUT_DIR}" \
    --conversion-mode reuse \
    --genome hg38 \
    --threads "${SLURM_CPUS_PER_TASK:-8}" \
    --wgbstools "$(command -v wgbstools)" \
    --samtools "$(command -v samtools)" \
    --tabix "$(command -v tabix)" \
    --reference-root "$(dirname "$(dirname "$(command -v wgbstools)")")/share/wgbs_tools"

for artifact in pat pat_manifest.tsv pat_provenance.json pat_execution.json PAT_COMPLETED PAT_STARTED; do
    mv -T "${WORK_DIR}/prepared/${artifact}" "${WORK_DIR}/${artifact}"
done

mapfile -t samples < <(tail -n +2 "${WORK_DIR}/pat_manifest.tsv" | cut -f1)
mapfile -t pat_relpaths < <(tail -n +2 "${WORK_DIR}/pat_manifest.tsv" | cut -f2)
[[ ${#samples[@]} -gt 0 && ${#samples[@]} -eq ${#pat_relpaths[@]} ]] || { echo "Invalid PAT manifest" >&2; exit 1; }
pats=()
for relpath in "${pat_relpaths[@]}"; do
    [[ "${relpath}" =~ ^pat/[A-Za-z0-9._-]+[.]pat[.]gz$ ]] || { echo "Unsafe PAT manifest path: ${relpath}" >&2; exit 1; }
    pats+=("${WORK_DIR}/${relpath}")
done
for pat in "${pats[@]}"; do
    [[ -s "${pat}" && -s "${pat}.csi" ]] || { echo "Missing PAT or CSI: ${pat}" >&2; exit 1; }
done

echo "UXM commit: ${actual_commit}"
echo "UXM: $(uxm --version)"
echo "wgbstools: $(wgbstools --version)"
printf 'PAT input: %s\n' "${pats[@]}"

uxm deconv \
    --atlas "${ATLAS}" \
    --tmp_dir "${INTERMEDIATE_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK:-8}" \
    --verbose \
    "${pats[@]}" \
    "${UXM_DECONV_ARGS[@]}" \
    --output "${RESULT_CSV}"
uxm plot "${RESULT_CSV}" --outpath "${RESULT_PDF}" --min_rate 0

expected="$(IFS=,; echo "${samples[*]}")"
awk -F',' -v expected="${expected}" '
    NR == 1 {
        observed = ""; for (i = 2; i <= NF; i++) observed = observed (i == 2 ? "" : ",") $i
        if (observed != expected) { print "Unexpected sample columns: " observed > "/dev/stderr"; exit 1 }
        n = NF; next
    }
    {
        for (i = 2; i <= NF; i++) {
            if ($i !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ || $i < 0) {
                print "Invalid lineage fraction at line " NR > "/dev/stderr"; exit 1
            }
            total[i] += $i
        }
    }
    END { for (i = 2; i <= n; i++) if (total[i] < 0.999 || total[i] > 1.001) exit 1 }
' "${RESULT_CSV}"
[[ -s "${RESULT_PDF}" ]] || { echo "UXM plot was not created" >&2; exit 1; }

python3 "${OUTPUT_CONTRACT}" publish "${FINAL_DIR}" --source "${WORK_DIR}"
execution_dir="${FINAL_DIR}"
cp "${FINAL_DIR}/pat_provenance.json" "${PROJECT_DIR}/runs/uxm/${LINEAGE_RUN_NAME}_${JOB_ID}.pat.json"
record_execution completed
python3 "${OUTPUT_CONTRACT}" complete "${FINAL_DIR}"
echo "Lineage workflow completed: $(date --iso-8601=seconds)"
echo "Results: ${FINAL_DIR}/results"
