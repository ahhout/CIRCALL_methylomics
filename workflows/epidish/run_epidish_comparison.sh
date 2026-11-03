#!/bin/bash
#SBATCH --job-name=emseq_epidish_cmp
#SBATCH --output=/path/to/hpc-home/analysis/CIRCALL_methylomics/runs/epidish/slurm_%j.out
#SBATCH --error=/path/to/hpc-home/analysis/CIRCALL_methylomics/runs/epidish/slurm_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=mediumq
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

set -euo pipefail
readonly STARTED_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

readonly PROJECT_DIR="/path/to/hpc-home/analysis/CIRCALL_methylomics"
readonly OUTPUT_CONTRACT="${PROJECT_DIR}/workflows/output_contract.py"
readonly EPIDISHR_DIR="${PROJECT_DIR}/workflows/epidish"
readonly SCRIPT_DIR="${EPIDISHR_DIR}/scripts"
readonly ATLAS="${PROJECT_DIR}/workflows/uxm/reference/Atlas.U25.l4.hg38.full.tsv"
readonly METADATA="${EMSEQ_METADATA:-${PROJECT_DIR}/metadata/metadata.tsv}"
readonly SAMPLESHEET="${EMSEQ_SAMPLESHEET:-${PROJECT_DIR}/samplesheets/samplesheet.csv}"
readonly UXM_DIR="${PROJECT_DIR}/workflows/uxm/software/UXM_deconv"
readonly EXPECTED_ATLAS_SHA256="e58fbec4fe6200a311540c6e1a16be901a2d1e9ed9d74b6a57888f2feb635f51"
readonly EXPECTED_UXM_COMMIT="8d0bb456742ece21802bc8e651e080869c67c944"
readonly WORK_BASE="/path/to/hpc-scratch/tmp/em_seq/downstream/epidishr"
readonly RESULTS_BASE="/path/to/hpc-scratch/nextflow/nf-core-methylseq/results/epidishr"
readonly REPOSITORY_RESULTS="${EPIDISHR_DIR}/results"

: "${METHYLSEQ_RESULTS_DIR:?METHYLSEQ_RESULTS_DIR must name the source shared by both UXM runs}"
: "${EPIDISHR_RUN_NAME:?EPIDISHR_RUN_NAME is required}"
: "${UXM_FULL_RUN_DIR:?UXM_FULL_RUN_DIR is required}"
: "${UXM_HEMATOPOIETIC_RUN_DIR:?UXM_HEMATOPOIETIC_RUN_DIR is required}"
readonly EPIDISHR_RPC_MAXIT="${EPIDISHR_RPC_MAXIT:-75}"
readonly MATCHED_BASELINE_DIR="${EPIDISHR_MATCHED_BASELINE_DIR:-}"
readonly EPIDISHR_RUN_NAME UXM_FULL_RUN_DIR UXM_HEMATOPOIETIC_RUN_DIR EPIDISHR_RPC_MAXIT
readonly JOB_ID="${SLURM_JOB_ID:-manual}"
readonly RECOVER_JOB_ID="${EPIDISHR_RECOVER_JOB_ID:-}"
if [[ -n "${RECOVER_JOB_ID}" ]]; then
    [[ "${RECOVER_JOB_ID}" =~ ^[0-9]+$ ]] || { echo "EPIDISHR_RECOVER_JOB_ID must be numeric" >&2; exit 1; }
    readonly WORK_DIR="${WORK_BASE}/.${EPIDISHR_RUN_NAME}.running-${RECOVER_JOB_ID}"
else
    readonly WORK_DIR="${WORK_BASE}/.${EPIDISHR_RUN_NAME}.running-${JOB_ID}"
fi
readonly FINAL_DIR="${RESULTS_BASE}/${EPIDISHR_RUN_NAME}"
readonly MATCHED_UXM_TMP="${WORK_DIR}/.matched_uxm_tmp"
readonly RUN_LOG="${PROJECT_DIR}/runs/epidish/${EPIDISHR_RUN_NAME}_${JOB_ID}.log"

source /path/to/hpc-home/analysis/CIRCALL_methylomics/scripts/lib/result_names.sh
if [[ -n "${RECOVER_JOB_ID}" && -d "${WORK_DIR}" ]]; then
    require_date_first_result_name "${EPIDISHR_RUN_NAME}" legacy-existing
else
    require_date_first_result_name "${EPIDISHR_RUN_NAME}"
fi
[[ "${EPIDISHR_RUN_NAME}" =~ ^[A-Za-z0-9._-]+$ ]] || {
    echo "Unsafe EPIDISHR_RUN_NAME: ${EPIDISHR_RUN_NAME}" >&2; exit 1;
}
[[ "${EPIDISHR_RPC_MAXIT}" =~ ^[1-9][0-9]*$ ]] || {
    echo "EPIDISHR_RPC_MAXIT must be a positive integer" >&2; exit 1;
}
[[ -s "${ATLAS}" ]] || { echo "Missing canonical atlas: ${ATLAS}" >&2; exit 1; }
[[ -s "${METADATA}" ]] || { echo "Missing metadata: ${METADATA}" >&2; exit 1; }
[[ -d "${UXM_FULL_RUN_DIR}" ]] || { echo "Missing full UXM run: ${UXM_FULL_RUN_DIR}" >&2; exit 1; }
[[ -d "${UXM_HEMATOPOIETIC_RUN_DIR}" ]] || {
    echo "Missing haematopoietic UXM run: ${UXM_HEMATOPOIETIC_RUN_DIR}" >&2; exit 1;
}
[[ -d "${UXM_DIR}/.git" ]] || { echo "Missing pinned UXM checkout: ${UXM_DIR}" >&2; exit 1; }
for script in prepare_epidish.R finalize_comparisons.R reuse_matched_uxm.R \
    run_unweighted_nnls.py recovery_contract.py epidish_common.R plot_epidish_standalone.R; do
    [[ -s "${SCRIPT_DIR}/${script}" ]] || { echo "Missing workflow script: ${SCRIPT_DIR}/${script}" >&2; exit 1; }
done
if [[ -n "${MATCHED_BASELINE_DIR}" ]]; then
    [[ -s "${MATCHED_BASELINE_DIR}/EPIDISHR_COMPLETED" ]] || {
        echo "Matched-UXM baseline is not a completed result: ${MATCHED_BASELINE_DIR}" >&2; exit 1;
    }
fi
[[ -L "${REPOSITORY_RESULTS}" ]] || { echo "Missing repository results symlink: ${REPOSITORY_RESULTS}" >&2; exit 1; }
[[ "$(readlink -f "${REPOSITORY_RESULTS}")" == "${RESULTS_BASE}" ]] || {
    echo "Repository results symlink does not resolve to ${RESULTS_BASE}" >&2; exit 1;
}
if [[ -z "${RECOVER_JOB_ID}" ]]; then
    [[ ! -e "${FINAL_DIR}" && ! -L "${FINAL_DIR}" ]] || { echo "Final output already exists: ${FINAL_DIR}" >&2; exit 1; }
fi
if [[ -n "${RECOVER_JOB_ID}" ]]; then
    [[ -d "${WORK_DIR}" ]] || { echo "Recovery work directory is missing: ${WORK_DIR}" >&2; exit 1; }
    [[ -s "${WORK_DIR}/provenance/recovery_contract.json" ]] || {
        echo "Missing recovery contract; legacy incomplete work requires a fresh run directory" >&2; exit 1;
    }
    [[ ! -e "${WORK_DIR}/EPIDISHR_COMPLETED" ]] || {
        echo "Recovery work directory is already complete: ${WORK_DIR}" >&2; exit 1;
    }
else
    [[ ! -e "${WORK_DIR}" && ! -L "${WORK_DIR}" ]] || { echo "Work directory already exists: ${WORK_DIR}" >&2; exit 1; }
fi

command -v conda >/dev/null 2>&1 || { echo "conda is unavailable" >&2; exit 1; }
command -v sha256sum >/dev/null 2>&1 || { echo "sha256sum is unavailable" >&2; exit 1; }
echo "${EXPECTED_ATLAS_SHA256}  ${ATLAS}" | sha256sum --check -
actual_uxm_commit="$(git -C "${UXM_DIR}" rev-parse HEAD)"
[[ "${actual_uxm_commit}" == "${EXPECTED_UXM_COMMIT}" ]] || {
    echo "Unexpected UXM commit: ${actual_uxm_commit}" >&2; exit 1;
}

# Validate source equality and shared PAT payloads before creating comparison work.
conda run --no-capture-output -n nextflow python3 \
    "${PROJECT_DIR}/workflows/pat/validate_pat_manifest.py" \
    --results-dir "${METHYLSEQ_RESULTS_DIR}" --pat-dir "${UXM_FULL_RUN_DIR}" \
    --compare-with "${UXM_HEMATOPOIETIC_RUN_DIR}" \
    --samplesheet "${SAMPLESHEET}"

# A permanent lock inode beside the final destination serializes every job ID.
command -v flock >/dev/null 2>&1 || { echo "flock is unavailable" >&2; exit 1; }
mkdir -p "${WORK_BASE}" "${RESULTS_BASE}"
exec {publication_lock}>"${FINAL_DIR}.lock"
flock -n "${publication_lock}" || { echo "Run destination is locked: ${FINAL_DIR}" >&2; exit 1; }
if [[ -z "${RECOVER_JOB_ID}" ]]; then
    [[ ! -e "${FINAL_DIR}" && ! -L "${FINAL_DIR}" ]] || { echo "Final output already exists: ${FINAL_DIR}" >&2; exit 1; }
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate epidishr
for executable in Rscript sha256sum; do
    command -v "${executable}" >/dev/null 2>&1 || { echo "Missing executable in epidishr: ${executable}" >&2; exit 1; }
done
Rscript -e 'stopifnot(as.character(packageVersion("EpiDISH")) == "2.26.0")'
readonly R_IDENTITY="$(Rscript -e 'cat(as.character(getRversion()), as.character(packageVersion("EpiDISH")))')"
conda activate nextflow
readonly CONTRACT_PYTHON="$(command -v python3)"
conda activate epidishr
ownership_args=("${FINAL_DIR}" --workflow "${EPIDISHR_DIR}" --input "${METADATA}" --input "${SAMPLESHEET}" --work-root "${WORK_DIR}")
if [[ -n "${RECOVER_JOB_ID}" ]]; then
    "${CONTRACT_PYTHON}" "${OUTPUT_CONTRACT}" resume "${ownership_args[@]}"
else
    "${CONTRACT_PYTHON}" "${OUTPUT_CONTRACT}" reserve "${ownership_args[@]}"
fi
contract_args=(reserve --work-dir "${WORK_DIR}" --atlas "${ATLAS}" --metadata "${METADATA}"
    --samplesheet "${SAMPLESHEET}" --full-run "${UXM_FULL_RUN_DIR}"
    --heme-run "${UXM_HEMATOPOIETIC_RUN_DIR}" --uxm-checkout "${UXM_DIR}"
    --source-results "${METHYLSEQ_RESULTS_DIR}" --r-identity "${R_IDENTITY}"
    --rpc-maxit "${EPIDISHR_RPC_MAXIT}" --baseline "${MATCHED_BASELINE_DIR}"
    --started-utc "${STARTED_UTC}" --attempt-id "${JOB_ID}")
if [[ -n "${RECOVER_JOB_ID}" ]]; then contract_args+=(--recover); fi
# This freezes new inputs, or validates all recovery evidence, before preparation writes.
"${CONTRACT_PYTHON}" "${SCRIPT_DIR}/recovery_contract.py" "${contract_args[@]}"
attempt_root="${WORK_DIR}"
attempt_finished=false
finish_attempt() {
    status=$?
    trap - EXIT
    if [[ "${attempt_finished}" != true ]]; then
        "${CONTRACT_PYTHON}" "${OUTPUT_CONTRACT}" fail "${FINAL_DIR}" || true
        "${CONTRACT_PYTHON}" "${SCRIPT_DIR}/recovery_contract.py" finish \
            --work-dir "${attempt_root}" --attempt-id "${JOB_ID}" --status failed || true
    fi
    exit "${status}"
}
trap finish_attempt EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
mkdir -p "${MATCHED_UXM_TMP}"
# Completed finalization seals commands.log through the finalization receipt.
if [[ -s "${WORK_DIR}/provenance/receipts/finalization.json" ]]; then
    readonly COMMANDS_LOG="${WORK_DIR}/provenance/attempt_commands/${JOB_ID}.log"
    mkdir -p "$(dirname "${COMMANDS_LOG}")"
else
    readonly COMMANDS_LOG="${WORK_DIR}/provenance/commands.log"
fi

mkdir -p "${PROJECT_DIR}/runs/epidish"
exec > >(tee -a "${RUN_LOG}") 2>&1
echo "EpiDISH-UXM comparison started: $(date --iso-8601=seconds)"
echo "Run name: ${EPIDISHR_RUN_NAME}"
echo "EpiDISH RPC maximum iterations: ${EPIDISHR_RPC_MAXIT}"
if [[ -n "${RECOVER_JOB_ID}" ]]; then
    echo "Recovering incomplete work from job: ${RECOVER_JOB_ID}"
fi
echo "Full UXM source: ${UXM_FULL_RUN_DIR}"
echo "Haematopoietic UXM source: ${UXM_HEMATOPOIETIC_RUN_DIR}"

if [[ -n "${RECOVER_JOB_ID}" ]]; then
    printf '# Recovery job %q at %q\n' "${JOB_ID}" "$(date --iso-8601=seconds)" >> "${COMMANDS_LOG}"
else
    : > "${COMMANDS_LOG}"
fi
printf 'Rscript %q %q %q %q %q %q %q %q\n' \
    "${SCRIPT_DIR}/prepare_epidish.R" "${ATLAS}" "${UXM_FULL_RUN_DIR}" \
    "${UXM_HEMATOPOIETIC_RUN_DIR}" "${METADATA}" "${WORK_DIR}" \
    "${UXM_DIR}" "${EPIDISHR_RPC_MAXIT}" >> "${COMMANDS_LOG}"
Rscript "${SCRIPT_DIR}/prepare_epidish.R" \
    "${ATLAS}" \
    "${UXM_FULL_RUN_DIR}" \
    "${UXM_HEMATOPOIETIC_RUN_DIR}" \
    "${METADATA}" \
    "${WORK_DIR}" \
    "${UXM_DIR}" \
    "${EPIDISHR_RPC_MAXIT}"

if [[ -n "${MATCHED_BASELINE_DIR}" ]]; then
    printf 'Rscript %q %q %q\n' \
        "${SCRIPT_DIR}/reuse_matched_uxm.R" "${WORK_DIR}" \
        "${MATCHED_BASELINE_DIR}" >> "${COMMANDS_LOG}"
    Rscript "${SCRIPT_DIR}/reuse_matched_uxm.R" \
        "${WORK_DIR}" "${MATCHED_BASELINE_DIR}"
fi

mapfile -t samples < <(awk -F'\t' 'NR > 1 {print $1}' "${WORK_DIR}/provenance/pat_manifest.tsv")
mapfile -t pat_relpaths < <(awk -F'\t' 'NR > 1 {print $2}' "${WORK_DIR}/provenance/pat_manifest.tsv")
[[ ${#samples[@]} -gt 0 && ${#samples[@]} -eq ${#pat_relpaths[@]} ]] || {
    echo "Prepared PAT manifest is invalid" >&2; exit 1;
}
pats=()
for i in "${!samples[@]}"; do
    [[ "${samples[$i]}" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "Unsafe sample identifier" >&2; exit 1; }
    [[ "${pat_relpaths[$i]}" =~ ^pat/[A-Za-z0-9._-]+[.]pat[.]gz$ ]] || {
        echo "Unsafe PAT path: ${pat_relpaths[$i]}" >&2; exit 1;
    }
    pat="${UXM_FULL_RUN_DIR}/${pat_relpaths[$i]}"
    [[ -s "${pat}" && -s "${pat}.csi" ]] || { echo "Missing PAT/CSI: ${pat}" >&2; exit 1; }
    pats+=("${pat}")
done

conda activate nextflow
export PATH="${UXM_DIR}:${PATH}"
for executable in python3 wgbstools uxm; do
    command -v "${executable}" >/dev/null 2>&1 || { echo "Missing executable in nextflow: ${executable}" >&2; exit 1; }
done
python3 -c "import numpy, pandas, scipy" || { echo "UXM Python dependencies are unavailable" >&2; exit 1; }
printf 'python3 %q %q\n' \
    "${SCRIPT_DIR}/run_unweighted_nnls.py" "${WORK_DIR}" >> "${COMMANDS_LOG}"
python3 "${SCRIPT_DIR}/run_unweighted_nnls.py" "${WORK_DIR}"

while IFS=$'\t' read -r profile threshold label marker_count comparison_dir matched_atlas native_source; do
    [[ "${profile}" == "full" || "${profile}" == "hematopoietic" ]] || {
        echo "Invalid comparison profile: ${profile}" >&2; exit 1;
    }
    [[ "${threshold}" =~ ^(1|3|5)$ ]] || { echo "Invalid threshold: ${threshold}" >&2; exit 1; }
    [[ "${label}" == "coverage-ge${threshold}" ]] || { echo "Invalid threshold label" >&2; exit 1; }
    [[ "${marker_count}" =~ ^[0-9]+$ && "${marker_count}" -gt 0 ]] || {
        echo "Invalid marker count" >&2; exit 1;
    }
    [[ "${comparison_dir}" == "comparisons/${profile}/${label}" ]] || {
        echo "Unexpected comparison directory: ${comparison_dir}" >&2; exit 1;
    }
    [[ "${matched_atlas}" == "${comparison_dir}/uxm_matched_atlas.l4.tsv" ]] || {
        echo "Unexpected matched atlas path: ${matched_atlas}" >&2; exit 1;
    }
    comparison_dir="${WORK_DIR}/${comparison_dir}"
    matched_atlas="${WORK_DIR}/${matched_atlas}"
    [[ -s "${matched_atlas}" && -s "${native_source}" ]] || {
        echo "Missing comparison input for ${profile}/${label}" >&2; exit 1;
    }
    tmp_dir="${MATCHED_UXM_TMP}/${profile}/${label}"
    mkdir -p "${tmp_dir}"
    raw_output="${comparison_dir}/.uxm_weighted_matched.raw.csv"
    stage="matched/${profile}/${label}"
    raw_relative="comparisons/${profile}/${label}/.uxm_weighted_matched.raw.csv"
    receipt_status=0
    "${CONTRACT_PYTHON}" "${SCRIPT_DIR}/recovery_contract.py" check \
        "${WORK_DIR}" "${stage}" "${raw_relative}" || receipt_status=$?
    if [[ "${receipt_status}" -eq 0 ]]; then
        printf '# Reused completed matched UXM output: %q\n' "${raw_output}" >> "${COMMANDS_LOG}"
        echo "Reusing matched UXM output: ${profile}/${label}"
        continue
    fi
    [[ "${receipt_status}" -eq 3 ]] || { echo "Matched UXM receipt validation failed" >&2; exit 1; }
    printf 'uxm deconv --atlas %q --tmp_dir %q --threads %q' \
        "${matched_atlas}" "${tmp_dir}" "${SLURM_CPUS_PER_TASK:-8}" >> "${COMMANDS_LOG}"
    printf ' %q' "${pats[@]}" >> "${COMMANDS_LOG}"
    printf ' --output %q\n' "${raw_output}" >> "${COMMANDS_LOG}"
    uxm deconv \
        --atlas "${matched_atlas}" \
        --tmp_dir "${tmp_dir}" \
        --threads "${SLURM_CPUS_PER_TASK:-8}" \
        --verbose \
        "${pats[@]}" \
        --output "${raw_output}"
    [[ -s "${raw_output}" ]] || { echo "Matched UXM output is missing" >&2; exit 1; }
    "${CONTRACT_PYTHON}" "${SCRIPT_DIR}/recovery_contract.py" receipt \
        "${WORK_DIR}" "${stage}" "${raw_relative}"
done < <(tail -n +2 "${WORK_DIR}/provenance/comparison_jobs.tsv")

conda activate epidishr
printf 'Rscript %q %q %q\n' \
    "${SCRIPT_DIR}/finalize_comparisons.R" "${WORK_DIR}" "${EPIDISHR_RUN_NAME}" >> "${COMMANDS_LOG}"
Rscript "${SCRIPT_DIR}/finalize_comparisons.R" "${WORK_DIR}" "${EPIDISHR_RUN_NAME}"

[[ -s "${WORK_DIR}/EPIDISHR_COMPLETED" ]] || { echo "Completion marker is missing" >&2; exit 1; }
"${CONTRACT_PYTHON}" "${OUTPUT_CONTRACT}" publish "${FINAL_DIR}" --source "${WORK_DIR}"
attempt_root="${FINAL_DIR}"
cp "${FINAL_DIR}/provenance/run_manifest.json" \
    "${PROJECT_DIR}/runs/epidish/${EPIDISHR_RUN_NAME}_${JOB_ID}.manifest.json"
"${CONTRACT_PYTHON}" "${SCRIPT_DIR}/recovery_contract.py" finish \
    --work-dir "${FINAL_DIR}" --attempt-id "${JOB_ID}" --status completed
"${CONTRACT_PYTHON}" "${OUTPUT_CONTRACT}" complete "${FINAL_DIR}"
attempt_finished=true

echo "EpiDISH-UXM comparison completed: $(date --iso-8601=seconds)"
echo "Results: ${FINAL_DIR}"
