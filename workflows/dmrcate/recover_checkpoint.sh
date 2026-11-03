#!/usr/bin/env bash
set -euo pipefail
workflow="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ "${DMRCATE_RECOVERY_WORKER:-}" == 1 ]]; then
    # Slurm executes a spool copy; BASH_SOURCE does not identify our snapshot.
    workflow="${4:?worker attempt required}/code"
fi
project=/path/to/hpc-home/analysis/CIRCALL_methylomics
source "$project/scripts/lib/result_names.sh"
contract="$project/workflows/output_contract.py"
run="$(realpath "${1:?usage: recover_checkpoint.sh RUN CONTRAST VALIDATION}")"
key="${2:?contrast required}"
validation="$(realpath "${3:?validation receipt required}")"
[[ "$key" =~ ^[A-Za-z0-9]+_minus_[A-Za-z0-9]+$ ]]
require_date_first_result_name "$(basename "$run")"
for command in python sbatch sha256sum flock; do command -v "$command" >/dev/null; done
[[ -f "$run/config.json" && -f "$run/started_utc" && ! -e "$run/COMPLETE" ]]
[[ -f "$run/comparisons/$key/cpgannotated.rds" && ! -e "$run/comparisons/$key/INFERENCE_READY.rds" ]]
for required in code/analysis.R code/prepare.py code.sha256 input_fingerprint.json prepared/count_directory.rds prepared/samples.tsv "comparisons/$key/eligible_rows.rds"; do
    [[ -r "$run/$required" ]] || { echo "Missing recovery input: $run/$required" >&2; exit 1; }
done
for result in dmrcate.rds regions.tsv regions.zero_based_half_open.bed regional_sample_methylation.tsv patient_effects.tsv regional_mean_effects.tsv gene_overlaps.tsv nearest_tss.tsv inference_status.txt; do
    [[ ! -e "$run/comparisons/$key/$result" ]] || { echo "Existing result requires review: $result" >&2; exit 1; }
done
runtime="$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["runtime"])' "$run/config.json")"
[[ -f "$runtime/dmrcate.sif" && -f "$runtime/STAGED" ]]
sha256sum -c "$run/code.sha256"
python - "$workflow" "$validation" "$run" <<'PY'
import hashlib,json,sys
from pathlib import Path
w,v,r=map(Path,sys.argv[1:])
receipt=json.loads(v.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
if receipt['recovery_sha256']!=sha(w/'recover_checkpoint.R') or receipt['analysis_sha256']!=sha(r/'code/analysis.R'):
    raise SystemExit('Recovery or frozen analysis differs from tested source')
if receipt['sif_record']!=(Path(json.loads((r/'config.json').read_text())['runtime'])/'sif.sha256').read_text():
    raise SystemExit('Recovery runtime differs from tested runtime')
PY
python "$contract" check "$run"
if [[ "${DMRCATE_RECOVERY_WORKER:-}" != 1 ]]; then
    attempt_base="$(date -u +%Y-%m-%d_%H-%M)-checkpoint-$key"
    attempt_name="${DMRCATE_RECOVERY_NAME:-$attempt_base-v1}"
    require_date_first_result_name "$attempt_name"
    attempt="$run/recovery/$attempt_name"
    mkdir -p "$run/recovery"
    python "$contract" reserve "$attempt" --workflow "$workflow" --input "$validation" --input "$run/config.json" --work-root "$run"
    mkdir "$attempt/code"
    cp "$workflow/recover_checkpoint.R" "$workflow/recover_checkpoint.sh" "$attempt/code/"
    cp "$validation" "$attempt/validation.json"
    sha256sum "$attempt"/code/* > "$attempt/code.sha256"
    git -C "$project" rev-parse HEAD > "$attempt/source_revision.txt"
    date -u +submitted_utc=%FT%TZ > "$attempt/execution.tsv"
    job="$(sbatch --parsable --job-name=dmrcate-recovery --partition=longq --cpus-per-task=8 \
      --mem=294G --time=96:00:00 --export=ALL,DMRCATE_RECOVERY_WORKER=1 \
      --output="$attempt/slurm-%j.out" "$attempt/code/recover_checkpoint.sh" \
      "$run" "$key" "$attempt/validation.json" "$attempt")"
    printf '%s\n' "$job" > "$attempt/job_id"
    printf '%s\t%s\tcheckpoint-recovery\t%s\n' "$(date -u +%FT%TZ)" "$job" "$attempt" >> "$run/submission.tsv"
    printf 'Submitted recovery job %s\n%s\n' "$job" "$attempt"
    exit
fi
[[ -n "${SLURM_JOB_ID:-}" ]]
attempt="${4:?worker attempt required}"
exec 9>"$run/.lock"
flock -n 9 || { echo 'Run already active' >&2; exit 1; }
python "$contract" check "$run"
python "$contract" check "$attempt"
sha256sum -c "$attempt/code.sha256"
started="$(date -u +%FT%TZ)"
printf 'started_utc=%s\n' "$started" >> "$attempt/execution.tsv"
printf '%s\t%s\trecovery-started\n' "$SLURM_JOB_ID" "$started" >> "$run/attempts.tsv"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then python "$contract" fail "$attempt" || true; fi; printf "%s\t%s\trecovery-exit=%s\n" "$SLURM_JOB_ID" "$(date -u +%FT%TZ)" "$rc" >> "$run/attempts.tsv"; printf "exit_utc=%s\nexit_code=%s\n" "$(date -u +%FT%TZ)" "$rc" >> "$attempt/execution.tsv"' EXIT
module load containers/apptainer/1.3.6
[[ "$(apptainer --version)" == 'apptainer version 1.3.6' ]]
sha256sum -c "$runtime/sif.sha256"
# Retain the frozen input checks, including runtime/config/code/report hashes.
/usr/bin/time -v -o "$attempt/preparation-resources.txt" python "$run/code/prepare.py" \
  --config "$run/config.json" --run "$run" --resume
mkdir -p "$attempt/tmp/cache"
binds=(--bind "$project:$project:ro"
  --bind /path/to/hpc-home/.conda/envs/nextflow/share/wgbs_tools/references:/path/to/hpc-home/.conda/envs/nextflow/share/wgbs_tools/references:ro
  --bind /path/to/hpc-scratch/nextflow/nf-core-methylseq/results:/path/to/hpc-scratch/nextflow/nf-core-methylseq/results:ro
  --bind /path/to/hpc-scratch/nextflow/nf-core-methylseq/references:/path/to/hpc-scratch/nextflow/nf-core-methylseq/references:ro
  --bind "$runtime:$runtime:ro" --bind "$run:$run:rw")
# Each stage gets its own R process; smoothing never loads the count/model data.
for stage in smooth finish; do
  /usr/bin/time -v -o "$attempt/$stage-resources.txt" \
    apptainer exec --cleanenv --no-home --containall "${binds[@]}" \
    --env "TMPDIR=$attempt/tmp,TZ=UTC,OMP_NUM_THREADS=8,OPENBLAS_NUM_THREADS=1,XDG_CACHE_HOME=$attempt/tmp/cache" \
    "$runtime/dmrcate.sif" Rscript "$attempt/code/recover_checkpoint.R" "$stage" "$run" "$key" "$attempt"
done
date -u +completed_utc=%FT%TZ >> "$attempt/execution.tsv"

python "$contract" complete "$attempt"
