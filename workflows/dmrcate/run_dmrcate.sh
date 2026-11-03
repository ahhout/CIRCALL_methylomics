#!/usr/bin/env bash
set -euo pipefail
workflow="${DMRCATE_WORKFLOW:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
project="$(cd "$workflow/../.." && pwd)"
contract="$project/workflows/output_contract.py"
source "$project/scripts/lib/result_names.sh"
config="${1:?usage: run_dmrcate.sh CONFIG RUN_NAME [--resume]}"
name="${2:?date-first UTC run name required}"
mode="${3:-}"
[[ -z "$mode" || "$mode" == --resume ]] || { echo 'Unsupported mode' >&2; exit 2; }
require_date_first_result_name "$name"
[[ -f "$config" && -f "$workflow/analysis.R" ]] || { echo 'Missing configuration or analysis script' >&2; exit 1; }
config="$(realpath "$config")"
python - "$config" "$workflow" <<'PY'
import csv,hashlib,json,sys
from pathlib import Path
c=json.load(open(sys.argv[1]))
validation=Path(c['runtime_validation'])
if not (validation/'VALIDATED').is_file():
    raise SystemExit('The frozen runtime has not passed its fixture validation')
if (validation/'sif.sha256').read_text().split()[0] != (Path(c['runtime'])/'sif.sha256').read_text().split()[0]:
    raise SystemExit('Validated runtime checksum does not match this configuration')
with (validation/'r/validated_code.tsv').open() as f:
    tested={r['file']:r['md5'] for r in csv.DictReader(f,delimiter='\t')}
for name in ('analysis.R','prepare.py','model_worker.R','recover_checkpoint.R','run_dmrcate.sh'):
    if tested.get(name) != hashlib.md5((Path(sys.argv[2])/name).read_bytes()).hexdigest():
        raise SystemExit('Current code differs from validated fixtures; validate it before submission: '+name)
PY
root="$(python -c 'import json,sys;print(json.load(open(sys.argv[1]))["scratch_root"])' "$config")"
run="$root/$name"
[[ "$root" == /path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate ]]
if [[ "${DMRCATE_WORKER:-}" != 1 ]]; then
    command -v sbatch >/dev/null
    # The bounded smoke analysis still scans full source reports and loads
    # annotations. Its allocation is separate from the full-cohort request.
    profile="$(python -c 'import json,sys; print("smoke" if json.load(open(sys.argv[1]))["smoke_interval"] is not None else "cohort")' "$config")"
    if [[ "$profile" == smoke ]]; then
        partition=mediumq; memory=32G; walltime=24:00:00
    else
        partition=mediumq; memory="${DMRCATE_MEMORY:-64G}"; walltime=24:00:00
    fi
    [[ "$memory" =~ ^[1-9][0-9]*G$ ]] || { echo 'DMRCATE_MEMORY must be a positive integer followed by G' >&2; exit 2; }
    if [[ "$mode" == --resume ]]; then
        [[ ! -e "$run/COMPLETE" ]] || { echo 'Completed run is immutable' >&2; exit 1; }
        python "$contract" resume "$run" --workflow "$workflow" --input "$config" --work-root "$run"
    else
        python "$contract" reserve "$run" --workflow "$workflow" --input "$config" --work-root "$run"
        cp "$config" "$run/config.json"
        mkdir -p "$run/code"
        cp "$workflow/prepare.py" "$workflow/analysis.R" "$workflow/model_worker.R" "$workflow/recover_checkpoint.R" "$workflow/environment.yml" "$workflow/container.def" "$workflow/run_dmrcate.sh" "$run/code/"
        git -C "$workflow" diff HEAD > "$run/code/tracked_changes.patch"
        git -C "$workflow" status --short > "$run/code/git_status.txt"
        sha256sum "$run"/code/* > "$run/code.sha256"
        date -u +naming_utc=%FT%TZ > "$run/submission.tsv"
    fi
    cmp "$config" "$run/config.json"
    cmp "$workflow/prepare.py" "$run/code/prepare.py"
    cmp "$workflow/analysis.R" "$run/code/analysis.R"
    cmp "$workflow/model_worker.R" "$run/code/model_worker.R"
    cmp "$workflow/recover_checkpoint.R" "$run/code/recover_checkpoint.R"
    date -u +submitted_utc=%FT%TZ >> "$run/submission.tsv"
    git -C "$workflow" rev-parse HEAD >> "$run/submission.tsv"
    printf '%s\tprofile=%s\tcpus=8\tmem=%s\ttime=%s\tpartition=%s\n' \
        "$(date -u +%FT%TZ)" "$profile" "$memory" "$walltime" "$partition" >> "$run/resource_requests.tsv"
    mkdir -p "$project/runs/dmrcate"
    sbatch --job-name=dmrcate --partition="$partition" --cpus-per-task=8 --mem="$memory" --time="$walltime" \
        --export="ALL,DMRCATE_WORKER=1,DMRCATE_WORKFLOW=$workflow" --output="$project/runs/dmrcate/${name}-%j.out" \
        "$workflow/run_dmrcate.sh" "$config" "$name" "$mode" | tee -a "$run/submission.tsv"
    exit
fi
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'Analysis must run through sbatch' >&2; exit 1; }
cmp "$config" "$run/config.json"
config="$run/config.json"
exec 9>"$run/.lock"
flock -n 9 || { echo 'Run already active' >&2; exit 1; }
[[ ! -e "$run/COMPLETE" ]]
python "$contract" check "$run"
sha256sum -c "$run/code.sha256"
started="$(date -u +%FT%TZ)"
[[ -f "$run/started_utc" ]] || printf '%s\n' "$started" > "$run/started_utc"
printf '%s\t%s\tstarted\n' "$SLURM_JOB_ID" "$started" >> "$run/attempts.tsv"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then python "$contract" fail "$run" || true; fi; printf "%s\t%s\texit=%s\n" "$SLURM_JOB_ID" "$(date -u +%FT%TZ)" "$rc" >> "$run/attempts.tsv"' EXIT
module load containers/apptainer/1.3.6
[[ "$(apptainer --version)" == 'apptainer version 1.3.6' ]]
runtime="$(python -c 'import json,sys;print(json.load(open(sys.argv[1]))["runtime"])' "$config")"
mkdir -p "$run/tmp/cache"
export TMPDIR="$run/tmp"
args=(--config "$config" --run "$run")
[[ "$mode" != --resume ]] || args+=(--resume)
# Bind the full project/reference/source trees read-only, with only this run writable.
# No mutable R library, user home, installation or network annotation helper is used.
binds=(--bind "$project:$project:ro"
 --bind /path/to/hpc-home/.conda/envs/nextflow/share/wgbs_tools/references:/path/to/hpc-home/.conda/envs/nextflow/share/wgbs_tools/references:ro
 --bind /path/to/hpc-scratch/nextflow/nf-core-methylseq/results:/path/to/hpc-scratch/nextflow/nf-core-methylseq/results:ro
 --bind /path/to/hpc-scratch/nextflow/nf-core-methylseq/references:/path/to/hpc-scratch/nextflow/nf-core-methylseq/references:ro
 --bind "$runtime:$runtime:ro" --bind "$run:$run:rw")
# Run this before scanning twelve whole-genome reports. --cleanenv requires
# explicit propagation of TMPDIR; host export alone does not reach R.
apptainer exec --cleanenv --no-home --containall "${binds[@]}" \
 --env "TMPDIR=$run/tmp,TZ=UTC,DMRCATE_TEST_MODE=1" \
 "$runtime/dmrcate.sif" Rscript -e \
 'a<-commandArgs(TRUE); source(a[1]); runtime_preflight(a[2],a[3])' \
 "$run/code/analysis.R" "$config" "$run"
/usr/bin/time -v -o "$run/preparation-resources-$SLURM_JOB_ID.txt" python "$run/code/prepare.py" "${args[@]}"
run_r() {
    local resources="$1"
    shift
    /usr/bin/time -v -o "$resources" \
      apptainer exec --cleanenv --no-home --containall "${binds[@]}" \
      --env "TMPDIR=$run/tmp,OMP_NUM_THREADS=8,OPENBLAS_NUM_THREADS=1,TZ=UTC,http_proxy=http://127.0.0.1:9,https_proxy=http://127.0.0.1:9,XDG_CACHE_HOME=$run/tmp/cache" \
      "$runtime/dmrcate.sif" Rscript "$@"
}
while IFS= read -r key; do
    comparison="$run/comparisons/$key"
    [[ ! -f "$comparison/INFERENCE_READY.rds" ]] || continue
    phase="$run/phases/$key/attempt-$SLURM_JOB_ID"
    mkdir -p "$phase"
    if [[ ! -f "$comparison/cpgannotated.rds" ]]; then
        run_r "$phase/model-resources.txt" "$run/code/model_worker.R" "$config" "$run" "$key" "$phase"
    fi
    # A contrast with fewer than two eligible CpGs is finalized by the model worker.
    if [[ ! -f "$comparison/INFERENCE_READY.rds" ]]; then
        run_r "$phase/smooth-resources.txt" "$run/code/recover_checkpoint.R" smooth "$run" "$key" "$phase"
        run_r "$phase/finish-resources.txt" "$run/code/recover_checkpoint.R" finish "$run" "$key" "$phase"
    fi
done < <(python -c 'import json,sys; c=json.load(open(sys.argv[1])); print("\n".join("_minus_".join(x) for x in c["contrasts"]))' "$config")
run_r "$run/validation-resources-$SLURM_JOB_ID.txt" -e \
 'Sys.setenv(DMRCATE_TEST_MODE="1"); a<-commandArgs(TRUE); source(a[1]); validate_outputs(a[2],a[3])' \
 "$run/code/analysis.R" "$config" "$run"
date -u +completed_utc=%FT%TZ > "$run/COMPLETE"
python "$contract" complete "$run"
