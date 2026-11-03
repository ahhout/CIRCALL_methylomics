#!/usr/bin/env bash
#SBATCH --job-name=dmrcate-fixtures
#SBATCH --partition=shortq
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
set -euo pipefail
workflow="${1:?workflow path required}"
runtime="${2:?staged runtime required}"
dest="${3:?date-first validation bundle required}"
source "$workflow/../../scripts/lib/result_names.sh"
require_date_first_result_name "$(basename "$dest")"
[[ -f "$runtime/STAGED" && ! -e "$dest" ]]
python "$workflow/../output_contract.py" reserve "$dest" --workflow "$workflow" --input "$runtime/sif.sha256"
mkdir "$dest/cache" "$dest/tmp"
date -u +started_utc=%FT%TZ > "$dest/execution.tsv"
# Validate a frozen source snapshot, including the selected output contract.
mkdir -p "$dest/code"
cp "$workflow/analysis.R" "$workflow/prepare.py" "$workflow/test_fixtures.R" "$workflow/test_prepare.py" "$workflow/model_worker.R" "$workflow/recover_checkpoint.R" "$workflow/test_checkpoint.R" "$workflow/run_dmrcate.sh" "$dest/code/"
cp -r "$workflow/config" "$dest/code/config"
contract="$workflow/../output_contract.py"
workflow="$dest/code"
module load containers/apptainer/1.3.6
[[ "$(apptainer --version)" == 'apptainer version 1.3.6' ]]
sha256sum -c "$runtime/sif.sha256"
# Preserve complete dependency artifact URLs and hashes from the digest-pinned base.
apptainer exec --cleanenv --no-home "$runtime/dmrcate.sif" python -c \
 'import glob,json; print(json.dumps([json.load(open(p)) for p in sorted(glob.glob("/usr/local/conda-meta/*.json"))],indent=2))' > "$dest/conda-packages.lock.json"
cp "$runtime/packages.lock.tsv" "$runtime/sif.sha256" "$dest/"
python "$workflow/test_prepare.py" > "$dest/python-fixtures.txt" 2>&1
apptainer exec --cleanenv --no-home --containall --bind "$workflow:$workflow:ro" \
  --bind "$runtime:$runtime:ro" --bind "$dest:$dest:rw" \
  --env "TMPDIR=$dest/tmp,http_proxy=http://127.0.0.1:9,https_proxy=http://127.0.0.1:9,TZ=UTC,XDG_CACHE_HOME=$dest/cache" \
  "$runtime/dmrcate.sif" Rscript "$workflow/test_fixtures.R" "$workflow" "$dest/r"
apptainer exec --cleanenv --no-home --containall --bind "$workflow:$workflow:ro" \
  --bind "$runtime:$runtime:ro" --bind "$dest:$dest:rw" \
  --env "TMPDIR=$dest/tmp,http_proxy=http://127.0.0.1:9,https_proxy=http://127.0.0.1:9,TZ=UTC,XDG_CACHE_HOME=$dest/cache" \
  "$runtime/dmrcate.sif" Rscript "$workflow/test_checkpoint.R" "$workflow" "$runtime" "$dest/fresh-process"
date -u +completed_utc=%FT%TZ >> "$dest/execution.tsv"
touch "$dest/VALIDATED"

python "$contract" complete "$dest"
