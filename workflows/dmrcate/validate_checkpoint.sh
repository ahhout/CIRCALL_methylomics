#!/usr/bin/env bash
#SBATCH --job-name=dmrcate-checkpoint-test
#SBATCH --partition=shortq
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
set -euo pipefail
workflow="${1:?workflow required}"
runtime="${2:?runtime required}"
dest="${3:?validation destination required}"
source /path/to/hpc-home/analysis/CIRCALL_methylomics/scripts/lib/result_names.sh
require_date_first_result_name "$(basename "$dest")"
[[ -n "${SLURM_JOB_ID:-}" && ! -e "$dest" ]]
contract="$workflow/../output_contract.py"
python "$contract" reserve "$dest" --workflow "$workflow" --input "$runtime/sif.sha256"
date -u +started_utc=%FT%TZ > "$dest/execution.tsv"
mkdir -p "$dest/code/config" "$dest/tmp"
cp "$workflow/analysis.R" "$workflow/recover_checkpoint.R" "$workflow/model_worker.R" "$workflow/test_checkpoint.R" "$workflow/prepare.py" "$workflow/test_prepare.py" "$workflow/run_dmrcate.sh" "$dest/code/"
cp "$workflow/config/cohort.json" "$workflow/config/samples.tsv" "$dest/code/config/"
sha256sum "$dest"/code/*.R > "$dest/code.sha256"
module load containers/apptainer/1.3.6
[[ "$(apptainer --version)" == 'apptainer version 1.3.6' ]]
sha256sum -c "$runtime/sif.sha256"
python "$dest/code/test_prepare.py" > "$dest/python-fixtures.txt" 2>&1
apptainer exec --cleanenv --no-home --containall --bind "$runtime:$runtime:ro" \
  --bind "$dest:$dest:rw" --env "TMPDIR=$dest/tmp,TZ=UTC,OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1" \
  "$runtime/dmrcate.sif" Rscript "$dest/code/test_checkpoint.R" "$dest/code" "$runtime" "$dest/fixture"
python - "$dest" "$runtime" <<'PY'
import hashlib,json,sys
from pathlib import Path
d,r=map(Path,sys.argv[1:])
v=json.loads((d/'fixture/checks.json').read_text())
for name,key in [('recover_checkpoint.R','recovery_sha256'),('analysis.R','analysis_sha256')]:
    v[key]=hashlib.sha256((d/'code'/name).read_bytes()).hexdigest()
v['sif_record']=(r/'sif.sha256').read_text()
(d/'VALIDATED.json').write_text(json.dumps(v,indent=2)+'\n')
(d/'r').mkdir()
names=['analysis.R','prepare.py','model_worker.R','recover_checkpoint.R','run_dmrcate.sh']
(d/'r/validated_code.tsv').write_text('file\tmd5\n'+''.join(n+'\t'+hashlib.md5((d/'code'/n).read_bytes()).hexdigest()+'\n' for n in names))
(d/'sif.sha256').write_text((r/'sif.sha256').read_text())
PY
date -u +completed_utc=%FT%TZ >> "$dest/execution.tsv"
touch "$dest/VALIDATED"

python "$contract" complete "$dest"
