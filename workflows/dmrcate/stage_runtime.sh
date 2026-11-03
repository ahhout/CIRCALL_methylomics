#!/usr/bin/env bash
# Submit through sbatch. Public downloads happen here, never in analysis.
#SBATCH --job-name=dmrcate-stage
#SBATCH --partition=shortq
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=06:00:00
set -euo pipefail
workflow="${1:?workflow directory required}"
root="${2:?date-first staging directory required}"
source "$workflow/../../scripts/lib/result_names.sh"
require_date_first_result_name "$(basename "$root")"
resume="${4:-}"
if [[ "$resume" == --resume-annotations ]]; then
    [[ -f "$root/sif.sha256" && ! -e "$root/STAGED" ]]
    sha256sum -c "$root/sif.sha256"
    date -u +annotation_resume_started_utc=%FT%TZ >> "$root/execution.tsv"
elif [[ -e "$root" ]]; then
    echo "Refusing existing stage without --resume-annotations: $root" >&2; exit 1
fi
mkdir -p "$root/sources" "$root/annotations" "$root/tmp" "$root/cache"
if [[ "$resume" != --resume-annotations ]]; then
date -u +started_utc=%Y-%m-%dT%H:%M:%SZ > "$root/execution.tsv"
module load containers/apptainer/1.3.6
[[ "$(apptainer --version)" == 'apptainer version 1.3.6' ]]
export APPTAINER_CACHEDIR="${3:-$root/cache}"
# Container extraction is metadata-heavy. Keep transient build files on the
# allocated node; the pulled immutable image and all provenance stay on scratch.
export APPTAINER_TMPDIR="${SLURM_TMPDIR:-/tmp}/dmrcate-build-${SLURM_JOB_ID:?}"
mkdir -p "$APPTAINER_TMPDIR"
printf 'build_temporary_directory\t%s\n' "$APPTAINER_TMPDIR" >> "$root/execution.tsv"
cp "$workflow/container.def" "$workflow/environment.yml" "$root/sources/"
cd "$root/sources"
image="$(sed -n 's/^From: //p' container.def)"
[[ "$image" == *@sha256:* ]] || { echo 'Container image must be digest-pinned' >&2; exit 1; }
apptainer pull "$root/dmrcate.sif" "docker://$image"
printf '%s\t%s\t%s\n' "$(date -u +%FT%TZ)" "$image" "$(sha256sum "$root/dmrcate.sif")" > downloads.tsv
apptainer exec --cleanenv --no-home "$root/dmrcate.sif" Rscript -e 'p<-installed.packages(); write.table(p[,c("Package","Version","Built")],stdout(),sep="\t",row.names=FALSE,quote=FALSE)' > "$root/packages.lock.tsv"
sha256sum "$root/dmrcate.sif" > "$root/sif.sha256"
apptainer exec --cleanenv --no-home "$root/dmrcate.sif" python -c \
 'import glob,json; print(json.dumps([json.load(open(p)) for p in sorted(glob.glob("/usr/local/conda-meta/*.json"))],indent=2))' > "$root/conda-packages.lock.json"
fi
cd "$root/annotations"
name=gencode.v50.basic.annotation.gtf.gz
url=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz
    [[ ! -e "$name" ]] || { echo "Refusing an unreviewed existing annotation: $name" >&2; exit 1; }
    # EBI's HTTPS endpoint can refuse connections from Flamingo. Its official
    # HTTP endpoint serves the same release-specific file; record the actual URL.
    if ! curl --fail --location --retry 2 --connect-timeout 15 "$url" -o "$name.part"; then
        [[ "$url" == https://ftp.ebi.ac.uk/* ]] || exit 1
        url="http://${url#https://}"
        curl --fail --location --retry 2 --connect-timeout 15 "$url" -o "$name.part"
    fi
    gzip -t "$name.part"
    mv "$name.part" "$name"
    printf '%s\t%s\t%s\n' "$(date -u +%FT%TZ)" "$url" "$(sha256sum "$name")" >> sources.tsv
sha256sum ./*.gz > SHA256SUMS
date -u +completed_utc=%Y-%m-%dT%H:%M:%SZ >> "$root/execution.tsv"
touch "$root/STAGED"
