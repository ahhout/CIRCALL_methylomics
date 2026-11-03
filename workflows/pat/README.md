# BAM → PAT workflow

## Workspace transition and output ownership

This maintained workflow now belongs to `analysis/CIRCALL_methylomics`.
Use the original `/path/to/hpc-home/em_seq` entrypoints for historical resumes.
New runs reserve a previously nonexistent UTC-dated result directory; even an
empty directory or dangling symlink is a collision. Completed results are not
reused. Failed reservations and staging directories are retained for diagnosis.
The collection roots, references and scientific defaults have not changed.

See [the downstream migration audit](../../docs/reorganization/downstream_entrypoints.md)
for supported entrypoints, resume boundaries and fixture-only validation.

**Results:** [eight-sample baseline PAT bundle](../uxm/results/2026-09-09_21-50-emseq-integrated-8samples-trim-8-3-3-10-v1-pat-v1/)
· [newest cohort: four new HEH samples](../uxm/results/2026-09-24_20-56-emseq-integrated-heh-4new-trim-8-3-3-10-v1-pat-v1/).
These are separate cohorts. [Matching upstream and downstream results](../../RESULTS.md).

Supported for the current Flamingo setup, exact wgbstools hg38 dictionary and
Bismark methylseq inputs. Reuse and conversion preserve wgbstools 0.3.0,
samtools 1.22.1, hg38, LBETA and the existing bam2pat defaults: MAPQ 10,
excluded flags 1796, paired-end included flags 3, minimum one CpG and mate merge.
Use the existing `nextflow` Conda environment defined by the repository
[environment](../../environment.yml). No new environment is required.

## Prepare or reuse

From the repository root, choose a new run base name and explicit mode. This
example uses the frozen four-sample HEH metadata with the current production
sheet. The historical eight-sample cohort requires its own matched snapshots:

```bash
sbatch /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/pat/run_prepare_pat.sh \
  --results-dir /absolute/path/to/methylseq-result \
  --metadata /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/pat/runs/2026-09-16-heh-4new-inputs/metadata.tsv \
  --samplesheet /path/to/hpc-home/analysis/CIRCALL_methylomics/samplesheets/samplesheet.csv \
  --run-name 2026-09-29_14-35-cohort-analysis-v1 \
  --conversion-mode complete
```

`complete` generates all selected samples when no reuse bundle is supplied.
To generate only missing samples, supply repeated `--reuse-dir BUNDLE` and
`--reuse-sample SAMPLE` arguments naming exactly which samples to reuse.
`reuse` requires one or more `--reuse-dir BUNDLE` arguments and reuses the
entire selected cohort; omit `--reuse-sample` in that mode. Each reused sample
must have exactly one matching source. The supported launcher does not accept
legacy bundles or a forced-rebuild mode. Existing helper CLI compatibility is
retained for explicitly scoped historical consumers.

Metadata and samplesheet must contain identical unique sample IDs in the same
order and match the source BAM cohort. Reference, source identity, conversion
settings and version checks are enforced by the existing helpers. Preparation
runs with eight CPUs, 32 GB and 24 hours on `mediumq`.

Outputs remain under:

```text
/path/to/hpc-scratch/nextflow/nf-core-methylseq/results/uxm_deconvolution/2026-09-29_14-35-cohort-analysis-v1-pat-v1/
```

This historical physical root is shared with UXM; PAT source ownership is here.
The bundle contains `pat/` (PAT/CSI/LBETA), ordered `pat_manifest.tsv`,
`pat_provenance.json`, `PAT_COMPLETED`, checksummed submission snapshots and
validation logs. Existing sources are symlinked, never duplicated or overwritten.
Publication occurs on the same filesystem under an exclusive lock, after
validation. Success requires a successful Slurm exit and
`validation_after_publish.log`, not just the helper's `PAT_COMPLETED` marker.
`pat_provenance.json` and `pat_execution.json` record preparation start and
successful completion. `publication_execution.json` separately records the
launcher start, publication status, and end time; failures have `ended_utc`
without `completed_utc`. These are execution times, independent of folder naming.
A late publication collision retains work and leaves the destination unchanged.

Failures retain the hidden `.2026-09-29_14-35-cohort-analysis-v1-pat-v1.running-JOB` workspace. Inspect it
before taking further action; there is no automatic partial-run recovery or
Nextflow `-resume`. Existing final/work destinations are refused.

## Validate

Within a compute allocation:

```bash
python3 /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/pat/validate_pat_manifest.py \
  --results-dir /absolute/path/to/methylseq-result \
  --pat-dir /absolute/path/to/2026-09-29_14-35-cohort-analysis-v1-pat-v1 \
  --samplesheet /path/to/hpc-home/analysis/CIRCALL_methylomics/samplesheets/samplesheet.csv
```

Success prints `"status": "PASS"` and sample order. Checks cover source
identities, reference dictionary, conversion settings, indexed PAT integrity
and artifact checksums. Node-local device-ID changes trigger content hashing;
use `--hash-bams` through Slurm only when a fresh whole-BAM check is required.
Routine identity checks cannot detect content changes preserving all checked
metadata fields.

## Evidence

The [E31 bundle validation](runs/e31_pat_validation_20260910.md) records the
accepted eight-sample bundle. Its existing name must not be resubmitted.
The [original E31 migration recipe](analyses/e31_migration/README.md) preserves
its special three-sample legacy reconciliation unchanged. Source physical paths
inside frozen records remain historical; see the downstream migration manifest.
[Release evidence](../../docs/audits/downstream_structure_20260914/releases.md)
distinguishes existing scientific validation from new launcher-interface checks.
