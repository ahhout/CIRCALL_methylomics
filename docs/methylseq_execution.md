# Methylseq execution after reorganization

The CIRCALL project owns execution settings; `preprocessing/methylseq` owns
pipeline code. The project pins the exact pipeline commit in
`config/methylseq_pipeline.json`. Machine paths are declared once in
`config/methylseq_paths.json`; existing scratch work, container cache,
reference payload, result collection, and maintenance lock locations are
retained. Scientific settings remain in `config/params.emseq.integrated.yml`.

Commit reviewed changes in both repositories before running preflight. Keep
both at the same revisions throughout execution and any resume. The launcher
checks untracked source as well as tracked changes.

```bash
PROJECT=/path/to/hpc-home/analysis/CIRCALL_methylomics
RUN_NAME="$(date -u +%F_%H-%M)-emseq-integrated-circall-v1"
"$PROJECT/run_methylseq_integrated.sh" --preflight --run-name "$RUN_NAME"
# Save this exact name. Submission below is a separate operational action.
"$PROJECT/run_methylseq_integrated.sh" --submit --run-name "$RUN_NAME"
```

`--preflight` is read-only: it checks paths, input presence, committed revisions,
configuration, and the naming contract. It does not run Nextflow, reserve an
output, acquire the work-cache lock, or submit a job. The submission wrapper
passes both source commits to Slurm and assigns distinct submission/job logs
under `runs/methylseq/scheduler`. Direct `sbatch` without these recorded commit
identities is deliberately refused by the job entrypoint.

At execution the shell activates Conda `nextflow` and then loads module
`containers/apptainer/1.3.6`. The helper checks Nextflow 25.10.2 and Apptainer
1.3.6. It acquires the original shared work-cache lock before reserving output,
and keeps that lock through pipeline execution and control QC.

New outputs must not already exist, even as empty directories or dangling
symlinks. Atomic reservation establishes ownership before snapshots or payloads
are written. The original collection root is retained:

```text
/path/to/hpc-scratch/nextflow/nf-core-methylseq/results/bulk/<RUN_NAME>
```

Each new output starts a new Nextflow session without `-resume`, using the
existing shared work root and the new project runtime `runs/methylseq`.
Historical UUIDs and databases are not transplanted. Separate new sessions may
recompute tasks even though the work root is shared. Existing results are never
implicitly adopted by the new launcher.

To resume an incomplete new-layout run, retain its exact name and both commits:

```bash
"$PROJECT/run_methylseq_integrated.sh" --submit --run-name "$RUN_NAME" --resume-existing
```

Resume requires an owned incomplete output, intact original snapshots, identical
source/input identity, the recorded UUID, and matching history/cache in the new
runtime. A run that failed before enough session evidence existed cannot be
automatically resumed. Choose a new run name after investigating; no directory
is deleted automatically. Completed outputs are immutable through this launcher.

The closed command interface accepts no arbitrary Nextflow/config overrides.
Explicit `-C` configuration contains the pinned pipeline config, submitted
project config, and a per-attempt report-path config. Unrelated user/runtime
configs cannot redirect publication. Params cannot override `outdir` or remote
configuration selection. Existing upstream process publication rules use the
reserved `params.outdir`. The container cache is supplied from the path file.
The command disables the unused remote institutional-profile configuration
with `--custom_config_base ''`; the local Apptainer profile and committed
project Slurm settings remain selected. This prevents remote configuration
changes from redirecting publication.

`pipeline_info/submitted` contains immutable initial snapshots, source identities,
checksums, and run metadata. Each attempt has a unique job/attempt directory with
its exact argv, launcher snapshot, start/end records, Nextflow log, trace,
timeline, report, and graph. Run and attempt timestamps record actual execution
separately from naming/submission time. Control QC writes only under the owned
bundle. The project `outputs/<RUN_NAME>` navigation link is created exclusively;
a conflicting link is never repointed.

The control-QC script also enforces ownership when invoked directly. It requires
`--owned-output-root` and the launcher's `--attempt-id`, checks both clean source
commits and the submitted samplesheet checksum, and refuses historical or
completed bundles and symlink destinations. It writes its own execution
timestamps in that attempt's `control-qc.json`. Numerical summaries and the
lambda/pUC19 reporting thresholds are unchanged. Each attempt also records the
selected Nextflow/Apptainer executable paths, versions, module and environment.

The legacy `scripts/methylseq_cache_session.py` and source-owned old session JSON
remain historical recovery evidence. The new launcher does not import either.
Use unchanged original checkouts for historical resumes. Do not invoke legacy
cleanup/cache utilities against the new runtime; their assumptions concern the
original execution identity.

Initial migration validation used temporary fixtures, mocked execution, syntax
checks, and read-only preflight. Subsequent bounded execution is documented in
[the successful smoke record](methylseq_smoke_validation.md). Full-cohort
production, historical resume and cross-layout cache reuse remain unverified. No data payload was
copied, rehashed, removed, or regenerated by this reorganization.

## Bounded CIRCALL smoke execution

The subsequent authorized smoke test uses the first 10,000 paired FASTQ records
from the first current production samplesheet entry. It keeps that entry's
sample identifier, reference, EM-seq mode, clipping, extraction and QC settings.
It is an execution test, not a representative biological subsample. Genome-wide
cytosine reporting still scans the full production reference even with few reads.

Commit source first, then prepare the immutable input bundle and inspect it:

```bash
PROJECT=/path/to/hpc-home/analysis/CIRCALL_methylomics
RUN_NAME="$(date -u +%F_%H-%M)-emseq-integrated-layout-smoke-v1"
MANIFEST=$(python3 "$PROJECT/scripts/methylseq_smoke.py" --run-name "$RUN_NAME")
"$PROJECT/run_methylseq_integrated.sh" --preflight --run-name "$RUN_NAME" --smoke-inputs "$MANIFEST"
"$PROJECT/run_methylseq_integrated.sh" --submit --run-name "$RUN_NAME" --smoke-inputs "$MANIFEST"
```

The manifest, derived FASTQs and samplesheet stay ignored under
`runs/methylseq/smoke_inputs/<RUN_NAME>/`. Their hashes, original input identities
and production samplesheet hash are checked before execution. Originals are not
modified. The exact manifest, subset sheet, original sheet and resource config
are included in submitted snapshots. The main job uses shortq with a six-hour
limit; task resources are bounded by `config/methylseq_smoke_resources.config`.
The shared work root, maintenance lock and result parent are unchanged. The new
output starts without `-resume`; any later incomplete continuation must supply
the same manifest and explicit `--resume-existing`.

## Recorded smoke outcome

Job `6251849` completed successfully on 2 November 2026; all 12 tasks exited 0.
See [the smoke validation record](methylseq_smoke_validation.md) for frozen inputs,
source revisions, output paths, evidence hashes and validation limits.
