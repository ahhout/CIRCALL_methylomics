# CIRCALL methylomics

This repository contains CIRCALL methylation study configuration, preprocessing
invocation, downstream workflows, and run provenance. The scientific questions,
cohort definitions, and interpretation limits remain in [PROJECT.md](PROJECT.md).
[RESULTS.md](RESULTS.md) indexes the retained results; no analysis was rerun for
this reorganization.

## Layout and supported entrypoints

- Pipeline implementation: `/path/to/hpc-home/preprocessing/methylseq`, preserving
  the existing nf-core/methylseq 4.2.0 baseline and fork history.
- This project owns `config/`, `samplesheets/`, `metadata/`, `environment.yml`,
  reference definitions, and `run_methylseq_integrated.sh`.
- Downstream implementations are in `workflows/pat`, `workflows/uxm`,
  `workflows/epidish`, `workflows/deconflow`, and `workflows/dmrcate`. Component
  READMEs retain scientific scope and validation limits. Their location does not
  establish biological validation or acceptance as the default method.
- BCP-ALL atlas, BWA-MEM3 pilot, PCDH EpiClone, and W footprint are separate
  exploratory repositories under `/path/to/hpc-home/brainstorming`.
- `runs/<workflow>/` owns new execution state; `outputs/` provides navigation.
  Bulk inputs, references, results and work caches remain at existing HPC paths.
- `archive/preprocessing/` and `archive/downstream/` retain earlier source and
  evidence. Historical records retain their original paths and timestamps.

## Methylseq execution

Read [the execution contract](docs/methylseq_execution.md). The path configuration
is `config/methylseq_paths.json`; `config/methylseq_pipeline.json` pins pipeline
identity. Production remains Bismark with `em_seq: true`, the exact wgbstools hg38
sequence plus lambda/pUC19, Nextflow 25.10.2, and module Apptainer 1.3.6.

From any working directory, inspect a proposed future run without submitting it:

```bash
/path/to/hpc-home/analysis/CIRCALL_methylomics/run_methylseq_integrated.sh \
  --preflight --run-name 2026-11-02_14-30-emseq-integrated-example-v1
```

To execute a future authorized run, use `--submit` instead of `--preflight`.
Do not reuse the example date/name; choose an explicit UTC run name for that
execution. Full pipelines execute only through Slurm. The initial reorganization used fixtures and non-submitting checks. The
subsequent [10,000-pair smoke execution](docs/methylseq_smoke_validation.md)
completed successfully; full-cohort production and historical resume remain
unverified.

## Output ownership and resume

Unchanged results parents do not mean reused run directories. Every new bundle
requires an atomic reservation of a nonexistent dated child. Empty directories,
files, aliases and dangling symlinks are collisions. There is no overwrite or
automatic renaming mode. Failure retains the reservation and recovery evidence.

Only an explicitly owned incomplete run can resume after checking submitted
inputs, source identity, work root and session. Completed bundles are immutable
through supported entrypoints; reruns and figure refreshes create new bundles.
Result links cannot be repointed to another bundle. See the component contracts
for staged downstream execution.

The first future methylseq run uses new project-owned launch state; it does not
adopt the old session database. Cache reuse after source relocation is unverified.
Historical resumes remain at `/path/to/hpc-home/em_seq` and
`/path/to/hpc-home/methylseq_AH`, using their original launchers and snapshots.
Those directories, existing scratch results/work, references, and installed
environments have not been moved. Do not edit frozen records or select runs by
mtime. New run names use `YYYY-MM-DD_HH-MM-description-vN` in UTC; retain original
timestamps on resume and append attempt records.

## Reorganization evidence

[Migration records](docs/reorganization/README.md) document checkpoints, file
ownership, ignored external dependencies, source preservation, and fixture
validation. Git commits protect source and appropriate small records, not
external inputs, installed tools, outputs or caches. The public source snapshot
and local-only exclusions are described in [publication scope](docs/publication.md).

## Completed smoke validation

Job `6251849` completed on 2 November 2026 with all 12 tasks successful.
See [the execution evidence and limits](docs/methylseq_smoke_validation.md).
