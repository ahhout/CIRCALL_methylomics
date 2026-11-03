# CIRCALL public source snapshot

The public repository is https://github.com/ahhout/CIRCALL_methylomics.
This snapshot publishes supported project launchers, source, tests, environment
definitions, reference definitions, generic configuration and operational
documentation. Pipeline implementation and its upstream/fork history remain in
https://github.com/ahhout/methylseq_AH.

This is a source snapshot. It does not import private local Git history, archived
analyses, per-run patient tables, derived patient results, logs, runtime/session
databases, installed dependencies, bulk sequencing/reference payloads or caches.
Reference symlinks record their HPC destinations; their targets are not uploaded.
Runtime paths and external software still require the recorded HPC environment.

## Public path templates

Automatic approval review rejected the initial publication commit because the
snapshot still contained coded sample identifiers and HPC filesystem paths.
The safer export replaces coded identifiers in documentation examples with
`SAMPLE_01`-style placeholders and converts site-specific home/storage prefixes
to `/path/to/...` templates. Local execution source, scientific settings and
completed output records are unchanged. Configure the path templates for the
intended environment before execution. Exported-file hashes are recorded in
the public manifest; original-source hashes are retained separately there.

The file and mode manifest is [publication_manifest.json](publication_manifest.json).

## Cohort tables and local state

The approved `metadata/metadata.tsv` and consolidated 12-sample
`samplesheets/samplesheet.csv` were already published on 2 November. This source
upload preserves their existing GitHub blobs without changing or expanding table
publication. `metadata/comparisons.tsv`, the retired separate eight-sample sheet,
patient-bearing downstream sample/configuration tables and historical copies
remain outside this upload. No additional table approval is inferred here.

The active HPC checkout still has its four-row production input sheet. The
12-sample consolidation is committed in the separate maintenance worktree and
has not been applied as part of the requested ignore/documentation housekeeping.
Existing frozen submissions keep their original source identity and inputs.

## Results and historical documentation

Result navigation and archived provenance are primarily local HPC records.
Some documentation links refer to retained local results, omitted patient
configuration or archived records and therefore do not resolve in a standalone
public checkout. Supply the original reviewed downstream sample tables from the
HPC checkout before invoking cohort-specific workflows; publication does not
change their scientific settings or authorize execution.

The [successful smoke validation](methylseq_smoke_validation.md) records bounded
execution with 10,000 read pairs. Full-cohort production, historical resume and
cross-layout cache reuse remain unverified.
