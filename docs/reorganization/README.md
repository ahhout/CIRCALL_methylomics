# EM-seq staged reorganization — 2 November 2026

The project and pipeline are independent clones on `reorganize/emseq-layout`.
Four experiments are separate local snapshot repositories. Original checkouts,
refs, runtime databases, scratch payloads, references and environments remain in
place. No jobs, cleanup, publication, or biological changes were performed.

## Records

- `original_state.json`: original HEADs/refs, source hashes, and selected runtime
  metadata fingerprints.
- `original_preservation_check.json`: verification of unchanged originals, retained
  links, local dependency patches, and selected runtime metadata.
- `clone_verification.json`: preserved refs and independent Git object storage.
- `untracked_inventory.json`: reviewed untracked source/provenance and retained
  untracked artifacts; the path map records their final ownership.
- `path_map.json`: each original tracked/imported source location and its new or
  retained home.
- `external_dependencies.json`: retained reference/software/navigation targets.
- `DMRichR-local.patch`: local modifications in the ignored archived dependency.
- `methylseq_fork_archive.json`: exact archived CIRCALL integration files removed
  from the new pipeline checkout.
- `methylseq_settings_reconciliation.json`: comparison of old project and fork
  settings, including allowed path substitutions.
- `verification.json`: executed checks, results and limitations; detailed commands
  are in `methylseq_verification.json`, `downstream_verification.json`, and
  `experiment_verification.json`.
- `dmrcate_runtime_validation.json`: fresh source-matched synthetic R validation
  using the retained pinned runtime; the failed sandbox attempt is preserved.
- `downstream_entrypoints.md`: changed downstream invocation interfaces.
- Each experiment owns `MIGRATION.md`, `MIGRATION_SOURCE.json`, `ENTRYPOINTS.tsv`,
  `EXTERNAL_DEPENDENCIES.json`, and `PROJECT_DEPENDENCY.json`.

## Transition and rollback

Future source development uses the new layout. Historical resumes continue from
the unchanged old entrypoints. Never adopt their output directories or session
UUIDs through a new launcher. New outputs must be exclusively reserved and carry
ownership; completed outputs cannot be changed by supported launchers.

If new-layout validation fails, retain the new source for review and continue
with the old entrypoints. Original locations must not be retired until their
remaining reference/software dependencies and historical recovery needs have
been separately resolved. Output protections cover supported entrypoints, not
manual filesystem edits.

## Validation boundary

Only temporary fixtures and non-submitting checks are part of this migration.
Production execution, cross-layout cache reuse and historical resume from new
paths have not been validated. Sarek/Oncoanalyser, scratch storage relocation,
cache retirement, pipeline upgrades and BWA-MEM3 integration are separate work.

## Retained dependencies

The original checkouts remain required for historical sessions, ignored UXM and
DecoNFlow software, patched DMRichR/ggbiplot, public PAT reference payloads,
private copy-number metadata, and historical local EpiDISH exports. These are
recorded navigation/dependency links, not copied data. BCP-ALL reference payloads
remain at their recorded scratch paths. Existing Conda environments, containers,
module Apptainer, reference payloads and shared work/cache roots remain external.
Maintained workflows do not import implementation code from brainstorming.

The archived DMRichR patch, imported audit tables, and fork version receipt keep
their original whitespace and line endings. Active changed source passes the
whitespace check; preserved provenance is excluded from that formatting check.

## Subsequent execution validation

The original migration boundary above describes the staged reorganization.
The subsequently authorized 10,000-pair CIRCALL smoke run completed successfully
on 2 November 2026. See [the verified execution record](../methylseq_smoke_validation.md).
Full-cohort production and historical resume from the new layout remain unverified.
