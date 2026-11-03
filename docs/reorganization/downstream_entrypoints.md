# Downstream entrypoints after the workspace transition

Only active source was adapted. Frozen source under component `runs/`, `analyses/`
and the historical UXM target export stays unchanged. Full workflow executions
still require Slurm. This migration submitted no jobs and did not touch historical
runtime state, reference payloads or cache contents.

| Component / public entrypoint | Output and recovery contract |
|---|---|
| PAT `run_prepare_pat.sh` / `launch_pat.py` | Reserve final leaf before stage writes, hold persistent destination lock, publish each child without replacement, mark success after postpublication validation. No resume/adoption of failed or completed bundles. |
| PAT `prepare_pat.py` | Atomic absent leaf; standalone CLI requires a dated canonical destination. Internal prepared stages require `--owner-output` or `EXPERIMENT_OWNER_OUTPUT`/`CIRCALL_OUTPUT_OWNER`. Existing empty leaves are rejected. Submitted helper copies must match the clean committed project source; historical snapshots retain their original entrypoints. |
| UXM full/hematopoietic launchers | Reserve final leaf under destination lock, preserve staging on failure, publish without replacing existing children. No resume of completed or failed bundles. |
| UXM `export_u25_targets.py` | Dated exclusive owned reference-export bundle; source atlas and dictionary identities unchanged. |
| EpiDISH comparison launcher | Reserve final leaf, retain existing frozen input/receipt recovery contract, allow explicit incomplete recovery only with matching final owner and work contract. Completed output refused. |
| EpiDISH standalone plot / saved-estimate export | New owned dated bundles only. Standalone figures default to a sibling of the source result; they are not written into completed source bundles. |
| EpiDISH `refresh_composition_plots.R` | Requires `COMPLETED_RESULT_DIR NEW_OUTPUT_DIR`. Exports the same saved-estimate composition plots and new provenance separately; never rewrites source plots or manifests. |
| DecoNFlow analysis and CelFiE patient entrypoints | `common.new_run` reserves ownership; explicit incomplete resume verifies source/input/runtime fingerprints. Per-launch UUID receipt must match `.nextflow.log`, cache DB and history; resumes pass the exact UUID. Shared Nextflow work cache uses the existing maintenance lock. |
| DecoNFlow reference build, runtime/route fixtures and comparisons | New owned dated bundles through `common.new_run`; completion becomes immutable. Runtime staging is an explicit separate operation and was not run during migration. |
| DecoNFlow matched-reference/count preparation | New owned dated bundle; existing incomplete output requires explicit `--resume` plus unchanged settings/input identity. Completed or terminal blocked results cannot be reused. |
| DecoNFlow direct `pat_counts.py` CLI | New owned dated bundle only; implicit cached directory adoption is refused. Internal count-library reuse remains inside its parent launcher’s owned incomplete execution. |
| DecoNFlow matched-count validator / observed-selection exporter | New owned dated bundles only; output and source provenance retained. Diagnostic observers are internal task helpers and exclusively create diagnostic files. |
| DecoNFlow `validate_pat_route.py` | Requires `--output NEW_DATED_BUNDLE` for read-only validation of completed results. Cache replay against a completed bundle is refused; cache checks remain part of an incomplete execution. |
| DMRcate run launcher | Reserve new output atomically or verify explicitly owned incomplete resume; config/code snapshots and input fingerprints remain required. Project-local job logs are separate from output reservation. |
| DMRcate checkpoint recovery | Requires an owned incomplete parent; reserves a distinct recovery attempt. `DMRCATE_RECOVERY_NAME` can supply the exact attempt name; collisions fail, without automatic version increments. |
| DMRcate runtime/checkpoint validation and healthy-WBC comparison | Exclusive owned dated bundles. Existing runtime-validation code hashes intentionally invalidate when source changes; rerun approved fixtures before any later real execution. |

Internal R/Python stages operate within their launcher's incomplete bundle and do
not authorize independent writes into completed results. Library functions used
by synthetic tests do not submit scheduler jobs. Sourceable plotting functions
are called only inside an already-reserved export's internal sections.

The ownership adapter hashes only small source/configuration manifests, pruning
runtime/output/dependency directories before traversal. The existing pinned atlas,
UXM/DecoNFlow software, wgbstools dictionary, BCP-ALL alignment reference and scratch
result parents remain external dependencies. All changes are path/ownership and
provenance changes; this migration does not grant new scientific validation.

Publication uses exclusive no-replace operations. BeeGFS does not support Linux
`renameat2(RENAME_NOREPLACE)` on this mount; the shared guard falls back to exclusive
directory creation and same-filesystem hard links/symlinks, removing only source
links after transfer. It never copies large payloads across filesystems. Empty
staging directories are retained as recovery evidence.

Supported execution requires a clean committed project repository. Staged, tracked
and nonignored untracked changes fail before new reservation or continued execution;
ignored output payloads do not count as source edits. Resume requires the original
project revision as well as the existing source/input/session checks.
