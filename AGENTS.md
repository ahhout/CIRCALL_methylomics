# Repository instructions

## Project context

- This repository is CIRCALL_methylomics: project configuration, preprocessing invocation, and downstream analysis in the Flamingo/HPC environment.
- The intended pipeline release is `nf-core/methylseq` `4.2.0`.
- The intended runtime is Slurm with Apptainer containers.
- The Flamingo SSH login is `ssh a_houtman@flamingo`.
- The Nextflow Conda environment is `nextflow`.
- The Apptainer module currently used by the launchers is `containers/apptainer/1.3.6`.
- Although `environment.yml` also records an Apptainer package, methylseq
  launchers must load and use the module-provided Apptainer 1.3.6 executable.
- Do not create or reference virtualenv/venv or `.venv`.
- Always consult `env.yml` / `environment.yml` if present.
- If adding dependencies, update the environment file and avoid unnecessary packages.

## Canonical entry points

- `run_methylseq_integrated.sh` is the sole supported production launcher.
- `config/params.emseq.integrated.yml` is the production methylseq parameter
  source.
- `samplesheets/samplesheet.csv` is the production input samplesheet.
- `config/nextflow.config` is the Slurm and Apptainer configuration source.
- `references/emseq_hg38_lambda_puc19` exposes the production reference: the
  exact wgbstools hg38 sequence plus lambda and pUC19.
- `README.md` owns operational usage, `PROJECT.md` owns scientific scope and
  interpretation, and `references/README.md` owns reference provenance.
- Files under `archive/` and `archive/downstream/` are provenance records, not
  supported entry points. Restoring the nf-core/iGenomes GRCh38-only route
  requires an explicit workflow-policy decision.

## Parameter and decision handling

- Treat methylation workflow settings and biological parameters as design choices.
- Do not silently choose or change `genome`, `fasta`, `bismark_index`, `aligner`, `em_seq`, `taps`, `pbat`, `rrbs`, sample sheet structure, trimming behavior, methylation extraction options, deduplication behavior, target-region settings, or sample merging.
- Confirm choices with the user before changing anything that affects interpretation or reproducibility.
- Explicit user authorization in the current conversation satisfies confirmation
  requirements for the specified action and scope; do not request the same
  approval again. Proceed with routine implementation and read-only investigation
  within the authorized scope. Ask before making an unapproved scientific choice
  or performing an action covered by the confirmation requirements below. While
  awaiting clarification, continue independent work.
- Keep production parameters in the canonical params file unless there is a
  clear reason to expose a setting in the launcher.
- Keep cluster/runtime settings in `config/nextflow.config`, not in the params file.
- For the built-in nf-core test profile, prefer letting `-profile test` provide its own test inputs and parameters.
- Clearly state assumptions when giving commands or making changes.

## Reproducibility

- Prefer deterministic behavior and explicit versions.
- Keep the methylseq release pinned to `4.2.0` unless the user asks to change it.
- Expose hidden defaults when they matter for reproducibility.
- Avoid overwriting or deleting raw data by default.
- Keep pipeline, config, params, sample sheet, and reference-data changes easy to audit.
- Production integrated EM-seq currently uses Bismark with `em_seq: true`.

## Launcher and script requirements

- Run full pipeline executions only via `sbatch`.
- Run small scripts, checks and quick debugging commands directly on the current
  node when they need at most 16 GB of memory and at most 8 CPUs; no `sbatch` or
  `sinteractive` request is needed for these tasks.
- For scripts and checks, request Slurm resources only when the task needs more
  than 16 GB of memory or more than 8 CPUs. Full pipeline executions still use
  `sbatch` as required above.
- Do not run full methylseq pipeline commands directly on a login node.
- Bash launchers and local helper scripts must use `set -euo pipefail`.
- Do not inject `set -o pipefail`, `set -euo pipefail`, or other shell prologues into nf-core/methylseq internals, generated `.command.sh` task scripts, or process definitions.
- Quote shell variables.
- Validate required files and commands early, including:
  - the script directory
  - work and output directories
  - `config/nextflow.config`
  - the canonical production params file
  - the canonical production samplesheet
  - `conda`
  - `nextflow`
  - `apptainer`
  - local reference paths before production runs if custom references are added
- Fail loudly with meaningful errors.
- Python helper scripts should provide a CLI interface with `argparse`; this does not apply to Slurm launcher scripts.

## Bioinformatics data conventions

- Assume datasets may be large, especially FASTQ, BAM, CRAM, bedGraph, methylation calls, and MultiQC outputs.
- Prefer streaming or indexed access when writing helper scripts.
- Preserve headers in BAM/CRAM-derived and tabular outputs where applicable.
- Explicitly specify genome build, such as GRCh38, for production runs.
- Be aware of contig naming differences, such as `chr1` versus `1`.
- Avoid coordinate conversions, target-region restrictions, sample merges, or methylation-call filtering without confirmation.

## Data and pipeline safety

- Assume clinical or sensitive data may exist.
- Scratch is finite. Keep one shared Nextflow work cache for resumability and
  do not create a separate work directory for each sequencing run.
- Do not duplicate large FASTQ, BAM, or PAT files when an authoritative path
  or symlink is sufficient.
- Do not delete individual files from a Nextflow task directory intended for
  `-resume`. No automatic cleanup is enabled. If cleanup is approved later,
  prioritize abandoned tasks, obsolete validated work caches, then unused
  container images.
- Require confirmation before:
  - data deletion or overwriting
  - changing genome build or reference resources
  - changing aligner or library mode
  - enabling GPU/Parabricks mode
  - enabling TAPS, PBAT, RRBS, or targeted sequencing modes
  - merging samples or changing sample identifiers
  - changing production output or work directories
  - changing methylation extraction, deduplication, or filtering behavior

## Documentation

- All new run-level output/result directories must start with the date, hour
  and minute in UTC: `YYYY-MM-DD_HH-MM-<analysis>-<cohort-or-purpose>-vN`.
  Validate both the calendar date and 24-hour time before creating outputs.
  Record the run-start timestamp and keep it unchanged on resume; never guess
  the latest result from mtime. Older date-only names may remain for existing
  results/recovery when no trustworthy recorded time is available. Do not
  invent `00-00` or use the rename time to fill missing historical times.
- Apply this rule to methylseq, every downstream workflow/project, and new
  figure/export bundles, including direct helper invocations. Collection roots
  (`outputs`, `results`, workflow names), internal sections (`summary`,
  `comparisons`, `prepared`), per-sample directories, references and recovery
  workspaces are not run names and do not need date prefixes.
- Every custom pipeline, experimental launcher and standalone analysis/plotting
  script must save an explicit UTC execution-start timestamp in its run
  manifest/provenance with at least hour-and-minute precision. Record completion
  time on success too; seconds may be retained. Use timezone-labelled fields
  such as `started_utc` and `completed_utc` with ISO 8601 values. Folder names,
  filesystem mtimes and scheduler history alone are not sufficient provenance.
  Distinguish submission/naming time from actual execution start. On resume,
  retain the original run timestamp and append timestamps for each attempt.
- Keep analysis IDs/sample IDs, scientific settings, original run manifests
  and historical commands unchanged when renaming results. Record an exact
  old/new path map and date source; update active navigation and retain
  compatibility links where old paths are embedded in provenance or consumers.
- Follow the naming contract and examples in `README.md` under
  "Output and result folder names" when adding or changing a launcher.

- Write project READMEs in plain language so someone returning a year later
  can quickly understand the scientific question, project context and biology
  without the conversation history. Define specialist terms on first use.
- Keep the main thresholds and analysis decisions in one easy-to-scan table.
  State each rule, why it was chosen and its source. Distinguish published
  methods from project-specific adaptations, and explain how these choices
  determine the reference, script behaviour and generated outputs.
- Provide explicit commands with expected inputs and outputs.
- Prefer practical, lab-runnable examples.
- Do not imply validation that was not performed.
- When documenting runs, include the pipeline version, params file, config file, samplesheet, output directory, and whether `-resume` was used.

## Git

- Before significant documentation changes or repository reorganization, create
  a local Git checkpoint commit of the current reviewed, tracked state. Do the
  same before deleting or moving any files under `/path/to/hpc-home/`, including
  files outside this repository. If the tracked working tree is already clean,
  record the existing HEAD as the checkpoint; an empty commit is unnecessary.
- Verify that the checkpoint commit succeeded before starting those changes.
  If Git is unavailable or committing is blocked, stop the dependent changes
  and explain the blocker. Do not claim an uncommitted snapshot is a checkpoint.
- Review `.gitignore` and the staged file list before every checkpoint. Keep
  source code, configuration, environment definitions, documentation and small
  provenance/approval manifests tracked. Exclude large sequencing/PAT/matrix
  payloads, generated result trees, runtime logs/caches, containers, credentials
  and nested downloaded software checkouts. Do not force-add ignored payloads.
  Check representative ignore rules with `git check-ignore`; remember that
  `.gitignore` does not remove files already tracked by Git.
- Git checkpoints protect tracked files only. Ignored data, external scratch
  targets and files outside this repository are not backed up by a commit.
  Before moving/deleting such files, record the exact scope and provenance and
  retain any required recovery evidence separately. A checkpoint does not
  replace the existing explicit deletion/overwrite approval requirements.
- Keep checkpoints local unless pushing/publishing is explicitly authorized.
  Account for Git object storage within the 1 TB home quota; do not duplicate
  large data into Git history as a preservation mechanism.
- Write concise imperative commit messages.
- Describe only what actually changed.
- Prefer a one-line subject under about 72 characters.
- Do not mention tests unless they were actually run.
- Avoid vague messages like "fix stuff" or "updates".

## Pull requests

- Write PRs for a technical reviewer.
- Include summary, why the change was needed, main files touched, validation performed with exact commands, and known limitations or follow-ups.
- Call out config, environment, parameter, reference, or sample-sheet changes that affect reproducibility.

## Philosophy

- Act as a cautious scientific collaborator.
- Keep validation proportional to concrete risks and observed issues. Reuse
  existing QC and successful checks; do not add experiments or repeat expensive
  validation solely for housekeeping when no specific concern warrants it.
  Retain required provenance, dependency and approved-deletion safeguards.
- Prefer transparency and parameter exposure.
- Avoid overengineering.
- Avoid confident biological claims without confirmation.

## Approved staged layout (2 November 2026)

- Project: `/path/to/hpc-home/analysis/CIRCALL_methylomics`; pipeline source:
  `/path/to/hpc-home/preprocessing/methylseq`. Project configuration pins source.
- Original `/path/to/hpc-home/em_seq` and `/path/to/hpc-home/methylseq_AH` are
  historical recovery locations; do not edit, move, or clean them during migration.
- Supported project components are in `workflows/`; four exploratory repositories
  live under `/path/to/hpc-home/brainstorming`. Keep scientific status and caveats.
- New execution state belongs to `runs/<workflow>/`; do not transplant old sessions.
- Source relocation does not validate cache reuse or historical resume.
- New outputs require exclusive atomic reservation, including rejection of empty
  existing directories and symlinks. Completed bundles cannot be modified through
  supported entrypoints. Use explicit ownership and input/source checks for resume.
- Fixtures/preflight only are authorized for this migration. Do not submit jobs,
  change scientific settings, copy bulk payloads, execute cleanup, or publish.
- User authorization covers creating the new layout and local commits. Existing
  recorded clean HEADs satisfy the original-state checkpoint requirement.
