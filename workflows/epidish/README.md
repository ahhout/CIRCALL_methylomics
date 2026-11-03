# UXM atlas + EpiDISH workflow

## Workspace transition and output ownership

This maintained workflow now belongs to `analysis/CIRCALL_methylomics`.
Use the original `/path/to/hpc-home/em_seq` entrypoints for historical resumes.
New runs reserve a previously nonexistent UTC-dated result directory; even an
empty directory or dangling symlink is a collision. Completed results are not
reused. Failed reservations and staging directories are retained for diagnosis.
The collection roots, references and scientific defaults have not changed.

See [the downstream migration audit](../../docs/reorganization/downstream_entrypoints.md)
for supported entrypoints, resume boundaries and fixture-only validation.

**[Open all-12-sample EpiDISH and UXM displays, with method comparisons kept separate](RESULTS.md).**

This is a user-accepted workflow, and the user's preferred method for blastic
samples. It preserves RPC, matched UXM and unweighted NNLS comparison/QC outputs.
The official atlas is shared with the supported UXM workflow. The experimental
BCP-ALL atlas is not selected implicitly.

## Run the accepted analysis

Requires completed full/heme [UXM runs](../uxm/README.md), both derived from the
same validated [PAT bundle](../pat/README.md), and their methylseq source.
From the repository root:

```bash
METHYLSEQ_RESULTS_DIR=/absolute/path/to/methylseq-result \
UXM_FULL_RUN_DIR=/absolute/path/to/completed-full-uxm \
UXM_HEMATOPOIETIC_RUN_DIR=/absolute/path/to/completed-heme-uxm \
EPIDISHR_RUN_NAME=2026-09-29_14-35-epidish-cohort-v1 \
  sbatch /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/epidish/run_epidish_comparison.sh
```

New runs use regional coverage thresholds 1, 3 and 5, as approved on
2026-09-14. Historical 1/5/10 results remain unchanged; validation and plot
refresh read the thresholds recorded in each run.

The workflow default RPC iteration limit is 75 (approved 2026-09-14). This command retains the existing full/heme
UXM input requirement; it is not a new standalone fitting route. Expected
outputs, recovery and validation commands are documented below. Accepted E31
execution and its eight nonconvergence warnings are retained in the
[release evidence](../../docs/audits/downstream_structure_20260914/releases.md).

This workflow compares EpiDISH Robust Partial Correlations (RPC), unweighted
NNLS, and the existing coverage-weighted UXM NNLS implementation. It uses the
same canonical hg38 `Atlas.U25.l4` reference and the same cohort-wide marker
set for each matched comparison. Native UXM results are included as context
because their default marker handling differs.

The workflow does **not** estimate a validated BCP-ALL blast fraction. The
atlas contains normal reference cell types, so disease-state DNA can be
assigned to the closest normal profiles. Results measure agreement and
coverage-threshold robustness between tools, not absolute accuracy.

The unweighted NNLS implementation uses `scipy.optimize.nnls`, gives every
retained marker equal influence, and normalizes the non-negative coefficients
to sum to one. This matches UXM's solver and normalization while omitting
UXM's fragment-count weighting.

## Native versus matched UXM

Both outputs use UXM's coverage-weighted NNLS algorithm. They differ only in
marker handling:

- **Native UXM** is imported from the original UXM run. UXM drops markers
  without usable coverage separately for each sample, so individual samples
  can be fitted using different marker sets.
- **Matched UXM** uses a fixed cohort-wide marker mask. A marker must have at
  least 1, 3, or 5 fragments in every included sample, and exactly the same
  marker rows are supplied to EpiDISH, unweighted NNLS, and weighted UXM.

Use EpiDISH versus matched UXM to compare fitting methods, weighted versus
unweighted NNLS to isolate the effect of fragment-count weighting, and native
versus matched UXM to assess the effect of marker filtering. Native UXM is
contextual and is not the controlled cross-method comparison.

## Runtime

Create the isolated, pinned R environment once:

```bash
conda env create --file /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/epidish/environment.yml
```

The environment pins R 4.5 and EpiDISH 2.26.0. The launcher switches to the
existing `nextflow` environment when it invokes the pinned UXM checkout and
module-independent wgbstools executable.

Run the synthetic validation and temporary recovery fixtures before the first analysis:

```bash
conda run --no-capture-output -n epidishr \
  Rscript /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/epidish/scripts/test_synthetic.R
conda run --no-capture-output -n epidishr \
  python3 /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/epidish/scripts/test_recovery_contract.py
```

## Evaluation status — 2026-09-14

The three completed `*-pilot-integrated-epidish-comparison*` results are
four-sample E20 comparisons (baseline, RPC-200 and RPC-200 with NNLS), retained
as method-evaluation evidence. Both E31 UXM profiles completed. E31 RPC-50
comparison retry 5685252 replaces job 5684633, which failed a node-local
filesystem device-ID check before analysis; the validator has been corrected.
Retry 5685252 completed and independent validation job 5685170 passed. Eight
RPC-50 non-convergence warnings remain recorded; structural acceptance does
not establish biological accuracy or convergence of those fits. See the
[migration record](../../docs/audits/e31_migration_20260910/README.md).

The historical eight-sample E31 sheet is a recorded cohort snapshot; it differs
from the current production sheet. Use explicit `EMSEQ_SAMPLESHEET` and
`EMSEQ_METADATA` snapshots when supplying historical UXM cohorts. New submissions
require `METHYLSEQ_RESULTS_DIR` in addition to the two UXM run directories and a
new valid UTC output name. The launcher validates both complete
cohort manifests against the selected `EMSEQ_SAMPLESHEET` (the canonical
`samplesheets/samplesheet.csv` is the default for matching current cohorts), their source BAMs and
the same resolved PAT/CSI/LBETA payloads before creating comparison work.
Historical outputs remain immutable; the legacy command below is a record,
not a command to rerun into an existing directory with the current launcher.

For new completed results, run the small input-driven validator directly:

```bash
conda run --no-capture-output -n epidishr \
  Rscript workflows/epidish/scripts/validate_comparison.R \
  /absolute/path/to/completed/comparison /absolute/path/to/matching-samplesheet.csv
```

It derives marker counts from the run records and cross-checks masks, matrices,
coverage thresholds, sample order, fractions and checksums. It does not impose
four-sample marker counts on E31. Its read-only check passed on the completed
RPC-200/NNLS historical result without the optional current-cohort samplesheet.
The four-sample validator was archived on 2026-09-09 with an exact helper
snapshot: [historical acceptance validator](../../archive/downstream/epidish_four_sample_validation_20260908/README.md).
Use `scripts/validate_comparison.R` for current comparisons.

## Historical runs and sensitivity settings

Historical four-sample commands and manifests remain in `runs/` and the frozen
[reconciliation record](../../docs/audits/reconciliation_status_20260914.md).
They are execution evidence, not current submission examples.

This workflow defaults to a maximum of 75 RPC robust-regression iterations.
Historical RPC-50 results retain their original settings. Sensitivity runs can
set a different positive integer explicitly, for example:

```bash
METHYLSEQ_RESULTS_DIR=/absolute/path/to/source-methylseq-run \
EPIDISHR_RPC_MAXIT=200 \
EPIDISHR_RUN_NAME=2026-10-28_15-30-epidish-rpc-maxit200-v1 \
UXM_FULL_RUN_DIR=/absolute/path/to/full-uxm-run \
UXM_HEMATOPOIETIC_RUN_DIR=/absolute/path/to/hematopoietic-uxm-run \
  sbatch /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/epidish/run_epidish_comparison.sh
```

The selected limit is recorded in `provenance/input_manifest.json`.
Because the matched UXM estimates do not depend on EpiDISH's iteration limit,
an iteration-limit sensitivity run can reuse their validated fraction tables:

```bash
EPIDISHR_MATCHED_BASELINE_DIR=/absolute/path/to/completed-epidish-result \
METHYLSEQ_RESULTS_DIR=/absolute/path/to/source-methylseq-run \
EPIDISHR_RPC_MAXIT=200 \
EPIDISHR_RUN_NAME=2026-10-28_15-30-epidish-rpc-maxit200-v1 \
UXM_FULL_RUN_DIR=/absolute/path/to/full-uxm-run \
UXM_HEMATOPOIETIC_RUN_DIR=/absolute/path/to/hematopoietic-uxm-run \
  sbatch /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/epidish/run_epidish_comparison.sh
```

Baseline reuse requires a completed run whose input and comparison manifests
verify every requested fraction table, marker mask, reference matrix and regional
sample table. Ordered samples and cell types, atlas checksum, UXM commit,
regional-measure definition, profile, coverage threshold and recorded PAT/CSI/LBETA
hashes and sizes must agree. Marker masks and reference/sample tables must match
exactly, comparing decompressed gzip contents. Relocated paths and filesystem
timestamps do not determine compatibility. The RPC iteration limit may differ in
this new sensitivity run because it does not enter the reused UXM calculation.
All six comparisons are checked before any estimate is imported. Missing historical
evidence requires recomputation in a fresh run. Reuse evidence is saved in
`provenance/matched_uxm_reuse.json`.

New runs freeze `provenance/recovery_contract.json` before preparation writes
outputs. It records scientific inputs, metadata, samplesheet, configuration,
script/environment identities, RPC iteration limit and baseline selection.
Recovery requires exact agreement with that immutable contract and valid completion
receipts for every existing prepared stage and matched estimate. Use the original
launcher revision and unchanged inputs when recovering; a revised launcher requires
a fresh run directory. Legacy incomplete runs without the contract cannot resume.
Do not reconstruct their missing timestamps or adopt unreceipted partial tables.

Preparation and matched/NNLS estimates have receipts binding their input fingerprint
to output checksums. Small raw matched CSVs and their receipts remain available
through finalization. Derived comparisons, summaries and figures are generated in
separate attempt staging; interrupted generation preserves reusable estimates.
Partially published finalization outputs without a finalization receipt block
recovery. Previous staging is retained rather than overwritten or cleaned up.

For a recoverable new run, repeat the original submission with the failed job ID
and the same RPC limit and baseline selection, if set:

```bash
METHYLSEQ_RESULTS_DIR=/absolute/path/to/source-methylseq-run \
EPIDISHR_RECOVER_JOB_ID=123456 \
EPIDISHR_RPC_MAXIT=200 \
EPIDISHR_MATCHED_BASELINE_DIR=/absolute/path/to/completed-epidish-result \
EPIDISHR_RUN_NAME=2026-10-28_15-30-epidish-rpc-maxit200-v1 \
UXM_FULL_RUN_DIR=/absolute/path/to/full-uxm-run \
UXM_HEMATOPOIETIC_RUN_DIR=/absolute/path/to/hematopoietic-uxm-run \
  sbatch /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/epidish/run_epidish_comparison.sh
```

For the historical four-sample pilot only, the original acceptance command remains:

```bash
conda run --no-capture-output -n epidishr \
  Rscript /path/to/hpc-home/analysis/CIRCALL_methylomics/archive/downstream/epidish_four_sample_validation_20260908/scripts/validate_pilot.R \
  /path/to/hpc-scratch/nextflow/nf-core-methylseq/results/epidishr/pilot_integrated_20260720-epidish-comparison
```

## Inputs

The full UXM run supplies the ordered PAT manifest, PAT/CSI files, and cached
regional `U/X/M` counts. The haematopoietic UXM run supplies its native
seven-cell result. Both manifests must contain exactly the same samples in the
same order.

For each marker and sample, the workflow calculates:

- `U / (U + X + M)` when the atlas direction is `U`;
- `M / (U + X + M)` when the atlas direction is `M`.

The full profile retains all 40 atlas columns. The haematopoietic profile
retains `Blood-B`, `Blood-T`, `Blood-NK`, `Blood-Mono+Macro`, `Blood-Granul`,
`Eryth-prog`, and `Megakaryocytes`, together with markers targeting those
types. Rows with an invalid reference value are removed. Separate common
marker masks then require at least 1, 3, or 5 fragments in every included
sample.

## Outputs

Incomplete work is kept under:

```text
/path/to/hpc-scratch/tmp/em_seq/downstream/epidishr/
```

A stable lock beside the final destination is held from before work reservation
through publication. Existing files, directories and dangling symlinks are refused.
Publication uses `mv -nT` and verifies that the source disappeared; a late collision
fails and preserves work. Lock files remain in place so concurrent submissions use
the same lock inode.

`provenance/execution_attempts.json` retains the original `started_utc` and appends
one record for each distinct Slurm job ID. Entry timestamps describe execution,
separately from the timestamp in the run name. Successful attempts record
`completed_utc`; failed attempts record `ended_utc` and failure status. No historical
runs are backfilled. The existing `EPIDISHR_COMPLETED` marker continues to mean
successful analysis finalization; attempt timing also tracks final publication.

Completed results are published under:

```text
/path/to/hpc-scratch/nextflow/nf-core-methylseq/results/epidishr/<run-name>/
```

The repository `results` link points to that completed-results directory.
Each matched comparison is self-contained:

```text
comparisons/
├── full/
│   ├── coverage-ge1/
│   ├── coverage-ge3/
│   └── coverage-ge5/
└── hematopoietic/
    ├── coverage-ge1/
    ├── coverage-ge3/
    └── coverage-ge5/
```

Each folder contains its marker mask, coverage summary, regional sample and
reference matrices, EpiDISH fractions and captured warnings, unweighted NNLS
fractions, matched and native UXM fractions, comparison metrics, plots,
filtered UXM atlas, and manifest. The top-level
`comparisons/summary/` contains the cross-method tables, metrics and threshold
sensitivity report. `epidish_r/` contains only EpiDISH RPC plots, fractions,
warnings and marker counts. New runs no longer generate the redundant
`composition_overview.pdf`. Historical completed bundles retain their original
layout and manifests; readers support both layouts. `provenance/` records inputs,
checksums, commands, software versions, sample metadata, and R session
information. `provenance/receipts/` holds preparation, NNLS, matched-estimate,
auxiliary-provenance and finalization receipts. Run-relative table and receipt
paths remain valid after publication.

No BAM, PAT, or LBETA data are copied. Source data and existing UXM results
are read in place and are not modified.

## Standalone EpiDISH output

New runs write `epidish_r/epidish_full.pdf` and
`epidish_r/epidish_hematopoietic.pdf`, with PNG counterparts and
`epidish_annotated_fractions.tsv`. These figures show EpiDISH RPC only,
with patients in separate panels and recorded coverage thresholds in rows.
Samples are ordered BM, PB, then cfDNA within each patient. Each sample is
labelled with its compartment, sample ID and BM/PB tumor burden (from the
recorded blast percentage). Hematopoietic cell types use contrasting colors;
other reference types use muted gray tones and remain separate. Plasma tumor burden is shown as not recorded. Recorded
blasts are metadata annotations, not EpiDISH estimates. Asterisks flag fits
with captured warnings. All reference cell types remain unmerged.

Regenerate figures from saved estimates; no EpiDISH or UXM fitting is run:

```bash
conda run --no-capture-output -n epidishr Rscript \
  workflows/epidish/scripts/plot_epidish_standalone.R \
  /absolute/path/to/completed/comparison
```

By default, outputs go to `figures/<UTC-timestamp>/` inside the completed run,
for example `figures/2026-09-29_14-35-epidish-figures-v1/` (UTC). Each folder contains both PDFs, both
PNGs, `epidish_annotated_fractions.tsv` and `epidish_plot_provenance.json`.
An optional second argument selects a different new figure directory.
Append `--coverage 1` to show only the saved >=1 estimates; omitting it
shows all saved thresholds. Within each bar, cell types are stacked from
largest at the baseline to smallest at the top, retaining consistent colors.
The original `summary/`, comparison outputs and run manifests remain unchanged.
The command uses the run's saved metadata and thresholds, records input
checksums, and refuses to overwrite existing standalone outputs. It performs
no fitting. The September 14 preview uses the completed eight-sample 1/5/10
run; the 1/3/5 analysis and its validation were cancelled as requested
(jobs 5717845 and 5717865).

## Display all samples from existing runs

`results/` links to the original fitted analyses on scratch. `exports/` holds
timestamped presentation bundles made from their saved estimates. The old
`outputs/` path is a compatibility link to `exports/` for recorded commands.

Use `scripts/export_saved_deconvolution.R` to combine saved estimates for
presentation. Supply an explicit metadata file defining the expected sample
set and explicit completed source runs; the script refuses missing, duplicate
or overlapping cohorts. It checks saved fraction checksums, sample grids and
native UXM inputs before writing. It does not fit, filter markers or merge
samples. The output folder must be new and use a valid UTC timestamp.

```bash
conda run --no-capture-output -n epidishr Rscript \
  workflows/epidish/scripts/export_saved_deconvolution.R \
  "workflows/epidish/exports/$(date -u +%Y-%m-%d_%H-%M)-deconvolution-all12-display-v1" \
  metadata/metadata.tsv 1 \
  workflows/epidish/results/2026-09-10_15-01-emseq-integrated-8samples-trim-8-3-3-10-v1-epidish-comparison \
  workflows/epidish/results/2026-09-25_00-04-emseq-integrated-heh-4new-trim-8-3-3-10-v1-epidish-comparison
```

Run this small plotting command directly on the current node. The `1` selects the
saved coverage >=1 estimates for the EpiDISH figures; tables keep all source
thresholds. Native UXM keeps its original per-sample marker handling (the
export's numeric coverage value `0` is a display sentinel, not a filtering
threshold). No methylseq or deconvolution pipeline is run.

| Display rule | Reason | Source |
| --- | --- | --- |
| Expected samples equal the union of disjoint source cohorts | Catch missing samples and duplicate exports | Supplied metadata and saved manifests |
| Keep source masks, thresholds and RPC limits | Preserve fitted estimates | Each source input manifest |
| Coverage >=1 in current combined figures | Common saved threshold; consistent with the prior baseline display | Explicit export argument; presentation only |
| Keep saved annotations; record current metadata differences | Avoid silently changing historical labels | Saved and current metadata |
| Separate `epidish_r/`, `uxm/`, and `comparisons/` | Keep method-only results separate from cross-method QC | User request, 29 September 2026 |

`display_manifest.json` records execution start/completion, source settings,
input checksums, scripts and command. `sample_coverage.tsv` lists every sample.
Each fraction retains its source run. The comparison folders link to original
per-cohort comparisons and selected summary files; no large inputs are copied.
Original completed results and manifests remain unchanged.
