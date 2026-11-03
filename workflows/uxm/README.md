# wgbstools + UXM lineage deconvolution

## Workspace transition and output ownership

This maintained workflow now belongs to `analysis/CIRCALL_methylomics`.
Use the original `/path/to/hpc-home/em_seq` entrypoints for historical resumes.
New runs reserve a previously nonexistent UTC-dated result directory; even an
empty directory or dangling symlink is a collision. Completed results are not
reused. Failed reservations and staging directories are retained for diagnosis.
The collection roots, references and scientific defaults have not changed.

See [the downstream migration audit](../../docs/reorganization/downstream_entrypoints.md)
for supported entrypoints, resume boundaries and fixture-only validation.

## Open results

**All 12 samples:** [full-atlas PDF](../epidish/exports/2026-09-29_10-25-deconvolution-all12-display-v2/uxm/uxm_full.pdf),
[haematopoietic PDF](../epidish/exports/2026-09-29_10-25-deconvolution-all12-display-v2/uxm/uxm_hematopoietic.pdf),
and [UXM-only tables](../epidish/exports/2026-09-29_10-25-deconvolution-all12-display-v2/uxm/).
These display the saved native UXM estimates from both cohorts without refitting.
[EpiDISH and cross-method comparisons](../epidish/RESULTS.md) have separate folders.


| Cohort | Full atlas | Haematopoietic atlas |
| --- | --- | --- |
| Eight-sample baseline (August 31) | [Results](results/2026-09-10_11-19-emseq-integrated-8samples-trim-8-3-3-10-v1-uxm-full-atlas/results/) | [Results](results/2026-09-10_12-33-emseq-integrated-8samples-trim-8-3-3-10-v1-uxm-hematopoietic/results/) |
| Newest cohort: four new HEH samples (September 16) | [Results](results/2026-09-24_20-57-emseq-integrated-heh-4new-trim-8-3-3-10-v1-uxm-full-atlas/results/) | [Results](results/2026-09-24_22-38-emseq-integrated-heh-4new-trim-8-3-3-10-v1-uxm-hematopoietic/results/) |

Each directory contains `lineage_fractions.pdf` and `lineage_fractions.csv`.
For the corresponding EpiDISH plots, use the [EpiDISH output guide](../epidish/RESULTS.md).
These are separate cohorts; see the [complete results index](../../RESULTS.md).

## Workflow

This is the canonical lineage branch downstream of `nf-core/methylseq`. Two
Slurm launchers provide a full-atlas analysis and a hematopoietic-focused
analysis. Both call the same validated implementation, which validates the
selected methylseq BAMs, verifies and links an explicit shared official-wgbstools
hg38 PAT/LBETA bundle, and runs the pinned official UXM deconvolution and plotting
commands.

The workflow discovers BAM samples from the selected methylseq result and
orders them according to `EMSEQ_METADATA` (default `metadata/metadata.tsv`);
sample IDs are not hard-coded. Select a cohort-matched metadata snapshot: the
current production sheet has four samples, while the eight-sample historical
cohort has its own frozen inputs.
It does not generate CpG matrices, PCA, DMRs, or a composition-adjusted model.

## How UXM estimates cell-type fractions

UXM does not assign each individual fragment to one cell type. Instead, it
estimates the mixture of reference cell types that best explains the sample
across many cell-type-informative regions:

1. The published U25 hg38 atlas contains regions selected because one normal
   cell type is specifically unmethylated there. Its cell-type columns record
   the expected fraction of unmethylated fragments at every region.
2. A PAT file preserves the consecutive CpG methylation pattern observed on
   each sequenced fragment. For this `l4` atlas, UXM considers fragments with
   at least four informative CpGs and uses the published U-score definition:
   fragments with no more than 25% methylated CpGs are classified as
   unmethylated (`U`); intermediate and highly methylated patterns are retained
   in the underlying `U/X/M` summary.
3. For each atlas region, UXM calculates the sample's fraction of eligible
   fragments classified as `U`. Regions are weighted by their contributing
   fragment counts, so better-covered observations have greater influence.
4. UXM uses non-negative least squares to find the combination of atlas
   cell-type profiles that best reconstructs the sample's regional U scores,
   then normalizes the fitted coefficients to sum to one. These normalized
   coefficients are the reported cell-type fractions.

A lineage absent from the normal atlas, including a BCP-ALL blast population,
cannot be reported directly: its signal may be distributed among the most
similar reference profiles because the fitted fractions are normalized to
100%. Laboratory blast percentages must remain an independent measurement.

## Analysis profiles

`run_lineage_deconvolution.sh` uses all 40 normal cell types in the published
atlas. This is the reference profile for cfDNA and also shows when a sample is
poorly represented by the atlas.

`run_hematopoietic_deconvolution.sh` retains only `Blood-B`, `Blood-T`,
`Blood-NK`, `Blood-Mono+Macro`, `Blood-Granul`, `Eryth-prog`, and
`Megakaryocytes`, including only the markers targeted to those seven types.
Use this as a secondary, conditional normal-haematopoietic fit for the current
cellular PB and Ficoll-enriched BM samples. It is not a direct blast-fraction
measurement.

The implementation always appends the selected profile to the supplied run
base name. For example, `ANALYSIS_ID` becomes
`ANALYSIS_ID-uxm-full-atlas` or `ANALYSIS_ID-uxm-hematopoietic`. The Slurm job
and output names also identify the profile.

## Accepted scope and PAT dependency

Full/heme E31 runs passed technical validation. Evidence and pinned dependencies
are in the [release record](../../docs/audits/downstream_structure_20260914/releases.md).
Prepare or validate the shared bundle using the [PAT workflow](../pat/README.md).
PAT conversion parameters and validation instructions are owned there.
The pinned UXM installation resides in the ignored `software/` directory.
Its checkout and contents were preserved during physical relocation. The repository `nextflow` environment
supplies the existing Python and wgbstools dependencies.

## Run UXM after bundle completion

Both profiles require `METHYLSEQ_RESULTS_DIR`, `PAT_INPUT_DIR` (the bundle
root containing the manifest), and a new `UXM_RUN_BASE_NAME`. For example:

```bash
EMSEQ_METADATA=/path/to/hpc-scratch/nextflow/nf-core-methylseq/results/uxm_deconvolution/2026-09-11_10-07-emseq-integrated-8samples-trim-8-3-3-10-v1-uxm-bcpall-matched907-v0.1/submission/metadata.tsv \
METHYLSEQ_RESULTS_DIR=/path/to/hpc-scratch/nextflow/nf-core-methylseq/results/bulk/2026-08-31-emseq-integrated-8samples-trim-8-3-3-10-v1 \
PAT_INPUT_DIR=/path/to/hpc-scratch/nextflow/nf-core-methylseq/results/uxm_deconvolution/2026-08-31-emseq-integrated-8samples-trim-8-3-3-10-v1-pat-v1 \
UXM_RUN_BASE_NAME=2026-09-29_14-35-cohort-analysis-v1 \
  sbatch /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/uxm/run_lineage_deconvolution.sh
```

Use the same variables with `run_hematopoietic_deconvolution.sh` for the
seven-cell profile. Use a new name; the recorded E31 profiles already exist. Both profiles link the same validated PATs; `CONVERSION_MODE=rebuild`
is refused in analysis launchers. Build PATs once with the [PAT launcher](../pat/README.md).

Following the methylseq layout, incomplete work and final outputs have separate
locations:

- Working directories and intermediates:
  `/path/to/hpc-scratch/tmp/em_seq/downstream/lineage_deconvolution/`
- Completed UXM run directories:
  `/path/to/hpc-scratch/nextflow/nf-core-methylseq/results/uxm_deconvolution/`
- Repository entry point: [`results/`](results/), a symbolic link to the
  completed-results directory.

Each completed run contains `results/lineage_fractions.csv`,
`results/lineage_fractions.pdf`, its PAT manifest and provenance. PAT artifacts
may themselves be symbolic links to authoritative inputs, avoiding large data
duplication. A run is published to the completed-results directory only after
UXM output validation succeeds. A persistent `<final-directory>.lock` is held
from before work reservation through publication. Existing files, directories
and dangling symlinks are refused. Publication uses `mv -nT` and verifies that
the work source disappeared; a late collision fails and retains the work.
`execution.json` records execution `started_utc` before analysis, plus status
and `ended_utc` on exit; `completed_utc` appears only on success. An incomplete
standalone UXM run requires a fresh dated run name; there is no automatic resume.

The official atlas represents normal human lineages; its output is not a
validated BCP-ALL blast fraction. Laboratory blast percentages remain separate
metadata.

## Experimental BCP-ALL reference

The normal atlas cannot directly report a validated BCP-ALL blast fraction.
An independently trained disease-state profile is therefore being evaluated
as a separate, gated development project; the eight project samples are not
training data.

Scientific rationale, design, validation criteria and stop conditions live
only in
[`BCP_ALL_ATLAS_PROJECT_RATIONALE.md`](/path/to/hpc-home/brainstorming/bcpall_atlas/BCP_ALL_ATLAS_PROJECT_RATIONALE.md).
Current execution status and commands live only in the
[`bcp_all_atlas_development` README](/path/to/hpc-home/brainstorming/bcpall_atlas/README.md).
