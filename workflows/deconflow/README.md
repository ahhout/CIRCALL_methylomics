# DecoNFlow lineage integration

## Workspace transition and output ownership

This maintained workflow now belongs to `analysis/CIRCALL_methylomics`.
Use the original `/path/to/hpc-home/em_seq` entrypoints for historical resumes.
New runs reserve a previously nonexistent UTC-dated result directory; even an
empty directory or dangling symlink is a collision. Completed results are not
reused. Failed reservations and staging directories are retained for diagnosis.
The collection roots, references and scientific defaults have not changed.

See [the downstream migration audit](../../docs/reorganization/downstream_entrypoints.md)
for supported entrypoints, resume boundaries and fixture-only validation.

The selected comparison methods are UXM and CelFiE. `run_deconflow.sh` runs
native UXM by default and accepts an explicit `celfie` subcommand for patient
count inputs. Both routes use unchanged pinned upstream fitting code. EpiDISH
has been removed from the active execution route. Patient CelFiE requires
validated reference and patient counts before fitting; see
[the comparison record](ACTIVE_COMPARISON.md) for current acceptance status.
The reconstructed megakaryocyte component disagrees with the shipped UXM
profile at checked regions. Patient CELFIE outputs are exploratory and cannot
currently establish an exact shared-reference method comparison.

Current numerical findings and acceptance decisions are maintained in the
[comparison record](ACTIVE_COMPARISON.md), including the
[healthy-WBC diagnostics](ACTIVE_COMPARISON.md#healthy-wbc-diagnostic-check-27-october-2026)
and [scoring and convergence follow-up](ACTIVE_COMPARISON.md#scoring-and-convergence-repair-diagnostic-27-october-2026).
The [implementation record](IMPLEMENTATION_STATUS.md) retains dated execution
history; this README owns commands and input/output contracts.

The scientific question is how normal reference cell types contribute to bulk
EM-seq samples. The existing cohort contains twelve bone marrow (BM), peripheral
blood (PB), and plasma cell-free DNA (cfDNA) samples from four patients. The
completed UXM results remain available; this change does not require rerunning
them. PB is diagnostic material, not a healthy control. These estimates do not
provide validated blast fractions or establish tumour specificity.

## Scientific and reproducibility rules

| Rule | Reason | Source |
|---|---|---|
| DecoNFlow 2.5.0, commit `4195f65ad986ea5a68129d2da5b2b1effed25bc8`; unchanged upstream source and images | Avoid changes to native fitting code | `config/source.lock.json`, `config/images.lock.json`, `runtime.lock.json` |
| Original hg38 U25/l4 atlas: 979 regions and 40 cell types | Preserve the approved full reference | Official atlas; SHA256 `e58fbec4fe6200a311540c6e1a16be901a2d1e9ed9d74b6a57888f2feb635f51` |
| Indexed PAT inputs, containing methylation patterns along individual molecules | Reuse the shared PAT preparation workflow without BAM conversion here | [PAT workflow](../pat/README.md) |
| At least four informative CpGs; native l4 bins `0,0.251,0.75,1` and native boundary handling | Preserve the existing fragment classification | Pinned UXM image; no threshold override |
| Profiles selected explicitly in the plan | Preserve existing full and optional hematopoietic models | `config/plan.pats.json`; native UXM module |
| Hematopoietic profile: seven original cell columns and their 175 targeted markers | Preserve the conditional blood-cell model | Existing approved profile; `pat_inputs.py` |
| No EpiDISH fitting or EpiDISH matrix preparation in the active launcher | Apply the user's method-removal instruction | 26 October 2026 instruction; local orchestration change |
| CelFiE requires actual methylated/total counts; no invented depths or U-score conversion | Preserve the count model and reference provenance | [Comparison record](ACTIVE_COMPARISON.md) |
| CELFIE counts every called C/T CpG observation, weighted by PAT multiplicity, without UXM's l4 minimum | CpG observation counts and fragment-pattern classes are different inputs | `tools/pat_counts.py`; native counter validation |
| Preserve zero-depth entries and all 979 regions/40 cells | Avoid an undocumented coverage filter; native CelFiE applies its own pseudocounts | Unchanged pinned `bin/celfie.py`; [comparison record](ACTIVE_COMPARISON.md) |

The seven hematopoietic cell types are Blood-B, Blood-T, Blood-NK,
Blood-Mono+Macro, Blood-Granul, Eryth-prog and Megakaryocytes. That profile is
conditional on those cell types; the full atlas remains the primary cfDNA model.
UXM keeps its native per-sample coverage, missing-data handling and fitting
parameters. It computes its own intermediate counts. The separate count and
matrix preparation previously used for EpiDISH is no longer called.

## Inputs and execution

Use the existing `nextflow` Conda environment from the repository's
`environment.yml`. The launcher loads module-provided Apptainer 1.3.6.
The UXM reference, sample identifiers and native fitting options remain
unchanged. Count preparation uses existing system `g++` and `gzip`; the compiler
version and compiled source/binary checksums are recorded for each run. No
additional Conda package or image is installed.

Supply `--plan` and `--phase cohort`. The existing plan format remains unchanged:
`samples`, `profiles`, `cohorts`, and `fragment_rlen`. The samples TSV requires
`sample`, `cohort`, `pat`, and `csi`; IDs must be unique and each index must match
its PAT. Paths may be absolute or relative to the containing plan/TSV.
`fragment_rlen` remains 4. The historical cohort coverage thresholds are retained
in the plan for provenance and format compatibility; they do not filter UXM
markers and no longer produce EpiDISH matrices. Do not interpret them as UXM
coverage thresholds.

The existing [plan](config/plan.pats.json) and
[samplesheet](config/samples.pats.tsv) retain both historical profiles and all
twelve samples. Another explicit plan can select the intended existing profile
and sample set. Submit through Slurm:

```bash
w=/path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/deconflow
root=/path/to/hpc-scratch/tmp/em_seq/downstream/deconflow
runtime="$root/2026-10-05_14-36-deconflow-runtime-v2"
validation="$root/2026-10-05_15-11-deconflow-fixtures-v5"
analysis="$root/$(date -u +%Y-%m-%d_%H-%M)-deconflow-uxm-all12-v1"

sbatch --partition=shortq --time=02:00:00 --output="$root/uxm-%j.out" \
  "$w/run_deconflow.sh" "$analysis" --runtime "$runtime" \
  --validation "$validation" --plan "$w/config/plan.pats.json" --phase cohort
```

Outputs are native fractions at
`cohort/uxm/<profile>/deconvolution_results/reference_based/uxm/output.csv`,
with traces and exact parameters under the result bundle. `execution.json`
explicitly records `methods: ["uxm"]`, execution start and completion; each
attempt is recorded in `attempts.jsonl`. Completion requires the expected sample
and cell labels, finite bounded fractions and sums of one within output-rounding
precision. `COHORT_COMPLETE` refers to the requested UXM profiles only.

All tasks retain the shared work cache at
`/path/to/hpc-scratch/tmp/em_seq/work_emseq_integrated`.
PATs and indexes are linked. No automatic cleanup is enabled.

For an interrupted UXM run created by this version, retain its exact analysis path
and add `--resume`; keep its plan, inputs, code and configuration unchanged.
The original execution-start timestamp is retained. This edit does not repair
the audited run-local configuration verification issue. The interrupted
EpiDISH-matrix issue remains in the historical helper, which is no longer called
by this launcher. No recovery redesign is included.

## Patient CELFIE from matched counts

The local count helper streams each complete compressed PAT once. It counts
actual methylated and total CpG observations within the existing atlas regions
and separately feeds overlapping fragments to native wgbstools `homog` for
pooling checks. A fixed index look-back is avoided because it can miss long
spanning fragments. No genome-wide count matrix or large PAT copy is created.
The compressed input hash, gzip integrity, source identities, per-sample counts
and completion state are retained. Cache reuse checks source path, size and
mtime, settings, compiled binary and count-output checksums; it retains the
original compressed checksum without rereading cached raw PATs solely to
rehash them.

The reference pools published PAT specimens within each cell type and sums
integer methylated and total counts. Averaging percentages or inventing
effective depths is prohibited. The
[pooling evidence](ACTIVE_COMPARISON.md#matched-count-reference-and-patient-route-26-october-2026)
records specimen membership and the PAT-backed atlas checks. The separate
megakaryocyte count reconstruction remains subject to the
[reference-identity limitation](ACTIVE_COMPARISON.md#patient-results-and-interpretation-limit).

The reference gate must pass before a patient fit; a combined gate also binds
the patient matrix checksum. Preparation status and frozen commands are in
[`runs/2026-10-26_13-27-uxm-celfie-matched-reference-v1.json`](runs/2026-10-26_13-27-uxm-celfie-matched-reference-v1.json).
These gates check count integrity and PAT-pool agreement; they do not certify
megakaryocyte equivalence or biological accuracy. Directory existence does not
establish completion. The following command illustrates first-patient validation:

```bash
set -euo pipefail
w=/path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/deconflow
root=/path/to/hpc-scratch/tmp/em_seq/downstream/deconflow
reference=/path/to/hpc-scratch/tmp/2026-10-26_13-38-matched-reference-counts-v1
patients=/path/to/hpc-scratch/tmp/2026-10-26_13-39-matched-patient-counts-v1
analysis="$root/$(date -u +%Y-%m-%d_%H-%M)-deconflow-celfie-first-patient-v1"
sbatch --partition=shortq --time=06:00:00 --output="$root/celfie-%j.out" \
  "$w/run_deconflow.sh" celfie "$analysis" \
  --runtime "$root/2026-10-05_14-36-deconflow-runtime-v2" \
  --reference-counts "$reference/reference.counts.tsv" \
  --patient-counts "$patients/patient.counts.tsv" \
  --reference-validation "$reference/combined_validation_gate.json" \
  --samples SAMPLE_01 SAMPLE_03 SAMPLE_02 --diagnostics --verify-cache
```

This pilot selects PB, BM and cfDNA from the first existing patient. The count
tables retain 979 rows with exact atlas coordinates and paired
`<name>_METH`/`<name>_DEPTH` integers; the reference contains all 40 cell types in
atlas order. The adapter adds the reference coordinate triplet and lower-case
suffixes required by the native parser, preserving values and row order.

Fractions are written to `results/1_tissue_proportions.txt`. `VALIDATED.json`
and `COMPLETE.json` bind the checked output. `--verify-cache` requires native
tasks to be `CACHED` on recovery with identical fractions. `--diagnostics`
observes a replay of the same ten native starts and seeds; selected fractions
must match the ordinary native run exactly. This characterizes numerical
behavior, not biological accuracy. Native reference pseudocounts change `0/0`
to `1/2` internally; original zero-depth inputs remain recorded.

The route locks each run and verifies inputs, runtime, local source and frozen
configuration identities on resume. It uses the pinned local Python SIF and
retains original start times with separate attempt timestamps. Add `--resume`
to the exact original command to recover an interrupted run.

For an all-twelve candidate fit in a fresh bundle, retain the documented
sample identifiers and settings:

```bash
analysis="$root/$(date -u +%Y-%m-%d_%H-%M)-deconflow-celfie-all12-v1"
sbatch --partition=shortq --time=06:00:00 --output="$root/celfie-%j.out" \
  "$w/run_deconflow.sh" celfie "$analysis" \
  --runtime "$root/2026-10-05_14-36-deconflow-runtime-v2" \
  --reference-counts "$reference/reference.counts.tsv" \
  --patient-counts "$patients/patient.counts.tsv" \
  --reference-validation "$reference/combined_validation_gate.json" \
  --samples SAMPLE_01 SAMPLE_03 SAMPLE_02 \
            SAMPLE_04 SAMPLE_06 SAMPLE_05 \
            SAMPLE_07 SAMPLE_09 SAMPLE_08 \
            SAMPLE_11 SAMPLE_12 SAMPLE_10 --diagnostics --verify-cache
```

`tools/run_celfie_comparison.sh` exports all 480 signed differences in percentage
points, per-sample/per-cell summaries, stacked compositions and a difference
heatmap. It verifies the preserved UXM checksum and validated CELFIE fractions.
It does not refit UXM. For the completed cohort:

```bash
cohort="$root/2026-10-26_14-21-deconflow-celfie-all12-v1"
comparison="$root/$(date -u +%Y-%m-%d_%H-%M)-uxm-celfie-candidate-comparison-all12-v1"
sbatch --output="$root/uxm-celfie-comparison-%j.out" \
  "$w/tools/run_celfie_comparison.sh" --celfie "$cohort" \
  --fractions "$cohort/results/1_tissue_proportions.txt" --output "$comparison"
```

The [comparison record](ACTIVE_COMPARISON.md#patient-results-and-interpretation-limit)
owns the retained export selection, numerical findings and interpretation.
Reference incompatibility and unavailable patient ground truth prevent these
outputs from establishing a method accuracy ranking.

## Historical results and validation

Existing UXM/EpiDISH outputs, manifests, native parameters, images and locks are
preserved. The old no-plan/triplet launcher route is retired. To reproduce or
resume that historical integration, use its matching recorded code and runtime;
do not resume an old paired run with this changed launcher. Local checkpoint
`45aedc5941878fee90e8a589a26ce8c78e5e29ea` preserves the launcher before removal,
as well as the previously untracked CelFiE fixture and comparison records.
Historical acceptance evidence remains in
[IMPLEMENTATION_STATUS.md](IMPLEMENTATION_STATUS.md),
`config/pat_route_validation.json` and the original result bundles.
The retained historical validation utilities expect paired outputs and are not
acceptance tools for new UXM-only bundles.

The synthetic CelFiE launcher remains
[`tools/run_celfie_fixture.sh`](tools/run_celfie_fixture.sh). It does not submit
patient fits. Its prior acceptance does not validate patient count preparation
or biological accuracy. The historical method-removal change itself did not
execute a new pipeline; the subsequent patient runs are recorded separately.

Run local checks on an allocated compute node:

```bash
env -u SLURM_JOB_ID /path/to/hpc-home/.conda/envs/nextflow/bin/python \
  -m unittest discover -s "$w/tests" -v
bash -n "$w/run_deconflow.sh"
```
