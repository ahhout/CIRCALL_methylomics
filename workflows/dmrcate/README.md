# Paired bulk EM-seq with DMRcate

## Workspace transition and output ownership

This maintained workflow now belongs to `analysis/CIRCALL_methylomics`.
Use the original `/path/to/hpc-home/em_seq` entrypoints for historical resumes.
New runs reserve a previously nonexistent UTC-dated result directory; even an
empty directory or dangling symlink is a collision. Completed results are not
reused. Failed reservations and staging directories are retained for diagnosis.
The collection roots, references and scientific defaults have not changed.

See [the downstream migration audit](../../docs/reorganization/downstream_entrypoints.md)
for supported entrypoints, resume boundaries and fixture-only validation.

**Undergoing validation; not yet the accepted replacement for DMRichR.**
This workflow asks where bulk methylation differs between bone marrow (BM),
peripheral blood (PB) and plasma cell-free DNA (cfDNA), within four patients
with diagnostic hyperdiploid B-cell precursor acute lymphoblastic leukaemia.
All twelve samples are retained; each comparison uses eight samples and models
patient pairing. PB is diagnostic material, not a healthy control. A compartment
difference can reflect cell mixture and does not establish tumour specificity.

**On hold as of 27 October 2026:** comparing patient plasma with healthy
plasma and validating plasma-marker disease specificity, until suitable
healthy plasma cfDNA controls are available. The user does not expect these
controls soon; no availability date or dependent analysis is assumed.
Patient compartment comparisons and candidate exploration can continue.
See the [project control policy](../../PROJECT.md#controls-and-reference-material).

Public healthy white-blood-cell (WBC) profiles remain a separate exploratory
option: they can show whether a candidate region such as CDH4 is methylated
in ordinary healthy blood-cell DNA. They are not plasma controls or additional
paired samples for the current model. Their age/sex distribution, WGBS
chemistry and external processing limit interpretation of differences from
this EM-seq cohort. Any descriptive comparison first requires validated data;
this documentation change does not initiate that analysis.

The user selected minimal analysis and essential tables on 6 October 2026.
The workflow uses unmodified DMRcate for inference, bsseq for Bismark import
and methylation fractions, and edgeR for the methylation-channel design matrix.
Project code supplies sample/reference checks, coverage eligibility, Slurm
execution, compact gene annotations, paired effect tables and provenance.
No bespoke HTML/PDF report, exploratory PCA/clustering/correlation plots,
methrix wrapper, promoter/island/shore/shelf classification or report-only
recovery machinery is retained. Historical runs and their reports remain intact.

## Healthy-WBC candidate comparison — 27 October 2026

The user authorized a descriptive check of the two retained cfDNA DMRs against
three public healthy adult female white-blood-cell samples (GSM6810026–28,
ages 37, 66 and 74). This follow-up uses the exact original 3 and 12 CpGs;
it does not change genome-wide discovery, pool controls into the paired model,
or perform new significance tests. Healthy-plasma specificity work remains
on hold.

| Region, hg38 | Healthy WBC, n=3 | Patient PB, n=4 | Patient BM, n=4 | Patient cfDNA, n=4 |
| --- | --- | --- | --- | --- |
| MTND2P28, chr1:629824–629882 | 1.3–8.3% | 1.3–4.1% | 1.5–4.9% | 28.9–67.5% |
| CDH4, chr20:61252799–61252885 | 0–10.3% | 0–63.4% | 51.8–93.8% | 45.7–78.4% |

Values are per-sample arithmetic means of CpG methylation fractions, not
transformed model coefficients. Every predefined CpG was observed in every
sample and had at least 3× coverage. No extra coverage filter or imputation
was applied. Missing CpGs would produce an NA regional mean rather than a
changed target set. The healthy PAT pass counts all overlapping patterns,
including fragments starting before a target, and respects multiplicities.

CDH4 merits follow-up because all patient marrow and plasma values exceed
these three healthy-WBC values. This is an external reference observation,
not proof of leukaemia specificity: sample types, age/sex distributions,
WGBS versus EM-seq chemistry and processing differ. At the chr1 pseudogene,
healthy WBC resembles patient cellular DNA and also has high average depth
(about 1,240–1,490 observations per CpG); mapping/mitochondrial context remains
unresolved. The result alone cannot establish a disease-specific plasma signal.

[Figure (PDF)](/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-27_08-50-healthy-wbc-candidate-comparison-v1/candidate_methylation.pdf),
[region table](/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-27_08-50-healthy-wbc-candidate-comparison-v1/region_summary.tsv)
and per-CpG counts are retained with a timestamped manifest in the same bundle.
The [run record](runs/2026-10-27_08-50-healthy-wbc-candidate-comparison-v1.json)
records source/input hashes and checks. All three controls passed full-stream
gzip integrity checking and matched their staged SHA-256 hashes; the original
staging record was preserved. Extracted patient means reproduce the original
DMRcate output within 1e-12. Two overlap/integrity tests passed. Execution was
directly on node n09, taking 3m13.82s with about 1.5 GiB peak RSS.

Reproduce in a new UTC-named bundle, using the existing `nextflow` environment:

```bash
w=/path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/dmrcate
root=/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate
module load containers/apptainer/1.3.6
/path/to/hpc-home/.conda/envs/nextflow/bin/python "$w/compare_healthy_wbc.py" \
  --patient-run "$root/2026-10-26_12-48-dmrcate-cohort-all12-v3" \
  --controls /path/to/hpc-home/data/literature_controls/2026-10-27_08-30-healthy-wbc-fox-fisher-v1 \
  --output "$root/$(date -u +%Y-%m-%d_%H-%M)-healthy-wbc-candidate-comparison-v2" \
  --runtime "$root/2026-09-29_14-28-runtime-bioc322-v2" \
  --tabix /path/to/hpc-home/.conda/envs/nextflow/bin/tabix \
  --apptainer "$(command -v apptainer)" --threads 3
```

## Analysis decisions

| Rule | Reason | Source |
|---|---|---|
| BM−PB, cfDNA−BM, cfDNA−PB; four paired patients each | Three independent compartment questions | Approved project design |
| Explicit twelve source reports from August 31 and September 16 | Preserve original IDs and actual input provenance | `config/samples.tsv`, upstream parameter records |
| Existing hg38 + lambda/pUC19 sequence and CpG namespace | Preserve alignment compatibility | Existing reference manifest and FASTA SHA256 |
| bsseq strand collapsing; integer methylated/total counts in HDF5 | Represent each CpG dyad once without duplicating large matrices | bsseq reader; approved project design |
| Validate dictionary-listed chromosomes against the existing CpG namespace; remaining contigs against the pinned combined FASTA | Accept human alternate/unplaced sequence, chrM and controls present in whole-genome reports; reject unknown or invalid loci | Existing reference and Biostrings FASTA reader |
| chr1–22, X, Y; ≥3 reads in all eight samples per contrast | Independent complete tested universes; no sex exception for Y | Approved project design |
| `~ patient + compartment`, expanded by `edgeR::modelMatrixMeth` | Model paired methylated/unmethylated counts | Approved design; edgeR/DMRcate sequencing interface |
| CpG BH FDR 0.05; `pcutoff="fdr"`; lambda 1000; C 2; min.cpgs 2 | Preserve agreed genomic-distance DMR discovery | Approved design; DMRcate parameters |
| No effect, concordance, SNP, target or upper-depth filter | Preserve the discovery universe | Approved project design |
| Mean CpG fraction per sample/region, four paired differences and their mean | Report percentage-point effects separately from transformed coefficients | Existing project reporting definition |
| Gene overlaps and nearest transcript TSS from frozen GENCODE v50 basic | Essential gene context using the existing annotation release | GENCODE v50; minimal output scope selected 6 October |
| Coverage retention at 1×, 3×, 5× and 10× | Describe coverage support without changing the ≥3× inference threshold | Existing QC definition |
| Three separate exploratory comparisons; no pooled FDR claim | Each contrast has its own tested universe | Approved interpretation |

A transcription start site (TSS) is the position where a transcript begins.
Zero methylated reads at positive coverage is an observed zero; zero coverage
remains missing. Upstream lambda/pUC19 control QC is retained unchanged.
Saved UXM/EpiDISH estimates remain separate descriptive context in
[the results index](../../RESULTS.md); they are not an analysis dependency
or an adjustment in the statistical model. Enrichment and marker selection
remain outside this workflow.

## Runtime and validation

`environment.yml` and `container.def` identify the digest-pinned DMRcate 3.6.0
BioContainer. New staging uses `apptainer pull` directly, with no added package
installation or custom image build. Staging freezes the SIF checksum, package
inventory and GENCODE download provenance. Analysis performs no installation
or annotation download and uses module `containers/apptainer/1.3.6`.

The configurations still select the existing frozen runtime from September 29;
it contains the historical methrix installation, which this code no longer
loads or needs. `runtime.lock.json` remains the exact record of that existing
runtime. A future direct-pull SIF needs its own staging/validation record and
an explicit configuration update; it must not replace the existing SIF in place.
That new staging route has not yet been executed.

```bash
w=/path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/dmrcate
root=/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate

# Optional: stage a NEW frozen runtime. Never overwrite an existing runtime.
stage="$root/$(date -u +%Y-%m-%d_%H-%M)-runtime-bioc322-v1"
sbatch --output="$root/stage-%j.out" "$w/stage_runtime.sh" "$w" "$stage"

# Validate the selected runtime and a frozen snapshot of the current code.
# For the existing runtime, obtain its explicit path from config/smoke.json.
runtime="$(python -c 'import json,sys;print(json.load(open(sys.argv[1]))["runtime"])' "$w/config/smoke.json")"
validation="$root/$(date -u +%Y-%m-%d_%H-%M)-dmrcate-lean-fixtures-v1"
sbatch --output="$root/fixtures-%j.out" "$w/validate_runtime.sh" "$w" "$runtime" "$validation"
```

After fixtures succeed, `runtime_validation` in the selected configuration
must identify that completed bundle. The launcher checks the SIF checksum and
that `analysis.R` and `prepare.py` match the tested source hashes before it
creates outputs or submits analysis. Old fixture acceptance cannot silently
approve changed analysis code. See [implementation status](IMPLEMENTATION_STATUS.md)
for actual execution evidence and pending acceptance gates.

## Submit analysis

```bash
# Smoke: chr22:20,000,000–21,000,000, one-based closed, all twelve samples.
smoke="$(date -u +%Y-%m-%d_%H-%M)-dmrcate-smoke-all12-v1"
bash "$w/run_dmrcate.sh" "$w/config/smoke.json" "$smoke"

# Genome-wide cohort, after accepting the current implementation's checks.
cohort="$(date -u +%Y-%m-%d_%H-%M)-dmrcate-cohort-all12-v1"
bash "$w/run_dmrcate.sh" "$w/config/cohort.json" "$cohort"

# Resume an incomplete run with identical frozen code, configuration and inputs.
bash "$w/run_dmrcate.sh" "$w/config/cohort.json" "$cohort" --resume
```

The launcher submits through `sbatch`; contrasts run sequentially. Each
contrast's model, smoothing and table generation run in separate fresh R
processes. The model still fits its complete eligible genome-wide CpG set;
stock DMRcate still smooths that set and applies its original genome-wide
corrections. No chromosome/block splitting or package patch is used.
The launcher checks the complete output contract after all three contrasts
before writing `COMPLETE`.

Smoke uses
8 CPUs, 32 GB, 24 hours on `mediumq`; the full cohort defaults to 8 CPUs,
64 GB and 24 hours on `mediumq`. The completed fresh full-cohort run took
4h38m40s and peaked at about 43.2 GiB, supporting the reduced request.
`DMRCATE_MEMORY` remains available as an explicit cohort memory override.

```bash
cohort="$(date -u +%Y-%m-%d_%H-%M)-dmrcate-cohort-all12-v3"
bash "$w/run_dmrcate.sh" "$w/config/cohort.json" "$cohort"
```

This fresh invocation reads the original reports, imports new disk-backed
counts and recomputes QC and all three models. It does not resume or copy
earlier results. Existing runs retain their own frozen code and paths.

The earlier smoke peaked
at 5.36 GiB. The first cohort attempt failed during coordinate validation after
import; its 7.16 GiB peak does not measure full inference. Requests, submission times,
actual UTC execution start/completion and attempt history are recorded.

Count validation, coverage QC and eligibility read one million loci per block,
matching the native bsseq HDF5 chunks. This is an I/O setting; genome-wide
models and BH correction remain unsplit.

Inputs and references are mounted read-only. R receives an explicit scratch
`TMPDIR`; preflight checks temporary storage and reference access. Compressed
reference tables are parsed only after successful gzip decompression. Code and
input hashes prevent changed inputs or code on resume. Completed bundles cannot
be overwritten. Partial count stores remain as recovery evidence; imports use
a new internal count directory on retry. Existing runs with earlier code need
their original frozen implementation for recovery; the lean launcher does not
resume them or regenerate their reports.

All reference contigs are retained for import/count-conservation checks and
coverage QC. The wgbstools chromosome list includes chrM, whose CpGs also use
the existing namespace. Inference still uses only chr1–22, X and Y.
Biostrings indexes the existing FASTA and reads the remaining contigs in one
native call; the workflow creates no reference index sidecar and changes no
reference resource. The October 9 contig fix addresses
the first cohort attempt's incorrect assumption that every contig outside the
wgbstools chromosome list was lambda or pUC19. Its retry uses a fresh run and
import; the failed run's partial counts remain intact.

The focused process-boundary fixture is submitted with
`validate_checkpoint.sh WORKFLOW RUNTIME NEW_VALIDATION_DIRECTORY`. It
compares model CpG statistics and smoothing exactly with direct stock calls,
checks the standard tables and completion contract, and exercises preserved
checkpoints, completed-result protection and empty smoothing. It also runs
the unchanged Python input fixtures. Its `VALIDATED` marker, runtime checksum
and `r/validated_code.tsv` can be selected through `runtime_validation`.
The broader `validate_runtime.sh` suite additionally retains the native
import, reference-contig, paired-design and empty three-contrast checks and
now includes the process-boundary fixture.

## Recover smoothing from an existing CpG checkpoint

For an incomplete contrast that already contains `cpgannotated.rds` and
`eligible_rows.rds`, `recover_checkpoint.sh` skips the completed model and runs
stock DMRcate 3.6.0 smoothing in a fresh R process. A second fresh process
creates the usual region, effect and gene tables and `INFERENCE_READY.rds`.
The original run directory, start timestamp, configuration, tested CpG universe,
frozen code and CpG statistics are retained. No DMRcate package patch is applied.
The recovery source, validation receipt, timestamps, logs and resource reports
are saved in a new dated directory under the existing run's `recovery/` section.

```bash
validation="$root/$(date -u +%Y-%m-%d_%H-%M)-dmrcate-checkpoint-fixtures-v1"
sbatch --output="$root/checkpoint-fixtures-%j.out" \
  "$w/validate_checkpoint.sh" "$w" "$runtime" "$validation"

# After VALIDATED.json exists, resume one explicitly selected saved contrast.
run="$root/2026-10-20_08-12-dmrcate-cohort-all12-v2"
bash "$w/recover_checkpoint.sh" "$run" cfDNA_minus_BM "$validation/VALIDATED.json"
```

This entry point submits through Slurm with 8 CPUs, 294 GB and 96 hours on
`longq`, retaining the most recent cohort allocation. It reuses the run's
`config.json`, `prepared/samples.tsv`, runtime and original input-fingerprint
checks; it requires a receipt matching the recovery source, frozen analysis and
runtime. It refuses completed contrasts or existing region outputs. An
interrupted publication requires review before retrying. It does not fit a
contrast without saved statistics, and completing one contrast does not mark
the entire three-contrast run complete. The ordinary `--resume` entry point
continues to enforce its original frozen-code contract.

## Essential outputs

| Location | Contents |
|---|---|
| `prepared/` | Sample table, per-contig source totals and reference-locus validation, count-conservation checks, native disk-backed bsseq counts |
| `summary/` | Coverage retention by chromosome, sample, matched pair and contrast |
| `comparisons/<numerator>_minus_<denominator>/` | Designs, eligible sites, CpG statistics, DMR tables/BED, native DMRcate objects when estimable, sample methylation, patient/mean effects, gene overlaps, nearest TSS, inference status/checkpoint |
| `phases/<contrast>/attempt-<job>/` | Separate model/smoothing/reporting timestamps, session information, smoothing result and resource measurements; retained on failure |
| `provenance/`, `code/`, `input_fingerprint.json`, `sessionInfo.txt` | Frozen source/input/runtime provenance and upstream control QC |
| `attempts.tsv`, `resources-<job>.txt`, `validated_outputs.json`, `COMPLETE` | Execution trail, resource measurements and verified completion |

`ind.fdr` is the CpG BH-adjusted value. `min_smoothed_fdr`, `Stouffer`, `HMFDR`
and `Fisher` retain their DMRcate meanings; the complete region list is not
claimed to have guaranteed 5% region-level FDR. `diff`, `maxdiff` and `meandiff`
are transformed model effects. `difference_pp` and `mean_difference_pp` are
methylation differences in percentage points.

Internal coordinates are one-based closed; the explicitly named BED export
is zero-based half-open. Empty tables retain headers and are valid outcomes.
Completion requires every essential output to exist and be nonempty, including
header-only empty-result tables, followed by successful analysis exit.

The small local annotation adapter is retained deliberately:
`DMRcate::extractRanges()` in the pinned release retrieves an ExperimentHub
annotation, which would change the frozen annotation contract. We retain native
DMRcate objects and use the already approved GENCODE release without that fetch.

## Sources

- [DMRcate Bioconductor 3.22](https://bioconductor.org/packages/3.22/bioc/html/DMRcate.html)
- [bsseq Bioconductor 3.22](https://bioconductor.org/packages/3.22/bioc/html/bsseq.html)
- [GENCODE v50](https://www.gencodegenes.org/human/release_50.html)
