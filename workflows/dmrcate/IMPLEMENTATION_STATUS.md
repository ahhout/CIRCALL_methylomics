# DMRcate replacement: undergoing validation

Implementation began from clean tracked Git checkpoint
`f2db536f3c212924dfed3c2fdbd503242bb52977` on 2026-09-29.
The staged list and ignore rules were reviewed; representative BAM, SIF,
generated-result and downloaded-software paths were excluded. Unrelated
untracked work was left untouched.

This is not yet the accepted default. Container construction, fixtures, a
Slurm smoke run and the full cohort must succeed before that designation.
Existing DMRichR code, runtime, results and jobs remain in place.

## Runtime and validation record

- Runtime: `/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-09-29_14-28-runtime-bioc322-v2`.
- Build job `5920033` produced the SIF successfully, then failed to connect to
  GENCODE over HTTPS. Annotation recovery job `5920073` succeeded on retry using
  the same official HTTPS URL. The SIF was verified and reused, not rebuilt.
- `runtime.lock.json` records the SIF digest, 329 Conda artifact records, 237
  resolved R package versions, methrix source checksum and both annotation
  checksums/retrieval records.
- Fixture job `5920038` passed strand/count conservation, coordinate/promoter,
  sample permutation and reversed-contrast checks, then exposed a Bioconductor
  safeguard against serializing the lazy methrix view. The view now explicitly
  retains absolute references to the authoritative HDF5 count store through
  `base::saveRDS`; no count matrix is reconstructed or duplicated.
- Fixture job `5920090` tests that correction, the three synthetic contrasts,
  empty-region reports, independent coverage universes and inference recovery.
  It completed successfully in 6 minutes 43 seconds with 5.22 GiB peak memory.
  Every scientific fixture and report/recovery check passed. A nonfatal font
  cache warning occurred; figures rendered and passed validation. The production
  launcher now supplies a writable cache directory.
- Real-data smoke job `5920106` was submitted after those checks passed. It
  exercises chr22:20,000,000–21,000,000 across all twelve original samples and all
  three contrasts. Its outputs have not yet been accepted.
- No full-cohort run has been submitted yet. The default-workflow transition
  remains gated on smoke review and a successful full run.

`validation.json` records the exact validation bundle, checks, R source digest,
resource measurements and smoke-run identity. The report-recovery input helper
also passed its three targeted Python tests after its final change.

Earlier new-workflow staging jobs `5919976`, `5919981` and dependent fixture job
`5920013` were superseded before analysis. The first pending submission's naming
timestamp was corrected; the subsequent build was restarted with node-local
temporary extraction to avoid slow shared-filesystem metadata operations. Its
downloaded container cache was reused. These changes did not touch DMRichR jobs.

## Smoke allocation adjustment, 1 October 2026

The user approved separating smoke resources from the full cohort. Before
changes, the tracked tree and staged list were clean; verified checkpoint
`3b89a25cb364d42da83185d503ac2ab791eb3ede` protects the tracked state. DMRcate
source was previously untracked; its existing run snapshots retain the earlier
code. Ignore rules and representative BAM/SIF/result exclusions were reviewed.

The launcher now selects 8 CPUs, 32 GB and 24 hours on `mediumq` for a non-null
`smoke_interval`. Full-cohort requests remain 8 CPUs, 196 GB and 96 hours on
`longq`. The smaller request allows headroom over the 5.22 GiB fixture peak,
including annotation loading and scans of all twelve source reports. It is an
initial bounded request, not a measured real-data requirement. There is no
automatic escalation. Scientific parameters, runtime and sample selection are
unchanged. Future submissions append `resource_requests.tsv`.

Pending job `5920106` was updated in place, retaining its original submission,
run directory and frozen code. Slurm verified the new resources at
`2026-10-01T11:30:18.170825+00:00`; execution had not started at that check.
Exact before/after scheduler records and operation start/completion timestamps
are in `resource_adjustment_20261001_retry.json` under the original
`2026-09-29_14-51-dmrcate-smoke-all12-v1` run. An initial attempt made no change:
this Slurm command requires numeric MiB (`MinMemoryNode=32768`), not `32G`.
Its failed-attempt record is retained alongside the successful one.

Validation: `bash -n workflows/dmrcate/run_dmrcate.sh` passed.
The launcher's actual resource-selection block was executed against both
`config/smoke.json` and `config/cohort.json`; both matched the requests above.
No biological tests or full-cohort jobs were added for this runtime-only change.

## Temporary storage and reference visibility recovery, 5 October 2026

Smoke `5920106` failed after 41m50s with 1.72 GiB peak RSS, not an OOM. Its
container gave R a 64 MiB temporary filesystem: the host `TMPDIR` did not pass
through `--cleanenv`. gzip exhausted this filesystem; `fread(cmd=...)` then read
a truncated dictionary and failed coordinate validation. The reference manifest
and FASTA also lacked a container bind and would have failed subsequently.

From clean tracked checkpoint `f9c4e48a86e944d509e6a28d4df7e0d5458a0100`, the
launcher now passes `TMPDIR` explicitly, mounts existing reference files
read-only and checks temporary-directory writability and reference accessibility
inside the container before scanning sample reports. gzip exit status is checked
before either compressed reference table is parsed. Fixture validation also
receives the explicit temporary path. No biological settings, input references,
container versions or resource limits changed; the failed run is preserved.

Targeted checks passed in the pinned container: valid gzip, truncated gzip,
missing gzip, missing reference, scratch temporary storage and actual reference
visibility. Evidence (including frozen tested code, UTC check timestamps and
filesystem details):
`/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-05_14-00-dmrcate-runtime-check-v1`.
The three `test_prepare.py` tests, Bash syntax checks for both launchers and
`git diff --check` passed. Existing scientific fixtures were not rerun merely
for runtime housekeeping; the new real-data smoke remains the acceptance gate.

Fresh smoke **5991522** was submitted with:

```bash
bash workflows/dmrcate/run_dmrcate.sh \
  workflows/dmrcate/config/smoke.json \
  2026-10-05_14-02-dmrcate-smoke-all12-v2
```

Output:
`/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-05_14-02-dmrcate-smoke-all12-v2`.
It was RUNNING at the post-submission check with 8 CPUs, 32 GB and 24 hours on
`mediumq`. No `--resume` was used because the frozen analysis code changed;
original run fingerprints and partial count stores remain intact. Per-attempt
execution timestamps are in the new run. Full-cohort execution remains gated.

## Smoke acceptance review, 5 October 2026

Job **5991522** subsequently completed with Slurm state `COMPLETED`, exit `0:0`,
elapsed **01:00:07** and peak batch RSS **5,618,300 KiB** (5.36 GiB).
The recorded analysis execution was `2026-10-05T14:03:03Z` through
`2026-10-05T15:03:03Z`. This was the existing 12-sample chr22:20–21 Mb smoke,
using `config/smoke.json`, the pinned Bioconductor runtime and no `--resume`.

Acceptance review independently summed `prepared/smoke_expected_counts.tsv`
and matched methylated and total observations for all twelve samples against
`prepared/count_conservation.tsv`. The analysis also checks individual CpG
counts before writing this table. Every artifact listed in
`validated_outputs.json` exists and is nonempty. All 24 local HTML report links
resolve; 11 PDF reports were produced and passed the analysis PDF checks.
All three patient-matched contrasts have CpG statistics, region tables,
annotations, patient effects and inference-ready objects. The smoke returned
one BM−PB region, zero cfDNA−BM regions and one cfDNA−PB region; the empty
region result is valid, not a missing contrast. These are bounded smoke outputs,
not full-cohort genome-wide conclusions.

The temporary-storage smoke gate has passed. No full-cohort job was submitted
during this review. Full execution retains the approved matched design and
separate descriptive deconvolution context; PB remains diagnostic material.

## Minimal implementation, 6 October 2026

The user requested reducing custom overhead to the strictly necessary layer and
selected minimal analysis/essential tables, dropping bespoke reports and
exploratory outputs. Existing checkpoint `066defd85578645087e9e8dd355a70c9cbabb431`
protected the clean tracked starting state. Before significant documentation
changes, the reviewed code was committed locally as `611ae946c99e936b50b3a5bc22323326a6b692bc`,
followed by verified checkpoint `70bb90cf862546995bd1d78703f0ebbb1a65fd87`
for the frozen-fixture path correction and final staging simplification. Staged lists and ignore rules were reviewed;
representative BAM, SIF and result paths remained ignored. Unrelated untracked
files were excluded. No historical runs or data were moved or deleted.

The lean contract removes bespoke HTML/PDF reports, report-only recovery,
PCA/clustering/correlation/distribution output, methrix wrapping/installation,
and detailed promoter/exon/intron/island/shore/shelf categories. Compact gene
overlap/nearest-TSS and sample/patient/mean effect tables remain. All original
input IDs, inference chromosomes, coverage eligibility, pairing, DMRcate
thresholds and selected reference/runtime paths are unchanged. Saved
deconvolution results are separate context, no longer copied into each run.
New runtime staging directly pulls the same digest-pinned BioContainer and only
retrieves the retained GENCODE annotation. The existing SIF and its historical
lock record are reused unchanged; a new stock-image pull has not been executed.

Production code is reduced from 815 to 612 physical lines (about 25%), including
comments/blanks; R analysis falls from 419 to 268 lines (about 36%). Those counts
cover analysis, preparation, launcher, staging and image descriptor, excluding
tests, configurations, documentation and installed package source. This measures
maintenance surface, not the fraction of statistical computation implemented
locally. Remaining custom code preserves the approved scientific contract and
HPC/provenance checks. The pinned DMRcate `extractRanges()` was inspected; its
ExperimentHub annotation fetch is why a small adapter for the existing GENCODE
annotation remains.

The launcher now requires current analysis/preparation source hashes to match
a successful fixture bundle before submission. Validation tests use a frozen
source/config snapshot. Initial fixture job `6000056` failed before synthetic
cohort analysis because its config still pointed at the live sample table,
which was not mounted after switching to the snapshot. The fixture now points
at its identical snapshotted sample table; failed evidence is retained. Retry
job `6000087` uses
`/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-06_09-40-dmrcate-lean-fixtures-v2`.

Python input fixtures (two retained tests), Bash syntax checks and
`git diff --check` passed. A launcher check rejected old validation before
creating a run or submitting a job. Full-cohort execution remains pending;
the earlier accepted smoke is evidence for the prior output contract.

Lean acceptance evidence:

- Fixture job `6000087` completed (`0:0`) in 5m53s with 4,332,624 KiB
  peak RSS (4.13 GiB). Execution: `2026-10-06T09:40:53Z` to
  `2026-10-06T09:46:44Z`. All ten retained scientific, empty-table and
  inference-recovery checks passed. Both configurations now select this bundle.
  Their launcher validation gates pass; the old source hashes are rejected.
- A bounded regression reused the accepted smoke's authoritative HDF5 count
  store read-only, without copying counts or rescanning source reports. Bundle:
  `/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-06_09-44-dmrcate-lean-saved-smoke-check-v2`.
  Analysis execution: `2026-10-06T09:45:23Z` through
  `2026-10-06T09:48:22Z` as recorded in `validated_outputs.json`. Original configuration: `config/smoke.json`;
  native DMRcate 3.6.0; no `--resume`; source upstream methylseq/reference
  provenance is retained in the original accepted smoke bundle.
- Independent `code/numerical_review.py` passed **35 comparisons**. DMR
  coordinates/counts, significance flags, eligible indices, designs, coverage
  retention, sample/patient/mean effects and compact gene/TSS annotations match
  exactly. Native statistics/coefficient differences have maximum absolute
  magnitude `1.0658141036401503e-14`; review uses absolute/relative tolerances
  of `1e-12`. Near-zero coefficients can have large relative differences, so
  absolute differences and exact discovery decisions are also checked.
- Strict-comparison job `6000140` ended `FAILED` because its zero-tolerance
  assertion rejected floating-point roundoff after inference/exports completed.
  That failed scheduler state remains recorded; the numerical review reuses
  those completed fits and does not alter their values. An earlier temporary
  comparator had a syntax error; its superseded job `6000123` was cancelled
  and the original script/evidence remains intact.
- Every essential output exists and is nonempty, including header-only empty
  tables. Region counts remain 1 BM−PB, 0 cfDNA−BM and 1 cfDNA−PB. These
  are bounded-smoke results. No bespoke reports or exploratory plots were
  produced. The direct stock-image staging path is syntax-checked but unrun;
  the selected SIF and runtime lock remain unchanged. Full-cohort execution
  and final default-workflow acceptance remain pending.

`validation.json` preserves the earlier validation record and adds this lean
acceptance record, exact evidence locations, checks, resource measurements and
limitations.

## Full-cohort submission, 6 October 2026

The user authorized proceeding with the genome-wide run. Job **6000293**
was submitted at `2026-10-06T10:07:29.090838+00:00` using DMRcate 3.6.0,
`config/cohort.json`, `config/samples.tsv` (all twelve original samples),
the existing pinned runtime and the completed lean fixture bundle.
No `--resume` was used. Scientific settings and resource requests are unchanged:
8 CPUs, 196 GB, 96 hours on `longq`.

Output: `/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-06_10-07-dmrcate-cohort-all12-v1`.
The source revision is `b819ef100ac2ba823da10f64d92e02a27c1511c4`; code and configuration
were snapshotted before submission. Submission/naming UTC times are recorded
separately from actual execution start, which the worker will record when it runs.

At `2026-10-06T10:08:28.352067+00:00`, Slurm reported **PENDING (Priority)**
and verified the requested allocation. Execution had not started; completion
and cohort output acceptance remain pending. Exact submission metadata is in
`runs/2026-10-06_10-07-dmrcate-cohort-all12-v1_submission.json` and the run directory.


## Full-cohort contig failure and correction, 9 October 2026

Job `6000293` failed (`1:0`) after native count import. Actual worker execution
started at `2026-10-07T20:43:33Z` and ended at `2026-10-07T21:43:45Z`.
Slurm elapsed time was 1h00m20s, peak RSS 7,510,572 KiB (7.16 GiB).
No completed DMR comparisons were produced. The failed run and partial native
count store remain intact at the original path.

The custom coordinate check incorrectly treated every contig outside the
wgbstools chromosome-size list as lambda or pUC19. Whole-genome reports
also contain human alternate, random and unplaced contigs; the
primary-chromosome smoke interval did not exercise this assumption.

Initial checkpoint `8033851e01d8191ca75cc6710078a2701c06c022` fixed this check.
CpGs on dictionary-listed chromosomes, including chrM, still use the existing
namespace. Other contigs use native
Biostrings indexing and selective reads of the same pinned combined FASTA,
one contig at a time, with explicit rejection of unknown, out-of-bounds and
non-CpG coordinates. No reference sidecar is generated. The existing runtime
already includes Biostrings; its environment record now declares that use.
No runtime, reference, sample, coverage filter, model or DMRcate setting changes.
Inference remains restricted to chr1–22, X and Y; all contigs remain available
for count conservation and coverage QC. A small per-contig validation table
records how each imported sequence was checked.

The retry will use a fresh run and standard native import. This avoids adding
partial-store recovery machinery to the minimal workflow. Validation uses a
frozen fixture snapshot and a separate read-only audit of the failed count store
against the actual reference and original source count totals. No statistical
fits are repeated for the audit. Exact evidence and scheduler records are in
`validation.json`; final cohort acceptance remains pending.


Contig-fix fixture job `6035774` completed (`0:0`) in 16m48s, peak RSS
5,090,620 KiB (4.85 GiB), using the existing 1-CPU/8-GB/2-hour `shortq`
validation launcher. All twelve R checks and both Python fixtures passed.
The new checks accept primary, alternate, random/unplaced, mitochondrial and
control contigs present in a fixture FASTA; they reject unknown sequences,
non-CpG sites, out-of-bounds positions, missing primary namespace entries and
duplicate collapsed coordinates. The ten retained statistical/output checks
also pass. After this initial validation, both configurations selected:
`/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-09_10-10-dmrcate-contig-fixtures-v1`.
Fixture execution: `2026-10-09T10:10:08Z` to `2026-10-09T10:26:53Z`.
R checks: `2026-10-09T10:18:03Z` to `2026-10-09T10:26:51Z`.
Bash syntax and `git diff --check` passed. The inference/output source after
coordinate validation is byte-for-byte unchanged; the new code snapshots match
the current source. The separate actual-cohort audit is job `6035777`, with
2 CPUs, 32 GB and 2 hours on `shortq`; it runs no inference and has read-only
access to the failed count store and reference.


The first actual-cohort audit (`6035777`) was cancelled after 30m10s. It
exposed excessive reference I/O: after scanning the FASTA index, individual
indexed record reads scanned preceding sequence again on this runtime.
Read-only counters showed more than 5 GB read before locus validation finished,
and a new read progressing through primary sequence to reach the first
non-primary record. No data were altered and the original audit is preserved.
The correction requests all non-primary records in one native Biostrings call,
then checks each against only the non-primary loci. Primary sequences are not
retained in that DNAStringSet. This preserves the same checks and avoids
reopening the file for each record. A primary-only fixture also exercises the
empty-subset branch. Corrected fixture job: `6035918`; corrected actual-cohort
audit: `6035919`. The launcher continues to reject changed source until its
matching fixture bundle succeeds; cohort submission remains gated on both checks.


Corrected batch-reader fixture job `6035918` completed (`0:0`) in 14m51s,
peak RSS 5,009,640 KiB (4.78 GiB). All twelve R checks and both Python
fixtures pass, including the added primary-only case. Both configurations
select its matching frozen bundle:
`/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-09_10-41-dmrcate-contig-fixtures-v2`.
Actual fixture execution: `2026-10-09T10:42:07Z` to `2026-10-09T10:56:56Z`.
R checks: `2026-10-09T10:50:51Z` to `2026-10-09T10:56:54Z`.
Final reader source checkpoint: `3137e0905269ade25b22d0e715ec50c378745393`.


The corrected reference check passed **all 30,983,025 loci on 457 contigs**
in audit `6035919`. ChrM is already in the wgbstools dictionary; alternate and
unplaced human sequence plus controls use the combined FASTA. Inference remains
restricted to chr1–22, X and Y. The following count scan exposed a second
operational bottleneck. `h5dump -pH` showed signed 32-bit M/Cov datasets with
native chunk dimensions `(1, 1000000)` (sample × locus), while the custom reader
requested 100,000 loci at a time. This repeatedly decompresses the same chunks.
The audit was cancelled after 1h13m24s; its completed coordinate evidence and
all other files remain intact. No completed count-conservation result is claimed.

The production helper now reads one million loci per block, matching native
bsseq storage. It is shared with the count-only validation rather than copied.
Every locus remains included; integer sums, coverage thresholds, paired designs,
models and DMRcate settings are unchanged. Genome-wide fits and BH correction
are not split into blocks. Fixture job `6036293` validates the final source.
Count audit `6036422` reuses the passed coordinate check, explicitly verifies
that the coordinate-validation function is identical and that coordinate order
and sample identities match, then tests all integer M/Cov totals using the
production helper. The large reference scan is not repeated.


## Final validation accepted, 20 October 2026

The interrupted session had completed both checks on 9 October but had not
selected the final fixture bundle or submitted the cohort retry. On 20 October,
the current analysis/preparation hashes, sample table and scientific settings
were checked against the final frozen fixture snapshot and match. The clean
starting checkpoint for this bookkeeping is
`5e3ffc343532866511a82c09bb9bd77f1f4454a3`.

- Final fixture job `6036293`: **COMPLETED, 0:0**, 9m25s, peak RSS
  4,820,124 KiB. All twelve R checks and both Python fixtures passed.
  Actual fixture execution: `2026-10-09T11:55:40Z` to `2026-10-09T12:05:04Z`.
- Aligned-block count audit `6036422`: **COMPLETED, 0:0**, 5m40s, peak RSS
  2,565,060 KiB. Its R execution ran from `2026-10-09T11:59:30Z` to
  `2026-10-09T12:04:48Z`; its completion marker was written at
  `2026-10-09T12:04:49Z`. All twelve samples' integer methylated and total-read
  counts are conserved across 30,983,025 loci. The passed 457-contig coordinate
  check was reused with identical validation code and coordinate/sample order.
- Both configurations now select the final frozen fixture bundle:
  `/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-09_11-52-dmrcate-contig-fixtures-v3`.
  Runtime, reference, original sample IDs, biological thresholds and paired
  models are unchanged. The failed cohort and superseded audit evidence remain
  intact. Full-cohort result acceptance is still pending.


## Cohort retry submitted, 20 October 2026

Job **6160387** was submitted at `2026-10-20T08:12:45Z` using DMRcate 3.6.0,
`config/cohort.json`, `config/samples.tsv` (all twelve original samples), the
unchanged pinned runtime and final accepted fixture bundle. No `--resume` was
used. The failed run remains intact. Resources remain 8 CPUs, 196 GB and
96 hours on `longq`; scientific settings and sample identities are unchanged.

Output: `/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-20_08-12-dmrcate-cohort-all12-v2`.
Source revision: `ea200c4fa8a92653a481c6e255da6e220b481fa9`.
Submission/naming timestamps, frozen code/configuration, resource requests and
job ID are recorded in the run directory and `runs/2026-10-20_08-12-dmrcate-cohort-all12-v2_submission.json`.
Actual execution start and completion will be recorded by the worker; cohort
results and final replacement acceptance remain pending.

Slurm confirmed **RUNNING on n04**, with 8 CPUs, 196 GB and a 4-day limit.
Scheduler allocation began at `2026-10-20T08:12:45Z`. The worker recorded
actual execution start at `2026-10-20T08:12:52Z`.

## Saved-CpG smoothing recovery, 26 October 2026

The user authorized fresh-process recovery after review of the DMRcate 3.6.0
source and the two smoothing OOM failures (196G at chr7, then 294G at chr12).
The clean tracked starting checkpoint was
`83ced975e22c37d5af4c572cd95b1644a06966c7`; unrelated untracked files were left
alone. The original analysis, preparation code, configuration, environment,
reference, sample IDs, CpG model and smoothing settings are unchanged.

`recover_checkpoint.sh` submits a selected saved contrast through Slurm.
`recover_checkpoint.R` first loads `cpgannotated.rds` and invokes the frozen
run's `call_regions()` helper in a fresh R process. A separate process loads
the disk-backed counts and creates the standard tables using the original
annotation/publication calculations. No DMRcate package patch, chromosome
splitting, block smoothing or added inter-chromosome collection is included.
Existing checkpoints and completed contrasts are preserved. Recovery code,
validation receipt, session information, stage timestamps and resource reports
are retained under a new dated internal recovery attempt; the original run
timestamp and input fingerprint remain intact.

Validation job **6200859** completed with exit **0:0**, elapsed **7m02s** and
peak RSS **5,327,760 KiB**. Worker execution ran from
`2026-10-26T11:17:11Z` to `2026-10-26T11:24:11Z`. The receipt is
`/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-26_11-16-dmrcate-checkpoint-fixtures-v1/VALIDATED.json`.
Checks cover exact equality with stock DMRcate on a nonempty synthetic
contrast, fresh-process smoothing/publication, unchanged checkpoint bytes,
completed-contrast overwrite rejection and empty smoothing.

Commands executed:

```bash
bash -n workflows/dmrcate/recover_checkpoint.sh workflows/dmrcate/validate_checkpoint.sh
git diff --check
sbatch --parsable --output=/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/checkpoint-fixtures-%j.out \
  /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/dmrcate/validate_checkpoint.sh \
  /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/dmrcate \
  /path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-09-29_14-28-runtime-bioc322-v2 \
  /path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-26_11-16-dmrcate-checkpoint-fixtures-v1
```

The target is the existing `cfDNA_minus_BM` checkpoint in
`2026-10-20_08-12-dmrcate-cohort-all12-v2`, using its frozen `config.json`,
prepared sample table and DMRcate 3.6.0 runtime. This is checkpoint recovery;
the original preparation helper is invoked with `--resume`. Resources retain
the most recent allocation: 8 CPUs, 294G, 96 hours, `longq`. Production success
is not established by the synthetic fixtures. `cfDNA_minus_PB` has no saved
CpG checkpoint yet; completing this recovery does not complete the entire run.

The first submission, job `6200885`, exited before analysis in its source
preflight: Slurm's spool path was incorrectly used to locate the recovery R
worker. The launcher now resolves the worker from its explicit snapshot
directory and allocates a new attempt version if the minute already exists.
The failed attempt is retained; no CpG checkpoint or contrast output changed.
The R worker and its successful fixture receipt are unchanged.

Corrected recovery job **6200887** was submitted at `2026-10-26T11:26:19Z`;
the worker recorded `started_utc=2026-10-26T11:26:26Z`. Slurm confirmed
**RUNNING on n03**, 8 CPUs and 294G. It passed the frozen-code and runtime
preflight. Source revision: `7c03fe8edf2d9229628a33e1dc4b83be04e31e4b`.
The exact command, configuration, sample table, receipt, attempt paths and
initial status are recorded in
`runs/2026-10-26_11-26-dmrcate-cpg-recovery-6200887.json`. Production completion
and the memory outcome remain pending.

## Fresh full-cohort run with separate R processes, 26 October 2026

The user requested a fresh run from start to finish after successful recovery
of cfDNA−BM. Recovery job `6200887` completed at `2026-10-26T12:05:05Z`, exit
`0:0`, elapsed 38m45s. Its isolated smoothing process peaked at 17,646,192 KiB
(about 16.8 GiB), completed every chromosome and produced the expected tables
and inference checkpoint. No DMRcate package patch was needed.

The canonical launcher now runs each contrast's model, smoothing and reporting
in separate R processes. It starts from the original CpG reports, imports and
validates a new native disk-backed count store, computes QC, and fits all three
genome-wide paired models. Model-only selection changes process scheduling;
the scientific calculations and tested universes are unchanged. Smoothing and
reporting use the same recovery worker that completed cfDNA−BM. A final fresh
process checks every contrast's required outputs before publishing `COMPLETE`.
The default direct `main()` implementation remains available to fixtures.

Reviewed source checkpoint `1087ce561274cb291e2f7e706eade86007cec9ba` precedes
the documentation/configuration update. Unrelated reference-audit work was
left alone. The only configuration change is selection of the passed
process-boundary validation record in `config/cohort.json` and `config/smoke.json`;
runtime, references, samples, thresholds and all three contrasts are unchanged.
The full validation entry point now also includes the process-boundary tests.

Focused fixture job **6201055** completed with exit **0:0** in **9m41s**, peak
RSS **5,693,812 KiB**, from `2026-10-26T12:37:20Z` to
`2026-10-26T12:47:01Z`. Its frozen sources match the production launcher,
analysis, model worker, recovery worker and unchanged preparation helper.
Both Python input fixtures passed. The R fixture confirms exact model CpG
statistics and smoothing against direct stock calls, complete final outputs,
checkpoint-byte preservation, completed-result protection and empty smoothing.
This extends the earlier import/reference/paired-design evidence; those
unchanged expensive checks were not rerun in a separate audit job.

Validation executed:

```bash
bash -n workflows/dmrcate/run_dmrcate.sh workflows/dmrcate/validate_checkpoint.sh workflows/dmrcate/validate_runtime.sh
git diff --check
sbatch --parsable --output=/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/fresh-process-fixtures-%j.out \
  /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/dmrcate/validate_checkpoint.sh \
  /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/dmrcate \
  /path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-09-29_14-28-runtime-bioc322-v2 \
  /path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-26_12-37-dmrcate-fresh-process-fixtures-v1
```

The requested fresh cohort uses DMRcate 3.6.0, the canonical cohort config and
sample sheet, 8 CPUs, 294G and 96 hours on `longq`. No `--resume` is used.
Older result bundles and original run timestamps are retained. Production
completion remains pending and the synthetic fixture does not establish
biological validity or final default-workflow acceptance.

Fresh cohort job **6201112** was submitted at `2026-10-26T12:48:55Z` and
recorded worker start at `2026-10-26T12:49:10Z`. Slurm confirmed **RUNNING on
n04**, 8 CPUs, 294G and a 96-hour limit. Frozen-source checks passed.
Source revision: `654ec44521a9c08f35b02bec1abb265af2b217f1`.
Output: `/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-26_12-48-dmrcate-cohort-all12-v3`.
The exact invocation, source hashes, settings and initial status are recorded
in `runs/2026-10-26_12-48-dmrcate-cohort-all12-v3_submission.json`:

```bash
DMRCATE_MEMORY=294G bash /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/dmrcate/run_dmrcate.sh \
  /path/to/hpc-home/analysis/CIRCALL_methylomics/workflows/dmrcate/config/cohort.json \
  2026-10-26_12-48-dmrcate-cohort-all12-v3
```

## Resource defaults reduced, 27 October 2026

Fresh cohort job `6201112` completed successfully in 4h38m40s with peak Slurm
RSS 45,247,160 KiB (about 43.2 GiB). All three contrasts and the final output
contract completed. The user authorized `mediumq` and 64G as the new cohort
defaults. Slurm reports a 24-hour maximum for `mediumq`, so the cohort time
limit is now 24 hours. CPUs remain eight; the explicit `DMRCATE_MEMORY`
override remains available. Smoke stays on `mediumq` with 32G and 24 hours.
The standalone historical checkpoint-recovery launcher is unchanged.

The clean tracked pre-change checkpoint was
`3b8f7f50e3dcf50fe5a7dae11725f5c5256d5f17`. Only the cohort resource line changed
in the launcher. `bash -n` passed; shell evaluation verified the cohort default
(`mediumq 64G 24:00:00`), an 80G override and the unchanged smoke profile.
All four validated scientific/preparation source hashes match the passed
26 October fixture. Those expensive fixtures were carried forward rather
than rerun. Both configs select the derived resource-only acceptance record:
`/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-27_08-13-dmrcate-resource-validation-v1`.
Its source receipt retains the original scientific fixture evidence and
records the exact launcher diff and UTC validation times. The same record is
tracked in `runs/2026-10-27_08-13-dmrcate-resource-validation-v1.json`.
Scientific parameters, environment and completed run bundles are unchanged.
No new pipeline job was submitted for this resource-only change.
