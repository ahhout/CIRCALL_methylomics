# Where to find the current results

**All 12 samples, saved deconvolution estimates:**
[EpiDISH R](workflows/epidish/exports/2026-09-29_10-25-deconvolution-all12-display-v2/epidish_r/),
[native UXM](workflows/epidish/exports/2026-09-29_10-25-deconvolution-all12-display-v2/uxm/),
and [method comparisons by original cohort](workflows/epidish/exports/2026-09-29_10-25-deconvolution-all12-display-v2/comparisons/).
[Display settings and provenance](workflows/epidish/RESULTS.md).
This is a combined display of the existing 8+4 samples, not a new joint fit.

Navigation and recorded completion status checked on **2 November 2026 (UTC)**.
Start here when browsing results; this documentation check did not rerun analyses.

**The newest cohort and the established baseline are different datasets.**
The August 31 run (E31) contains the eight samples used by most downstream
projects. The September 16 run contains four new HEH samples and has its own
completed methylseq → PAT → UXM → EpiDISH chain. It does not replace or merge
the eight-sample baseline. Dates in directory names identify runs; they do
not necessarily give the date of completion or the latest figure refresh.

## Methylseq and shared workflows

For direct EpiDISH PDFs, fraction tables and reports, use the
[EpiDISH output guide](workflows/epidish/RESULTS.md).

Choose the column for the samples you want. Links open the explicit published
bundle, including its provenance. PAT stores molecule-level methylation
patterns; UXM and EpiDISH estimate cell-type contributions.

| Analysis | Eight-sample baseline — August 31 | Newest cohort — four new HEH samples, September 16 |
| --- | --- | --- |
| Methylseq alignment, methylation calls and QC | [E31 methylseq](outputs/2026-08-31_10-02-emseq-integrated-8samples-trim-8-3-3-10-v1/) | [HEH methylseq](outputs/2026-09-16_08-08-emseq-integrated-heh-4new-trim-8-3-3-10-v1/) |
| PAT inputs for downstream analysis | [E31 PAT](workflows/uxm/results/2026-09-09_21-50-emseq-integrated-8samples-trim-8-3-3-10-v1-pat-v1/) | [HEH PAT](workflows/uxm/results/2026-09-24_20-56-emseq-integrated-heh-4new-trim-8-3-3-10-v1-pat-v1/) |
| UXM — full atlas | [E31 full atlas](workflows/uxm/results/2026-09-10_11-19-emseq-integrated-8samples-trim-8-3-3-10-v1-uxm-full-atlas/) | [HEH full atlas](workflows/uxm/results/2026-09-24_20-57-emseq-integrated-heh-4new-trim-8-3-3-10-v1-uxm-full-atlas/) |
| UXM — haematopoietic atlas | [E31 haematopoietic](workflows/uxm/results/2026-09-10_12-33-emseq-integrated-8samples-trim-8-3-3-10-v1-uxm-hematopoietic/) | [HEH haematopoietic](workflows/uxm/results/2026-09-24_22-38-emseq-integrated-heh-4new-trim-8-3-3-10-v1-uxm-hematopoietic/) |
| EpiDISH comparison — preferred method for blastic samples | [E31 comparison](workflows/epidish/results/2026-09-10_15-01-emseq-integrated-8samples-trim-8-3-3-10-v1-epidish-comparison/) | [HEH comparison](workflows/epidish/results/2026-09-25_00-04-emseq-integrated-heh-4new-trim-8-3-3-10-v1-epidish-comparison/) |

The [workflow release record](docs/audits/downstream_structure_20260914/releases.md)
documents E31 acceptance. The [September 29 status record](docs/audits/reconciliation_status_20260929.md)
records completion of the four-new-sample chain. Completion alone does not
establish biological validation. See the [EpiDISH README](workflows/epidish/README.md)
for run-specific settings, historical comparisons and interpretation limits.

## Downstream projects: results to open first

These projects answer different questions, so each has its own result entry.
Their newer analysis or plotting dates do not imply newer sequencing inputs.

| Project | Result to open | Dataset, date and status |
| --- | --- | --- |
| PCDH EpiClone | [Figure guide](/path/to/hpc-home/brainstorming/pcdh_epiclone/FIGURES.md) · [current figure bundle](/path/to/hpc-home/brainstorming/pcdh_epiclone/outputs/2026-09-23_08-38-pcdh-figures-v1/) | E31-derived counts; accepted figure refresh September 23. Includes feasibility, coverage, marrow and plasma figures. |
| W footprint — compartment analysis | [Compartment histograms](/path/to/hpc-home/brainstorming/w_footprint/outputs/2026-10-09_10-10-w-footprint-gabbutt-12samples-4patients-v1/compartment_histograms.pdf) · [bundle](/path/to/hpc-home/brainstorming/w_footprint/outputs/2026-10-09_10-10-w-footprint-gabbutt-12samples-4patients-v1/) | All 12 diagnostic samples from four patients; October 9 run has a completion marker and successful provenance. Historical eight-sample runs remain in the project output index. |
| W footprint — PCDH comparison | [PCDH/Gabbutt histograms](/path/to/hpc-home/brainstorming/w_footprint/outputs/2026-09-23_09-32-w-footprint-pcdh-vs-gabbutt-8samples-3patients-v1/pcdh_vs_gabbutt_histograms.pdf) · [bundle](/path/to/hpc-home/brainstorming/w_footprint/outputs/2026-09-23_09-32-w-footprint-pcdh-vs-gabbutt-8samples-3patients-v1/) | E31; September 23. A separate comparison, alongside the compartment analysis. |
| DMRichR — historical paired marrow versus blood | [Completed E31 result](/path/to/hpc-scratch/tmp/em_seq/downstream/dmrichr/2026-09-21_12-12-diagnostic_bm_vs_pb_3patients_emseq_v1/) · [archive and interpretation](archive/downstream/dmrichr/README.md) | Six paired samples from E31; local completion September 21. Workflow archived October 26; results retained. GOfuncR, GREAT and enrichR remain deferred. |
| BCP-ALL atlas — matched-marker application | [Matched 40/41-cell comparison](workflows/uxm/results/2026-09-11_10-07-emseq-integrated-8samples-trim-8-3-3-10-v1-uxm-bcpall-matched907-v0.1/results/matched/) · [project guide](/path/to/hpc-home/brainstorming/bcpall_atlas/README.md) | E31; September 11 experiment. Published bundle has `MATCHED_PAIR_COMPLETED`; disease fractions remain exploratory. |
| DMRcate — paired compartment comparisons | [Completed cohort bundle](/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-26_12-48-dmrcate-cohort-all12-v3/) · [workflow guide](workflows/dmrcate/README.md) | All 12 diagnostic samples from four patients; all three contrasts completed October 26. Final default-workflow acceptance remains pending. |
| DMRcate — healthy-WBC candidate check | [Candidate figure](/path/to/hpc-scratch/tmp/em_seq/downstream/dmrcate/2026-10-27_08-50-healthy-wbc-candidate-comparison-v1/candidate_methylation.pdf) · [run record](workflows/dmrcate/runs/2026-10-27_08-50-healthy-wbc-candidate-comparison-v1.json) | Descriptive check of two retained candidate regions against three external healthy white-blood-cell profiles; October 27. These are not healthy-plasma controls. |

PCDH selection is recorded in the [figure refresh report](/path/to/hpc-home/brainstorming/pcdh_epiclone/runs/figure_refresh_20260923.md).
W-footprint selections are described in its [project README](/path/to/hpc-home/brainstorming/w_footprint/README.md).
The BCP-ALL entry is the later matched-marker application; reference construction
and earlier experiments remain documented in the project guide.

## Experimental comparisons and diagnostics

| Analysis | Evidence to open | Current status |
| --- | --- | --- |
| Native UXM and CELFIE | [Comparison record](workflows/deconflow/ACTIVE_COMPARISON.md) · [all-twelve comparison bundle](/path/to/hpc-scratch/tmp/em_seq/downstream/deconflow/2026-10-26_14-34-uxm-celfie-candidate-comparison-all12-v2/) | Computations and export completed October 26. CELFIE uses a reconstructed megakaryocyte component with unresolved reference incompatibility; fractions remain exploratory. |
| Healthy-WBC CELFIE scoring and convergence | [Repair diagnostic](docs/audits/2026-10-27_09-13-celfie-repair-diagnostics-v1/README.md) | All ten starts converged in the approved longer diagnostic. Reference compatibility and biological interpretation remain unresolved; this did not change the active launcher's defaults. |
| Corrected BWA-MEM3 pilot | [Current operation and handoff](/path/to/hpc-home/brainstorming/bwamem3_pilot/CORRECTED_PILOT.md) | All three sample stages have completion pointers. Final comparison and reconciliation of the shared run summary remain pending as checked November 2. |

The [current audit register](docs/audits/README.md) tracks these open gates and
separates them from historical failures and completed cleanup.

## Keeping this index useful

Update this page after a completed result or accepted figure refresh is
published. Record its cohort, purpose, status and explicit dated path. Add a
separate entry when the samples or scientific question differ. Do not select
results solely by filesystem modification time or a newer directory name.

This is a navigation index checked on the date above, not an automatically
updated status page. Analysis commands should continue to use explicit run
paths and saved provenance. Older run directories and hidden `.running-*`
workspaces are historical or recovery material; they are not the starting
points listed here. [Operational instructions](README.md) and the
[downstream workflow/project index](downstream/README.md) describe how to run
the analyses.
