# Longitudinal EM-seq study of cfDNA origin and clonal dynamics in ALL

## Document status

This document records the scientific purpose and current design of the
enzymatic methylation sequencing (EM-seq) work within the CIRCALL PhD project.
It complements [README.md](README.md), which documents how the current
nf-core/methylseq workflow is run on Flamingo.

The study context follows the investigator-supplied CIRCALL proposal. The aim
order places cell-free DNA origins before integrated clonal dynamics, as
requested by the investigator. The conservative cohort estimates and the hold
on cerebrospinal fluid work were confirmed on 21 September 2026 and are
recorded below.

The design below is a working specification. Items under **Decisions still
required** must be resolved before confirmatory analysis. They must not be
silently inferred from the pilot data.

## CIRCALL context and working hypothesis

CIRCALL is a prospective study of cell-free DNA (cfDNA) in children and
adolescents with newly diagnosed B- or T-cell acute lymphoblastic leukaemia
(ALL). Within PEDIACRIEX, the broader programme developing liquid-biopsy
approaches for paediatric cancers, it asks whether cfDNA measurements can
improve understanding of early treatment
response, residual disease and clonal dynamics. Measurable residual disease
(MRD) refers to leukaemia remaining during treatment that can be detected with
sensitive assays. The proposal plans to compare cfDNA measurements with
conventional bone-marrow MRD and integrate genetic, methylation and clinical
information.

The working hypothesis is that interpreting cfDNA requires both the amount of
leukaemia-derived DNA and the changing contribution of other cells. Combining
these measurements with genetic markers may help distinguish overall disease
clearance from changes in particular leukaemic populations. Associations with
treatment persistence, resistance or relapse remain research questions.

This repository provides EM-seq processing, quality control (QC), methylation
reference development and downstream analysis, chiefly supporting Aims 2 and
3. Molecular MRD assays and clinical outcome analyses belong to the broader
study; their measurements will provide independent comparisons for the
methylation work.

## Project rationale

Preliminary measurements showed that total circulating cfDNA remained abundant
at Day 29 despite low MRD. The observation concerns cfDNA abundance, not a
demonstrated ALL-specific cfDNA signal.

Because cfDNA has a short circulating half-life, cfDNA present at Day 29 cannot
represent molecules physically retained from diagnosis. It must reflect ongoing
DNA release. Short half-life does not, however, distinguish among continued
release from residual ALL cells, normal haematopoietic turnover, inflammation,
treatment-related tissue injury, or other systemic sources.

The project will combine longitudinal EM-seq with molecular MRD measurements to
identify the likely cellular or tissue origins of cfDNA and to investigate
clonal dynamics during treatment.

## Planned cohort and current scope

The following choices define the current methylation study and its relationship
to the broader proposal. Cohort figures are planning estimates; actual enrolment
and usable sample counts must be recorded separately.

| Item | Current rule | Reason and source |
| --- | --- | --- |
| Patient cohort | Approximately 80 patients with diverse ALL subtypes. | Retain the conservative project estimate, as confirmed by the investigator on 21 September 2026. |
| Healthy plasma controls | Approximately 30 healthy paediatric plasma cfDNA samples remain a longer-term planning estimate. Healthy-plasma comparison and plasma-marker disease-specificity validation are on hold until suitable controls are available. | Planning estimate confirmed on 21 September 2026; investigator deferred this work on 27 October 2026 because controls are not expected soon. |
| Active compartments | Bone marrow (BM), peripheral blood (PB) and plasma; cerebrospinal fluid (CSF) work is on hold. | Investigator clarification on 21 September 2026. CSF and central nervous system (CNS) questions remain part of the broader proposal. |
| Current longitudinal scope | Induction, Days 1–29, using the specimen matrix below. | Existing project design; the proposal's consolidation, maintenance and later follow-up are broader study scope, with methylation implementation still to be specified. |
| Current sequencing stage | Genome-wide EM-seq pilot for method development and feasibility. | Existing production workflow and accepted pilot runs. |
| Planned sequencing stage | Targeted EM-seq covering candidate ALL signatures, healthy cell/tissue markers and protocadherin (PCDH) loci. | Proposal design; marker selection, assay depth and performance remain to be established before implementation. |

Subtype must therefore be retained as a core variable in the metadata and
carried through QC, deconvolution, ALL-versus-healthy comparisons and clonal
analysis. The first analyses should show patient-level trajectories and
subtype-stratified summaries rather than assume that every ALL subtype shares
one methylation profile. Pooling or formal modelling across subtypes will
require a prespecified strategy informed by the number of patients in each
subtype.

## Scientific aims

The numbering below follows the revised PhD proposal: MRD monitoring, cfDNA
origins, then integrated genetic and epigenetic clonal dynamics.

### Aim 1: Evaluate the added value of cfDNA for early MRD monitoring

The broader study will compare serial plasma tumour cfDNA measurements with
conventional BM MRD to describe early disease clearance. Patient-specific
markers intended to track the whole leukaemic population include genomic
fusion breakpoints and immunoglobulin/T-cell receptor rearrangements, measured
using droplet digital PCR (ddPCR). Their suitability for representing the
whole population must be documented for each patient.

The proposal also includes CSF cfDNA for investigating CNS involvement; that
component is on hold. The repository contributes methylation-derived source
estimates and candidate disease patterns for comparison with independent MRD
trajectories. Clinical added value has not been established by the current
pilot.

### Aim 2: Dissect the cellular origins of cfDNA during early treatment

The primary question is whether the lineage composition of cfDNA changes
between Day 1 and Day 29, particularly when total cfDNA remains abundant but
ddPCR-based MRD is low.

The intended first analysis is methylation-based lineage and tissue
deconvolution, which estimates contributions from different cell or tissue
sources using reference methylation profiles. Intermediate time points will
establish the trajectory rather than treating Day 1 and Day 29 as isolated
observations.

Candidate explanations include:

- continued shedding from residual ALL cells;
- normal haematopoietic cell turnover;
- inflammation or immune-cell turnover;
- treatment-related injury to non-haematopoietic tissues;
- a mixture of these sources.

Methylation deconvolution estimates source proportions. It does not by itself
demonstrate inflammation or causality. Interpretation of a systemic or
inflammatory effect should incorporate available clinical and laboratory
metadata, such as infection, fever, inflammatory markers, treatment phase and
organ toxicity. A wide range of clinical data is expected to be available, but
formal association with these variables is deferred from the first analysis to
keep the initial workflow focused.

### Aim 3: Define ALL methylation signatures and integrate genetic and epigenetic clonal dynamics

The project will identify candidate ALL-associated methylation signatures and
investigate whether methylation-defined molecule patterns change across BM, PB
and plasma cfDNA during treatment. Reference development and regional
methylation comparisons support signature discovery, with subtype and cell
composition retained in the interpretation. Genetic markers and independent
MRD measurements will support evaluation of whether these signals track
overall disease burden or distinguish changing subpopulations.

PCDH methylation patterns are candidate clonal barcodes. The whole-genome
pilot assesses recovery from diagnostic BM and plasma; targeted EM-seq is
planned to obtain deeper measurements. Independent genetic and MRD measurements
will support biological comparisons.

A CpG is a cytosine followed by guanine in the DNA sequence, where methylation
can be measured. Following a diagnostic pattern means measuring the frequency of the same
ten-CpG pattern at the same locus in matched specimens and later samples.
The [PCDH rationale](/path/to/hpc-home/brainstorming/pcdh_epiclone/PCDH_EPICLONE_WORKFLOW_RATIONALE.md#interpretation-and-planned-follow-up--clarified-2026-09-21)
defines the evidence needed for compartment correspondence, burden tracking
and subclone discrimination, including the limits of high diagnostic clonality
and the selected pilot benchmark. Targeted-assay sensitivity remains to be
specified and evaluated.

The intended evidence includes:

- persistence, contraction, expansion or replacement of PCDH patterns;
- concordance of patterns among matched BM, PB and cfDNA;
- changes in pattern abundance across Day 1, Day 15 and Day 29, with
  additional plasma observations on Days 4 and 8 where available;
- comparison with multiplex ddPCR trajectories;
- distinction of the ALL methylome from healthy haematopoietic and cfDNA
  backgrounds.

Until supported by longitudinal, genetic or orthogonal evidence, molecule
patterns will be described as methylation-defined subpopulations or candidate
epiclones, not validated clones.

## Longitudinal sampling design

This matrix describes the current induction-focused methylation study. The
CSF component is on hold. Later time points in the broader CIRCALL
proposal do not yet have a specified methylation workflow here.

Day 1 is the day treatment begins and is treated as the diagnostic time point.
Day-based nomenclature will be used for now.

| Day | Clinical interpretation | BM cellular DNA | PB cellular DNA | Plasma cfDNA | Multiplex ddPCR |
| --- | --- | ---: | ---: | ---: | ---: |
| Day 1 | Diagnosis and treatment initiation | Yes | Yes | Yes | Baseline status to be documented |
| Day 4 | Early treatment | No | No | Yes | Yes |
| Day 8 | Early treatment | No | No | Yes | Yes |
| Day 15 | Interim response | Yes | Yes | Yes | Yes |
| Day 29 | End-of-induction response | Yes | Yes | Yes | Yes |

This is the intended matrix for every patient. Actual availability and failures
must be represented explicitly; missing samples must not be replaced by
cross-patient comparisons.

## Molecular MRD measurements

Samples from Day 4 through Day 29 will undergo multiplex ddPCR. Patient-specific
targets may include:

- gene fusions;
- immunoglobulin or T-cell receptor (IG/TR) rearrangements;
- secondary alterations with prognostic value;
- one or more reference genes used to estimate an approximate variant allele
  fraction (VAF).

The ddPCR results provide a molecular disease trajectory against which total
cfDNA abundance, methylation-derived source fractions and candidate clonal
patterns can be compared. Approximate VAF is not automatically equivalent to
blast fraction or to the fraction of cfDNA molecules derived from ALL; copy
number, zygosity, assay sensitivity and the biological status of each target
must be retained in the interpretation.

Total cfDNA abundance is measured as ng DNA per mL plasma. This normalized
quantity, rather than total extraction yield or eluate concentration, should be
used for longitudinal abundance comparisons. The plasma volume and any
dilution, extraction and elution factors needed to reproduce the calculation
must be retained.

## Controls and reference material

The anticipated reference material follows the conservative planning estimate
in the scope table:

- approximately 30 healthy paediatric plasma cfDNA samples;
- a small number of healthy paediatric BM samples.

These are planning estimates, not an inventory of available controls.
**Healthy-plasma comparison and validation of plasma-marker disease specificity
are on hold until suitable healthy plasma cfDNA controls become available.**
The investigator deferred this work on 27 October 2026; no availability date
is assumed. Patient compartment comparisons and candidate-region exploration
can continue, with disease specificity left unresolved.

Public healthy white-blood-cell (WBC) profiles can provide a separate,
exploratory normal-blood baseline at candidate regions such as CDH4. They
do not fill the healthy-plasma control requirement. Cross-study differences
in age, sex, sequencing chemistry and processing must remain explicit; these
profiles are external descriptive references, not additional paired patients
in the existing DMRcate model.

Healthy paediatric cfDNA is the principal reference for normal cfDNA source
composition and variability. It should be matched as closely as possible for
age, sex and pre-analytical handling. Collection tube, processing delay,
centrifugation, storage, extraction and library batch can affect cfDNA results
and should be recorded.

The small healthy BM set can provide biological and technical context but is
not, by itself, a powered healthy BM reference cohort. Bulk BM also contains a
complex mixture of haematopoietic populations; comparisons must not be
interpreted as cell-intrinsic differences without composition assessment.

## Methylation references and deconvolution

The intended starting point for lineage deconvolution is the modified
UXM-based approach described by Nicola et al. (2026). That study:

- targeted 977 informative CpG segments;
- began from a published 40-cell-type methylation atlas;
- replaced broad B- and T-cell references with naïve and memory B- and T-cell
  populations;
- produced a modified 44-cell-type atlas;
- calculated segment methylation with UXM_deconv;
- used EpiDish robust partial correlations for deconvolution.

This is a candidate method, not yet a fixed project parameter. The published
implementation used a targeted panel and GRCh37 coordinates, whereas this
repository currently produces genome-wide GRCh38 EM-seq data. Before use, the
project must validate coordinate mapping, marker coverage, cell-type
identifiability and fit in paediatric cfDNA and ALL samples.

The ALL methylome study by Hetzel et al. (2022) is a key reference for
distinguishing ALL-associated methylation from healthy haematopoietic
methylation. It showed that ALL has an unusually highly methylated genome, with
CpG-island hypermethylation and variation among genetic subtypes. Its purified
healthy haematopoietic cells and ALL subtype profiles may support marker
selection and interpretation, but they are not automatically a validated cfDNA
deconvolution atlas.

Because ALL methylation varies by subtype, subtype must be carried in all
metadata and analyses. The planned cohort contains diverse ALL subtypes, so a
single universal ALL methylation component must not be assumed without
validation.

## Computational analysis scope

Shared EM-seq processing and QC support both methylation aims. Source
deconvolution with UXM or EpiDISH addresses Aim 2. Experimental B-cell precursor
ALL (BCP-ALL) reference development connects source estimation in Aim 2 with
disease-associated signatures in Aim 3. Regional methylation comparisons,
molecule-level PCDH patterns and the exploratory W-footprint analysis
support Aim 3. The [downstream index](downstream/README.md) maps these components
to their scientific questions and current workflow status.

Regional methylation analysis describes averages across genomic regions and
cannot by itself distinguish cell-intrinsic change from composition. UXM estimates mixtures relative to its
reference atlas and cannot report a lineage absent from that atlas directly.
PCDH pattern analysis describes molecule-level methylation patterns and tests
their potential for disease and subclone tracking. A distinct pattern is not
automatically a distinct clone; lineage claims require orthogonal evidence.
The W-footprint project describes the distribution of methylation fractions
across selected fluctuating CpGs, sites whose methylation can change over time.
Its compartment and blast-content comparisons do not identify individual
clones. These analyses provide complementary evidence at different levels.

Operational entry points and current development status belong to the root and
downstream workflow READMEs. Reference identities and provenance belong to the
reference registry. This document governs only scientific scope and
interpretation.

The [regional methylation workflow](workflows/dmrcate/README.md)
implements paired compartment comparisons. Historical DMRichR methods and
results remain in the [workflow archive](archive/downstream/dmrichr/README.md).
The [downstream index](downstream/README.md) owns workflow classification;
[RESULTS.md](RESULTS.md) links completed analyses.

EpiDISH with the official UXM atlas is the user's preferred method for blastic
samples. Normal-atlas fits describe reference-relative contributions; an
experimental BCP-ALL coefficient requires normal-control specificity and mixture
calibration before it can be interpreted as a blast fraction. The
[BCP-ALL project](/path/to/hpc-home/brainstorming/bcpall_atlas/README.md) and
[PCDH project](/path/to/hpc-home/brainstorming/pcdh_epiclone/README.md) own their current
validation status and available inputs.

## Library preparation assumptions

The intended DNA input is 30 ng for both cellular DNA and cfDNA.

- Cellular DNA will be sheared before EM-seq library preparation.
- cfDNA will not be sheared because it is already fragmented.

Input mass alone does not make the two specimen types technically equivalent.
Fragment-size distribution, recoverable genome equivalents, library complexity
and CpG span may differ. Library protocol version, actual input mass, extraction
method, fragment distribution, PCR cycles and sequencing batch should be
recorded per sample.

## EM-seq fixed-clipping decision

Production fixed end-clipping was selected from an empirical M-bias comparison
of the same eight paired-end EM-seq samples processed with the nf-core EM-seq
10-bp preset and with zero fixed clipping. Adapter and quality trimming remained
enabled in both runs. The zero-clipping run showed the strongest bias at the
start of R1, a smaller effect at the start of R2, and increasing non-CpG
methylation near the end of R2. The 10-bp preset removed the prominent edge
artifacts but retained approximately 88.4% of the cytosine observations from
the zero-clipping run.

On 2026-08-31, the selected compromise was 8 bp from the 5-prime end of R1,
3 bp from the 5-prime end of R2, 3 bp from the 3-prime end of R1, and 10 bp
from the 3-prime end of R2. These asymmetric values target the observed bias
while retaining more sequence than clipping 10 bp from every end. The
executable values are maintained in `config/params.emseq.integrated.yml` with
`skip_trimming_presets: true`. The resulting eight-sample run
`2026-08-31-emseq-integrated-8samples-trim-8-3-3-10-v1` completed and passed the technical
PAT-conversion gate in the [2026-09-04 QC review](archive/runs/2026-09-04_emseq_trim_8_3_3_10_qc_review.md).
This factual status correction was recorded on 2026-09-09. Residual broad
positional gradients and the R2 terminal effect remain interpretation limits;
technical acceptance does not establish downstream biological validity.

## Interpretation limits

- Total cfDNA abundance is not a direct measure of ALL burden.
- A short cfDNA half-life excludes physical persistence of diagnostic
  molecules, not continued release by residual disease.
- Bulk BM and PB methylation reflect cell mixtures as well as cell-intrinsic
  methylation.
- General tissue atlases may not distinguish ALL blasts from closely related
  normal B-cell progenitors.
- Deconvolution fractions depend on atlas completeness, marker coverage and
  assay compatibility.
- Approximate ddPCR VAF is influenced by target biology and copy number.
- PCDH patterns can nominate candidate clones but require longitudinal,
  genetic, phased, single-cell or other orthogonal support for lineage claims.
- Subtype-associated methylation can be mistaken for a universal ALL signal.
- Missing longitudinal specimens and low cfDNA molecule depth can create
  informative censoring.

## Decisions still required

1. Confirm which intended specimens and time points are available for each
   patient.
2. Define subtype groupings and the minimum sample count required for
   subtype-stratified summaries or models.
3. Evaluate the experimental hyperdiploid BCP-ALL column, including
   normal-control specificity and mixture calibration. Accepted PAT and
   normal-atlas workflow operation is documented in the downstream index.
4. Decide whether later UXM references should include other subtype-specific
   ALL components after the hyperdiploid experiment is reviewed.
5. Define the targeted PCDH assay and confirmatory analysis, including panel
   content, informative molecule depth, marker specificity and criteria for
   genetic/MRD concordance. The existing whole-genome pilot targets, filters
   and eligibility rules are already documented in the PCDH rationale; they
   do not establish targeted-assay sensitivity.
6. Define sequencing depth separately for cellular DNA and cfDNA using pilot
   complexity and marker-coverage results.
7. Define the analysis-specific statistical summaries for DMRichR, UXM and
   PCDH EpiClone before applying them to the longitudinal cohort.
8. Define a later plan for selecting and testing clinical correlates without
   expanding the first longitudinal analysis.

## Key references

- Hackett SF et al. *In vivo lineage tracing across human tissues using
  methylation barcodes in the protocadherin gene cluster*. bioRxiv (2026),
  preprint. <https://doi.org/10.64898/2026.02.23.707349>
- Ma X et al. *Rise and fall of subclones from diagnosis to relapse in
  pediatric B-acute lymphoblastic leukaemia*. Nature Communications (2015).
  <https://doi.org/10.1038/ncomms7604>
- Nicola PA et al. *Polyclonal selection of immune checkpoint mutations in
  thyroid autoimmunity*. Nature (2026).
  <https://www.nature.com/articles/s41586-026-10493-9>
- Hetzel S et al. *Acute lymphoblastic leukemia displays a distinct highly
  methylated genome*. Nature Cancer (2022). PMID: 35590059.
  <https://pubmed.ncbi.nlm.nih.gov/35590059/>
- UXM_deconv: <https://github.com/nloyfer/UXM_deconv>
- nf-core/methylseq 4.2.0: <https://nf-co.re/methylseq/4.2.0/>
