# CIRCALL methylomics

Project repository for CIRCALL methylomics, including PAT, UXM, EpiDISH, DecoNFlow and DMRcate workflows.

Publication is in progress. The repository currently contains this README, `.gitignore`, and the approved cohort tables:

- [samplesheets/samplesheet.csv](samplesheets/samplesheet.csv): all 12 samples, consolidating the former 4-sample and 8-sample sheets. Original sample identifiers, paired FASTQ paths, CSV schema and the first sample are preserved. The separate 8-sample sheet is retired.
- [metadata/metadata.tsv](metadata/metadata.tsv): the matching 12 sample records, published with investigator approval.

FASTQ paths refer to the HPC environment; sequencing files are not included. `metadata/comparisons.tsv` is retained locally for an archived Stage 4 pilot and is excluded from publication. No supported workflow consumes that table.

The remaining project source snapshot has not yet been uploaded. The preprocessing pipeline is maintained separately in [methylseq_AH](https://github.com/ahhout/methylseq_AH). Earlier project history and existing results remain in their original HPC locations.
