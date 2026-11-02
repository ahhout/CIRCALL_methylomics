# CIRCALL methylomics

Project configuration, anonymised samplesheets and metadata, and downstream workflow source for CIRCALL methylomics. Workflows include PAT, UXM, EpiDISH, DecoNFlow and DMRcate.

The preprocessing pipeline is maintained separately in [methylseq_AH](https://github.com/ahhout/methylseq_AH). This repository starts with a snapshot of the reorganised project; full earlier history remains in the original and local CIRCALL repositories.

## Running on Flamingo

See [execution instructions](docs/methylseq_execution.md). The supported methylseq entrypoint is `run_methylseq_integrated.sh`. Configuration records the pipeline revision, production reference, shared work cache and result collection root. The source retains HPC-specific absolute paths. Reference payloads, installed tools, runtime databases and results are external and are not distributed here. Recreate reference links from the project path configuration and reference documentation before use on another checkout.

New runs reserve unique output bundles; completed bundles cannot be reused by supported launchers. Historical runs continue through their original entrypoints. The bounded CIRCALL smoke test uses 10,000 paired reads and the production reference; successful production execution is not implied by fixture tests.

Apptainer must come from module `containers/apptainer/1.3.6`; the environment file also contains a Conda Apptainer package, which the launcher does not select.

`SOURCE_MANIFEST.json` identifies the source snapshot and file hashes. Anonymised study metadata and samplesheets are included with the owner's authorization; sequencing data and historical runtime records are not included.
