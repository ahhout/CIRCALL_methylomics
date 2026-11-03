# Successful CIRCALL smoke execution

Verified on 3 November 2026 against Slurm accounting, the frozen attempt end
record, output ownership, task trace and report files. Job **6251849** completed
with exit code **0:0** on **2 November 2026 at 16:14:10 UTC** (17:14:10 Paris).
Slurm elapsed time was **01:37:14**. All **12 Nextflow tasks** completed with
exit code 0; the output owner and attempt end record both mark completion.

The input was the first **10,000 paired records** of the first production
samplesheet entry. Its identifier and production scientific settings were
retained. This validates bounded execution from the new project/pipeline layout;
it does not establish full-cohort production performance or biological QC
acceptance. Historical resume and cross-layout cache reuse remain unverified.

| Setting | Frozen execution value |
| --- | --- |
| Pipeline | nf-core/methylseq 4.2.0 |
| Pipeline commit | `c6a960d78ce8316e26a42ed4c0291408114695ea` |
| Project commit | `7e88608ac65ab491fbbb728844cfc10461de06c7` |
| Runtime | Nextflow 25.10.2; module Apptainer 1.3.6; Slurm |
| Reference | Exact production wgbstools hg38 plus lambda and pUC19 |
| Params | `config/params.emseq.integrated.yml` |
| Config | `config/nextflow.config` plus `config/methylseq_smoke_resources.config` |
| Input sheet | Frozen `pipeline_info/submitted/samplesheet.csv` inside the output below |
| Resume | No; new output reservation and session |
| Run name | `2026-11-02_14-32-emseq-integrated-layout-smoke-v2` |

Output directory:

```text
/path/to/hpc-scratch/nextflow/nf-core-methylseq/results/bulk/2026-11-02_14-32-emseq-integrated-layout-smoke-v2
```

The completed stages include FastQC, Trim Galore, Bismark alignment and
deduplication, methylation extraction, sorting/indexing, Qualimap, cytosine
reporting, Bismark reports and MultiQC. The launcher then wrote
`control_qc/control_qc.tsv`. The MultiQC report is
`multiqc/bismark/multiqc_report.html`. Report existence is an execution check;
no clinical interpretation or cohort-wide QC acceptance is asserted here.

Frozen identity, submitted params/config/input sheets and the exact command are
under `pipeline_info/submitted/`. Attempt `6251849-f807a2b87bbc` retains
`started.json`, `ended.json`, `trace.txt`, `nextflow.log`, reports and control-QC
completion. The run name's 14:32 UTC timestamp is naming time; the launcher
actually started at 14:37:16 UTC. These times must not be conflated.

First attempt 6251800 failed before analysis because an empty Nextflow config
argument became `true`. The explicit boolean fix was tested, the failed output
was retained, and this retry used a fresh name without resume. Neither completed
output nor its source identity was changed by the subsequent housekeeping.

[Machine-readable validation and evidence hashes](methylseq_smoke_validation.json)
record the exact paths, timestamps, source revisions and input hashes.
