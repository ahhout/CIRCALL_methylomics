# Reference entry points

Large reference bundles remain outside the repository. The registry exposes
the supported alignment bundles and official human/CpG source through symbolic
links. The external legacy wgbstools dictionary and default, described below,
are additional historical resources without separate registry links. The small control-source
FASTAs under `emseq_controls/` remain regular tracked files.

## Reference choices

| Repository entry point | Contents | Use | Decision |
|---|---|---|---|
| `wgbstools_hg38` | wgbstools hg38 FASTA and CpG dictionary | UXM and PAT coordinate namespace; human source for derived alignment bundles | Canonical human reference because UXM must use this exact CpG namespace. |
| `wgbs_hg38_lambda` | exact wgbstools hg38 plus `phage_lambda`, with Bismark index | Public BCP-ALL sodium-bisulfite WGBS (`em_seq: false`) | Lambda was reported in the source libraries and permits conversion QC; it is excluded from PAT, UXM and atlas analyses. pUC19 is absent because it was not reported. |
| `emseq_hg38_lambda_puc19` | exact wgbstools hg38 plus `phage_lambda` and `plasmid_puc19c`, with Bismark index | Production EM-seq (`em_seq: true`) | Both controls are retained because production QC measures unmethylated lambda and methylated pUC19 independently. |
| `legacy_igenomes_hg38_lambda_puc19` | former NCBI/iGenomes-derived GRCh38 plus both controls | Provenance for already completed EM-seq runs only | Retained so historical outputs remain auditable; do not select it for new runs after migration. |
| `uxm_atlas_pat` | 207 primary reference cell-type hg38 PAT and CSI files from GSE186458 | Reference training data for the UXM normal cell-type atlas | Direct download approved 2026-09-16. Tracked via manifest; payload ignored by Git. |

The two current wgbstools-derived alignment bundles deliberately share the
same human sequence but differ in reported spike-in controls. Chemistry
remains a pipeline parameter: public WGBS uses `em_seq: false`, while
production EM-seq uses `em_seq: true`.

The archived GRCh38-only launcher instead selected nf-core's `GRCh38`
iGenomes bundle and did not include control contigs. Although hg38 and GRCh38
refer to the same assembly generation, these bundles are not assumed to be
sequence-identical. The exact wgbstools hg38 plus lambda and pUC19 bundle is
the sole supported production EM-seq reference; the iGenomes route is retained
only as historical provenance.

## Canonical human source

`wgbstools_hg38` currently points to:

```text
/path/to/hpc-home/.conda/envs/nextflow/share/wgbs_tools/references/hg38
```

Required source checksums:

```text
362592f9bfeda5e761db0ccf6d8fca9988e2a6a4e22071985e4c86e1a53f7b75  hg38.fa.gz
28180d18667dc0a6e66dd4250593e03a8bf10be50e9212aefd4e9282afea8590  CpG.bed.gz
```

The symlink isolates workflows from this environment-internal path. Retarget it
only to an identical, checksum-verified copy if the reference is later moved
to stable shared storage.

## WGBS bundle

`wgbs_hg38_lambda` points to the validated BCP-ALL Step 3 bundle. Its combined
FASTA checksum is:

```text
642b80dba142e23ffdab162d6a9646d261a6b5d66acfe2c1492f94ec611cdbcd  BismarkIndex/genome.fa
```

The bundle contains one lambda record and no pUC19 record. Its 12 Bismark
index-file checksums are recorded inside the target directory.

## Production EM-seq bundle

`emseq_hg38_lambda_puc19` points to the versioned production bundle built from
the same canonical human FASTA and `emseq_controls/lambda_puc19.fa`. Its target
contains a manifest plus portable FASTA and index checksum files. Production
configuration must use this repository symlink, not its external target.

## Compatibility and legacy dictionary — checked 2026-09-08

The current EM-seq alignment bundle has 457 contigs and the public WGBS bundle
456. They share all 455 exact wgbstools hg38 human sequences and the official
29,401,795-site CpG namespace. The older iGenomes-derived alignment bundle has
197 contigs; only 176 of the 455 official human sequences match exactly.
Assembly names alone do not establish PAT/atlas compatibility.

The external dictionary at
`/path/to/hpc-home/.conda/envs/nextflow/share/wgbs_tools/references/GRCh38_lambda_puc19`
contains 29,152,891 sites, 248,904 fewer than the official dictionary. Its
`references/default` alias still selects this legacy resource. Supported PAT
commands must explicitly pass `--genome hg38`; do not retarget the legacy
default or reinterpret old PAT indices as the official namespace. These are
[audit findings](../docs/audits/emseq_repository_inventory_20260908.md), not
newly regenerated references.

Public sodium-bisulfite WGBS uses `em_seq: false` and has no project requirement
for spike-ins or pUC19. Its existing lambda contig is optional conversion-QC
support; absent lambda calls are recorded as unavailable. Preserve that
reference and its approved WGBS clipping policy. The shared requirement with
EM-seq is exact human-sequence/CpG compatibility, not identical chemistry or
control content.

## Local control sources

`emseq_controls/lambda_puc19.fa` is the source for the two production controls:

```text
phage_lambda
plasmid_puc19c
```

Its SHA-256 is:

```text
f0d3c9c5cdbf9b52c54192e4db1031dd6753763808491d7e4251c66145657a84
```

These controls are alignment and QC contigs. They must never enter the human
CpG dictionary, PAT coordinate system, UXM marker matrices or atlas.

## UXM atlas reference PATs

`uxm_atlas_pat/` is the download target for the primary reference PAT and index
files from [GSE186458](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186458)
([Loyfer et al., *Nature* 2023](https://doi.org/10.1038/s41586-022-05580-6)):

- **Genome assembly**: GRCh38 / hg38 (`.hg38.pat.gz` and `.hg38.pat.gz.csi`).
- **Scope**: 207 distinct GSM records with 82 distinct manifest cell-type
  labels, excluding the later WBC and cfDNA controls. These labels are not a
  verified lineage mapping; exact membership in the official UXM training set
  has not been independently established.
- **Payload size**: 70.169 GiB total across 207 PAT and 207 CSI files.
- **Manifest**: `download_manifest.tsv` tracks all 207 samples, expected byte
  sizes, and direct HTTPS source URLs from NCBI GEO.
- **Verification**: exact expected file sizes only; neither content checksums
  nor PAT/index consistency are verified by this downloader.
- **Safety**: `.gitignore` excludes completed PAT/CSI files and their `.part`
  downloads. Existing final files with unexpected sizes cause errors and are
  preserved for explicit resolution. Dry-run mode does not write or delete files.
  Completed files are reused; partial downloads resume on resubmission.

The [26 October 2026 donor audit](../docs/audits/2026-10-26_12-33-uxm-reference-donor-audit-v2/audit.json)
matched all 205 published Supplementary Table S1 specimens to this collection
and verified the expected sizes of all 207 PAT/index pairs. These are specimens,
not 207 independent donors: the publication records 135 source-scoped PatientID
labels. The two extra cardiomyocyte specimens have prior numerical evidence
supporting exclusion from the shipped heart column. Subsequent complete PAT
streaming reproduced all 39 PAT-backed atlas columns at published precision,
including all 19 endothelial specimens, the four selected heart specimens and
all four ovary/endometrial specimens. The `Ovary-Ep` label includes one ovarian
and three endometrial specimens. Exact count-preparation and pooling evidence
is in the [matched-reference record](../downstream/workflows/deconflow/runs/2026-10-26_13-27-uxm-celfie-matched-reference-v1.json). Megakaryocyte inputs
originate outside this collection. The [resolution assessment](../docs/audits/2026-10-26_12-33-uxm-reference-donor-audit-v2/resolution.json)
records the outstanding library/pooling questions and evidence needed to settle
them. The existing l4/40-cell atlas and production membership files are unchanged.

The subsequent [public megakaryocyte count investigation](../docs/audits/2026-10-26_12-58-megakaryocyte-count-investigation-v1/README.md)
retrieved 71 BLUEPRINT libraries and reconstructed candidate donor count
profiles. Pooling all libraries by donor reproduces all 124 non-missing
percentages and all 11 missing values at 45 published CpGs after verified
hg19-to-hg38 coordinate reconciliation. The single-cell-only candidate fails
for donor 1. These checks support a documented CELFIE megakaryocyte count
component, but do not uniquely identify historical library membership or
recover UXM fragment patterns. The candidate covers 971 of the existing 979
regions in at least one donor, including all 25 megakaryocyte target regions;
the eight uncovered regions remain explicit. The count component is retained
on scratch and has not replaced any production reference. Subsequent patient
pilot checks found concrete incompatibilities between these counts and the
shipped megakaryocyte UXM profile at checked regions. Matching the published
donor fingerprints therefore does not establish the required atlas identity;
see the [comparison record](../downstream/workflows/deconflow/ACTIVE_COMPARISON.md).

The [27 October source reconciliation](../docs/audits/2026-10-27_08-17-megakaryocyte-source-reconciliation-v1/README.md)
identifies those 71 libraries as Farlik et al. 2016, GEO GSE87196,
GSM2324744–GSM2324814: 61 single-cell libraries and ten five-cell pools from
three donors. The publication chain supports a common original BLUEPRINT
study, but exact historical processing remains unresolved. The 45-CpG
fingerprint check spans eight small regions and excludes the two conflicting
markers; several library subsets share those fingerprints. Every retained
library has zero unmethylated observations at both conflicts, so changing pool
membership or nonnegative weights cannot restore positive strict regional
unmethylated-fragment counts. Published hg38 marker-selection metadata also
expects low methylation at these loci. No reference or coordinate convention
was changed, and a possible hg19/hg38 CpG-index mismatch remains an untested
hypothesis rather than an identified error.

### Download and verification

Submit the background download job via Slurm:

```bash
sbatch references/scripts/run_download_uxm_atlas_pat.sh
```

To run a dry-run or verify existing files interactively:

```bash
python3 references/scripts/download_uxm_atlas_pat.py --dry-run
python3 references/scripts/download_uxm_atlas_pat.py --verify-only
```

On 2026-09-16, initial job `5735174` was cancelled at the user's request
following review. Corrected launcher resubmitted as job `5735215`, using
`references/uxm_atlas_pat/download_manifest.tsv` and the same output directory,
`references/uxm_atlas_pat/`. Existing completed files are reused by size and
`.part` transfers resume with curl. This is a standalone download, not a
methylseq execution; pipeline version, params, config, samplesheet and Nextflow
`-resume` are not applicable. Submission does not establish download completion.
