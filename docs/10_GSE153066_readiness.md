# GSE153066 External NP Readiness Report

Audit date: 2026-08-14.

## Acquisition and integrity

The usable compact expression asset was downloaded from the GEO FTP location;
no SRA files were downloaded and the original compressed file was retained
unchanged under `data/raw/`.

| Asset | URL | Bytes | SHA-256 |
|---|---|---:|---|
| `GSE153066_AllSample.counts.tsv.gz` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE153nnn/GSE153066/suppl/GSE153066_AllSample.counts.tsv.gz` | 100,366,013 | `ca081538d81d671bafe4af4e1e161a8a8140fc3d9d35e44c66f42cade24decb1` |
| `GSE153066_family.soft.gz` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE153nnn/GSE153066/soft/GSE153066_family.soft.gz` | 2,793 | `974b09e054b3d3dcbb42ab932cc0bca23a7f09ae35058599be7787f77d1401b0` |

Both records are appended to `data/raw/download_manifest.ndjson`.  The
matrix passed a complete gzip stream read.

## Matrix structure audit

- Format: gzip-compressed dense TSV, genes in rows and cells in columns.
- Header first field: `gene`.
- Matrix dimensions: 36,847 gene rows by 85,205 cell columns.
- All 36,847 gene identifiers are unique.
- All 85,205 cell-barcode headers are unique.
- Every non-empty gene row has exactly 85,205 serialized value fields.
- Barcode prefixes are exactly `CTL1`-`CTL8` and `IDD1`-`IDD8` (16 prefixes).
- The prefix counts are recorded in `data/derived/GSE153066_donor_ledger.csv`;
  total cells are 85,205 (1,750-9,976 cells per prefix).

The value-type check is deliberately sampled, not a claim about every matrix
entry: the first and every 500th non-empty gene row were checked, using 258
evenly distributed columns per sampled row (74 rows and 19,092 serialized
values in total).  Every sampled value was a non-negative integer.  The
script does not materialize the dense matrix.

The machine-readable structure result is
`data/derived/GSE153066_count_matrix_audit.csv`.

## Barcode-to-sample and donor mapping

Each observed prefix exactly matches one unique GEO `Sample_title` in the SOFT
record and one GSM.  The dedicated audit therefore passes the
barcode-to-sample mapping check for all 16 samples.  All samples are labelled
nucleus pulposus (NP):

- relatively normal: 8 presumed sample-level donors (`CTL1`-`CTL8`);
- degenerated: 8 presumed sample-level donors (`IDD1`-`IDD8`).

GEO does not expose a separately named patient identifier in this record.
Consequently, the sample title/prefix is a **presumed donor and library key**,
not independently verified patient identity.  Cells must remain nested within
that key; they cannot be treated as independent biological replicates.

## Processing and interpretation limits

GEO states that UMI-tools was used to calculate raw counts and a cell
whitelist, followed by Seurat normalization and cell filtering (`MT% < 20%`,
detected genes between 200 and 5,000).  The compact matrix is therefore a
GEO-retained, processed-cell matrix rather than an untouched raw-read archive.
It may support donor-level count aggregation only with this prior-filtering
history reported; it does not support claims of a newly reprocessed raw-read
pipeline or recovery of excluded cells.  Per-cell QC metrics cannot be
recomputed from the public asset alone.

The public metadata report relatively normal samples from motor-vehicle
accident material and degenerated samples from lumbar-disc-herniation
material.  Clinical source is consequently completely confounded with the
normal/degeneration label.  Ages are 32-40 in the relatively normal group and
38-47 in the degenerated group, so age remains a relevant cohort covariate;
the public design does not justify an age-independent disease interpretation.

## Readiness decision

GSE153066 passes the archive, matrix-structure, and barcode-to-SOFT mapping
audits and is suitable as an independent human NP normal-versus-degenerated
support cohort.  The usable contrast is donor-level 8 versus 8, with the
sample-level donor assumption and clinical-source/age limitations above.
Results must be reported as association or directional reproducibility.  No
causal, age-independent, clinical-prediction, or treatment conclusion is
permitted from this cohort alone.

## Reproduction command

```powershell
& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\audit_gse153066_count_matrix.py" `
  ".\data\raw\GSE153066_AllSample.counts.tsv.gz" `
  ".\data\raw\GSE153066_family.soft.gz" `
  --matrix-audit-output ".\data\derived\GSE153066_count_matrix_audit.csv" `
  --sample-ledger-output ".\data\derived\GSE153066_donor_ledger.csv"
```
