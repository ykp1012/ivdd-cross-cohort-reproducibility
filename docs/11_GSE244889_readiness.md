# GSE244889 NP Data-Readiness Report

Audit date: 2026-08-14.

## Scope and cohort boundary

GSE244889 is a human nucleus pulposus (NP) study containing two distinct
assays.  This audit keeps them separate:

- seven scRNA-seq libraries with GEO-supplied 10x Cell Ranger v6.1.1 feature,
  barcode, and Matrix Market triples;
- six bulk RNA-seq samples represented only by a GEO-supplied processed FPKM
  table.

The scRNA component has four GEO-labelled MDD samples (grades 1, 1, 2, 2) and
three SDD samples (grades 3, 4, 4).  The SOFT titles encode ages 17, 55, 24,
30 years in the MDD group and 41, 62, 59 years in the SDD group; sex is also
title-encoded (two female/two male per group).  The age distributions overlap
but are not balanced, and cohort membership is still completely aligned with
the grade/status contrast.  This small observational design cannot identify
an age-independent degeneration effect.

The sample-title prefix is unique within the seven scRNA records and is used
as a **presumed donor/library key**.  GEO does not provide a separately named
patient identifier, so cells remain nested within that key and must never be
treated as independent biological replicates.

## Acquisition and integrity

All source files were retained unchanged in `data/raw/`.  The detailed
machine-readable provenance table is `data/derived/GSE244889_asset_ledger.csv`.

| Asset | Bytes | SHA-256 | Role |
|---|---:|---|---|
| `GSE244889_RAW.tar` | 546,119,680 | `3bf1617b85e2fa08e2074c0de8fc87aaeb97c42b29bbbce77634f2c43d35c9b1` | raw scRNA 10x archive |
| `GSE244889_family.soft.gz` | 3,329 | `877bb67d2d22e35e28b047857a08999c73355d6a69b733141d2f7e2b81f0bdcb` | GEO sample metadata |
| `GSE244889_FPKM.txt.gz` | 647,864 | `6fcccafc760c5b70e55c322a07c94ff5493e36cbb0c50cca9bd1a7ee209c7ea4` | processed bulk expression table |
| `GSE244889_filelist.txt` | 1,564 | `df08d97fc3417ce809927c65437356cfb1cbd9e562ff57d2b1792ba9f8883635` | GEO archive-member list |

The GEO `filelist.txt` and the downloaded TAR agree exactly: archive name and
byte count match, and all 21 listed file members have matching names and
uncompressed member sizes.  The detailed cross-check is
`data/derived/GSE244889_filelist_audit.csv`.

## scRNA archive structure

The TAR contains exactly 21 files: one barcode, one feature, and one Matrix
Market file for each of seven expected scRNA GSM records.  Every library has
37,487 feature rows; barcode counts range from 5,437 to 11,562 (55,264 total
pre-QC barcodes).  The 10x dimension audit confirms every feature and barcode
file agrees with its Matrix Market header.  The full stream audit then read
each compressed Matrix Market member to EOF, exercising gzip stream integrity,
and found that every observed coordinate-line count equals the declared NNZ.

| Audit | Result |
|---|---|
| Expected scRNA libraries represented | 7/7 |
| Complete 10x triples | 7/7 |
| Feature/barcode versus Matrix Market dimension checks | 7/7 pass |
| Full gzip stream/NNZ-line integrity checks | 7/7 pass |
| Feature identifiers per library | 37,487 |
| Total pre-QC barcodes | 55,264 |

The corresponding tables are:

- `data/derived/GSE244889_scrna_sample_ledger.csv`;
- `data/derived/GSE244889_archive_inventory.csv`;
- `data/derived/GSE244889_10x_matrix_audit.csv`;
- `data/derived/GSE244889_10x_stream_integrity_audit.csv`.

Thus, the scRNA component is eligible for independent, donor-aware raw-count
ingestion and later pseudobulk construction **only after** the project's
pre-specified QC, resident-cell annotation, exclusion-ledger, and cell-count
sensitivity gates have been applied.  This audit is not a biological result.

## Bulk FPKM asset is separate

The bulk table has 19,834 gene rows and six sample columns: three MDD
(grade 2) and three SDD (grades 3-4).  Its complete structural audit finds no
row-width mismatch, empty value, non-numeric value, or negative value.  It
contains 19,832 unique gene labels; `1-Mar` and `2-Mar` are duplicated labels
and must be handled by a pre-specified duplicate-gene policy before any
score-level use.  FPKM header columns were mapped to the six bulk SOFT titles
only through their shared `S11`-`S16` tokens; this is a transparent presumed
mapping rather than an explicit table-provided GSM mapping.

The bulk table is a processed FPKM asset, not raw counts and not single-cell
data.  It must not be merged with the 10x matrices, used for scRNA pseudobulk,
or supplied to edgeR/DESeq2 raw-count models.  It may later serve only as a
separately normalized, donor-level direction check after module scoring and
duplicate-gene handling rules are frozen.  Details are in
`data/derived/GSE244889_bulk_sample_ledger.csv` and
`data/derived/GSE244889_FPKM_audit.csv`.

## Readiness decision

GSE244889 clears provenance and raw-archive structural gates for the
seven-library NP scRNA component.  It is an independent small external NP
severity cohort with a 4-versus-3 donor contrast.  Because the comparison has
fewer than four severe donors, it cannot provide a confirmatory significance
claim under the project protocol.  Its permitted role is a donor-level,
directional reproducibility and leave-one-donor-out stability check after the
analysis plan is frozen.  Any results must state the small group size,
title-based presumed donor keys, age/grade structure, and the separate status
of the bulk FPKM table.

## Reproduction commands

```powershell
& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\parse_geo_soft_metadata.py" `
  ".\data\raw\GSE244889_family.soft.gz" `
  --dataset GSE244889 `
  --output ".\data\derived\GSE244889_soft_metadata.csv"

& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\audit_gse244889_assets.py" `
  ".\data\derived\GSE244889_soft_metadata.csv" `
  ".\data\raw\GSE244889_family.soft.gz" `
  ".\data\raw\GSE244889_RAW.tar" `
  ".\data\raw\GSE244889_FPKM.txt.gz" `
  ".\data\raw\GSE244889_filelist.txt" `
  --asset-output ".\data\derived\GSE244889_asset_ledger.csv" `
  --scrna-ledger-output ".\data\derived\GSE244889_scrna_sample_ledger.csv" `
  --bulk-ledger-output ".\data\derived\GSE244889_bulk_sample_ledger.csv" `
  --fpkm-audit-output ".\data\derived\GSE244889_FPKM_audit.csv" `
  --filelist-audit-output ".\data\derived\GSE244889_filelist_audit.csv"

& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\inspect_geo_archives.py" `
  ".\data\raw\GSE244889_RAW.tar" `
  --output ".\data\derived\GSE244889_archive_inventory.csv"

& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\audit_10x_tar.py" `
  ".\data\raw\GSE244889_RAW.tar" `
  ".\data\derived\GSE244889_scrna_sample_ledger.csv" `
  --output ".\data\derived\GSE244889_10x_matrix_audit.csv"

& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\verify_10x_tar_matrix_streams.py" `
  ".\data\raw\GSE244889_RAW.tar" `
  ".\data\derived\GSE244889_10x_matrix_audit.csv" `
  --output ".\data\derived\GSE244889_10x_stream_integrity_audit.csv"
```
