# GSE251686 External NP Readiness Report

Audit date: 2026-08-14.

## Scope and decision

GSE251686 is a human nucleus-pulposus (NP) single-cell RNA-seq series.  The
outer archive, TAR/GZIP structure, and five nested matrix payloads pass the
provenance and structural checks, but `GSM7986002` fails the independent text
integrity audit.  It is retained in the ledger for traceability and is
**excluded from all downstream analysis**; it must not be repaired, have NUL
bytes deleted, be interpolated, or be otherwise transformed for use.  The
remaining usable records are two mild and three severe presumed
sample/library keys.  With no independently named patient identifier and no
per-sample demographic or individual Pfirrmann-grade fields, the remaining
cohort is only an incomplete, non-balanced exploratory NP mild-versus-severe
**direction check** after the scoring rule and cell-selection rules are frozen
elsewhere.

It is not eligible for a confirmatory p-value claim, discovery of a signature,
cell-as-replicate inference, external validation, decisive meta-analysis, or
an age-independent degeneration conclusion.  It cannot support
normal-versus-degenerated claims: both groups comprise surgical
lumbar-disc-herniation NP material, and the linked article explicitly says
that normal NP was not included.

## GEO/SOFT design and sample mapping

GEO lists six human NP samples under `GPL24676` (Illumina NovaSeq 6000).  The
GEO sample titles define the only available severity groups:

| GEO title group | GSMs | Archive labels | Permitted analysis key |
|---|---|---|---|
| Mildly degeneration, replicate 1 | `GSM7986001` | `NP1` | included presumed sample/library key |
| Mildly degeneration, replicate 2 | `GSM7986002` | `NP3` | **excluded**: malformed Matrix Market text; do not repair or use |
| Mildly degeneration, replicate 3 | `GSM7986003` | `NP4` | included presumed sample/library key |
| Severely degeneration, replicate 1-3 | `GSM7986004`-`GSM7986006` | `NP5`, `NP6`, `NP9` | three included presumed sample/library keys |

All six SOFT records and linked BioSample records label the source as NP.  The
GEO series design states that samples were obtained from lumbar-disc-herniation
patients undergoing endoscopic discectomy at L4/L5 or L5/S1, but it does not
link a particular sample to either level.  Neither GEO SOFT nor the BioSample
records expose patient ID, age, sex, individual Pfirrmann grade, batch, or
clinical covariates.  Consequently, `GSM7986001` through `GSM7986006` are
**presumed sample/library keys**, not independently verified biological-donor
identities.  Only five keys pass the matrix text-integrity gate; cells remain
nested within those keys and must never be used as independent replicates.

The linked study reports six patients with Pfirrmann II, III, or IV and calls
the two title-defined groups mild and severe, but it does not provide an
auditable GSM-to-individual-grade mapping in the GEO/SOFT record.  This project
therefore retains only the published mild/severe title grouping and does not
impute grade labels from the article.

The machine-readable sample ledger is
`data/derived/GSE251686_donor_ledger.csv`; the corresponding summary is
`data/derived/GSE251686_design_summary.csv`.

## Acquisition and filelist integrity

The untouched GEO archive and SOFT file are kept under `data/raw/`.  The
retrieval records are appended to `data/raw/download_manifest.ndjson`.

| Asset | GEO URL | Bytes | SHA-256 |
|---|---|---:|---|
| `GSE251686_RAW.tar` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE251nnn/GSE251686/suppl/GSE251686_RAW.tar` | 287,856,640 | `7b95a5eec4aa249ebdc3eb93caa5a7f9b0ab0c4cebceaa440cfe4fdb3a31a4b4` |
| `GSE251686_family.soft.gz` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE251nnn/GSE251686/soft/GSE251686_family.soft.gz` | 2,498 | `a3248b36f0f56746e4cc0315abeefb851c97bcc86f4f0d8798c120360cdd4f23` |
| `GSE251686_filelist.txt` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE251nnn/GSE251686/suppl/filelist.txt` | 580 | `69cd8d9e1020923370e59dd05fbadf0524bba1a0442cc0b45c331f9fea6dc406` |

The GEO filelist contains one outer archive and six per-GSM nested archives.
All seven published names and byte sizes exactly match the downloaded outer
TAR and its members.  The audit is recorded in
`data/derived/GSE251686_filelist_audit.csv`; the unmodified member inventory
is `data/derived/GSE251686_archive_inventory.csv`.

## Matrix structure and identifier audit

The outer TAR contains six nested `*.tar.gz` payloads.  Each contains exactly
one `genes.tsv`, `barcodes.tsv`, and `matrix.mtx`; this is a matrix-exchange
layout, not the one-level Cell Ranger feature/barcode/MTX archive layout used
by the GSE230809 discovery project.  Dimensions match for all six libraries,
but the independent stream audit finds a single malformed payload in
`GSM7986002`.

| GSM | Group | Features | Pre-QC barcodes | Matrix nonzeros | Stream-audit status |
|---|---|---:|---:|---:|---|
| GSM7986001 | mild | 27,770 | 15,448 | 19,176,437 | pass |
| GSM7986002 | mild | 31,158 | 12,275 | 30,130,524 | **excluded**: 5,471,800 NUL bytes, one malformed line, text-integrity failure |
| GSM7986003 | mild | 27,567 | 5,932 | 11,399,727 | pass |
| GSM7986004 | severe | 30,035 | 7,356 | 19,591,322 | pass |
| GSM7986005 | severe | 20,847 | 2,433 | 3,626,556 | pass |
| GSM7986006 | severe | 58,336 | 7,995 | 15,162,098 | pass |

The total number of supplied barcodes is 51,439 across all six GEO records;
39,164 barcodes remain in the five included records.  These are library-yield
observations only, not biological cell-composition results.  The linked paper
reports 47,610 retained cells after its own filtering.  The differing counts
show that the public nested matrices are not simply the paper's final
retained-cell matrix; they are compatible with a post-cell-calling,
pre-paper-QC matrix set.  The five included records must undergo the project's
donor-stratified QC before any scoring.

Within every library, `genes.tsv` has two columns, Ensembl feature IDs are
unique, and barcodes are unique.  Gene symbols are unique in five libraries;
`GSM7986006` has 1,523 duplicate gene symbols across 58,336 unique Ensembl
features.  Later cross-library feature mapping must therefore use Ensembl ID
as the primary feature key, followed by a pre-specified version-stripping and
symbol-mapping rule; it must not silently collapse rows by symbol.  The
machine-readable checks are:

- `data/derived/GSE251686_nested_10x_matrix_audit.csv`
- `data/derived/GSE251686_matrix_identifier_audit.csv`
- `data/derived/GSE251686_nested_matrix_stream_audit.csv`
- `data/derived/GSE251686_GSM7986002_header_audit.csv`
- `data/derived/GSE251686_GSM7986002_header_audit.json`

## Platform and processing evidence conflict

The most specific provenance source is GEO/SOFT.  For every GSM it says:

- tissue dissociation used the Singleron PythoN automated dissociator and
  sCelLive reagents;
- libraries used the GEXSCOPE Single-Cell Sequencing Kit and Singleron Matrix
  automated single-cell system, with microwell loading;
- sequencing occurred on Illumina NovaSeq 6000;
- CeleScope v1.10.0 extracted/corrected barcodes and UMIs, then reads were
  aligned to GRCh38 with STAR and assigned with FeatureCounts;
- a gene-by-cell matrix was generated by barcode, UMI, and gene assignment.

The linked article's Methods section reports the same GEXSCOPE/Singleron and
CeleScope workflow.  Its Results section, however, describes the profiling as
performed on a "10X Genomics platform."  The outer member names also contain
`EmptyDrops_CR_matrix`.  These labels conflict with the detailed GEO/SOFT and
Methods workflow.  For this project, GEO/SOFT controls the data-provenance
classification: GSE251686 is described as a **Singleron GEXSCOPE/CeleScope
UMI-derived matrix dataset**.  It must not be called Cell Ranger output, a
10x Genomics dataset, or an unfiltered raw-read dataset without independent
sequencing-read reprocessing.

The supplied matrix is integer, UMI-derived, but it reflects the contributor's
CeleScope/EmptyDrops pipeline and is not raw FASTQ data.  Its eligibility for
the project's primary count-model path is therefore set to `False` pending a
separate pre-specified technical decision and full reproducibility assessment;
the present audit does not turn it into raw-count negative-binomial input.

## Analysis boundary and next gate

After the project's pre-specified modules, Ensembl mapping, QC thresholds,
resident-cell annotation, and score rule were frozen without reference to this
cohort's group results. The five included records now have a separately
generated audit-gated score table under
`data/derived/GSE251686_exploratory_scores/`; that directory is intentionally
isolated from the default `data/derived/donor_module_effect_summary/` result.
They may supply:

- per-presumed-sample QC and retained-cell reporting;
- a 2-versus-3, non-balanced presumed-key descriptive effect size and
  score-direction display;
- leave-one-key-out direction stability as a small-sample sensitivity check.

It may not supply confirmatory differential expression, an external validation
claim, a universal pathway claim, age-adjusted inference, grade-dose
inference, causal interpretation, or a clinical/treatment claim.  Any result
must use association language and report all six GEO records, the five
included keys, the explicit exclusion of `GSM7986002`, pre/post-QC cell counts,
identifier mapping loss, and the platform/metadata limitations above.

The separate GSE251686 score/effect package must never be added to the default
20-effect table, the default cross-cohort sign-alignment count, or the main
figures/tables without a newly pre-specified analysis decision and a full
re-audit. Its present role is an isolated descriptive sensitivity display only.

## Reproduction commands

Run the commands below from the project root.  The scorer's
`--min-retained-cells 30` primary eligibility threshold is 30 source-restricted
cells per presumed sample/library key.  The subsequent isolated effect
summarizer deliberately applies a stricter audit gate: every selected key must
also pass the 20-, 30-, and 50-cell thresholds (including the `>=50` gate)
before the 2-versus-3 descriptive contrast is emitted.  All five retained keys
pass the 50-cell gate in the recorded run.

```powershell
& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\inspect_geo_archives.py" `
  ".\data\raw\GSE251686_RAW.tar" `
  --output ".\data\derived\GSE251686_archive_inventory.csv"

& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\parse_geo_soft_metadata.py" `
  ".\data\raw\GSE251686_family.soft.gz" `
  --dataset GSE251686 `
  --output ".\data\derived\GSE251686_soft_metadata.csv"

& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\audit_gse251686_filelist.py" `
  ".\data\raw\GSE251686_filelist.txt" `
  ".\data\raw\GSE251686_RAW.tar" `
  --output ".\data\derived\GSE251686_filelist_audit.csv"

& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\audit_nested_10x_tar.py" `
  ".\data\raw\GSE251686_RAW.tar" `
  --output ".\data\derived\GSE251686_nested_10x_matrix_audit.csv"

& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\audit_gse251686_matrix_identifiers.py" `
  ".\data\raw\GSE251686_RAW.tar" `
  --output ".\data\derived\GSE251686_matrix_identifier_audit.csv"

& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\audit_nested_matrix_stream.py" `
  ".\data\raw\GSE251686_RAW.tar" `
  --output ".\data\derived\GSE251686_nested_matrix_stream_audit.csv" `
  --json ".\data\derived\GSE251686_nested_matrix_stream_audit.json"

& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\make_gse251686_ledger.py" `
  ".\data\derived\GSE251686_soft_metadata.csv" `
  ".\data\derived\GSE251686_nested_10x_matrix_audit.csv" `
  --stream-audit ".\data\derived\GSE251686_nested_matrix_stream_audit.csv" `
  --ledger-output ".\data\derived\GSE251686_donor_ledger.csv" `
  --summary-output ".\data\derived\GSE251686_design_summary.csv"

& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\score_gse251686_exploratory.py" `
  ".\data\raw\GSE251686_RAW.tar" `
  ".\data\derived\GSE251686_donor_ledger.csv" `
  ".\config\program_modules.json" `
  --stream-audit ".\data\derived\GSE251686_nested_matrix_stream_audit.csv" `
  --nested-audit ".\data\derived\GSE251686_nested_10x_matrix_audit.csv" `
  --identifier-audit ".\data\derived\GSE251686_matrix_identifier_audit.csv" `
  --panel-config ".\config\cell_marker_panels.json" `
  --output-dir ".\data\derived\GSE251686_exploratory_scores" `
  --min-genes 200 `
  --min-umi 500 `
  --max-mt-pct 20.0 `
  --min-retained-cells 30 `
  --min-mapped-fraction 0.80

& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\summarize_gse251686_exploratory.py" `
  --score-table ".\data\derived\GSE251686_exploratory_scores\GSE251686_exploratory_module_scores.csv" `
  --library-ledger ".\data\derived\GSE251686_exploratory_scores\GSE251686_exploratory_library_ledger.csv" `
  --score-parameters ".\data\derived\GSE251686_exploratory_scores\GSE251686_exploratory_score_parameters.csv" `
  --output-dir ".\data\derived\GSE251686_exploratory_scores" `
  --confidence-level 0.95 `
  --bootstrap-replicates 10000 `
  --random-seed 20260814 `
  --min-mapped-fraction 0.80
```
