# Data Provenance Ledger

## Acquisition policy

- Download only from the original GEO FTP/HTTPS location into `data/raw/`.
- Keep archive names unchanged.
- Record URL, retrieval time, byte size, SHA-256 checksum, and archive listing.
- Do not alter raw archives. Extraction happens under `data/derived/`.

## Planned data assets

| Accession | Resource role | Species | Assay | GEO-reported design | Raw archive expected | Status |
|---|---|---|---|---|---|---|
| GSE229711 | Exploratory discovery, low-grade AF/NP | Human | scRNA-seq | Three young male donors, paired AF and NP; ages 21-27 | `GSE229711_RAW.tar` | downloaded; archive and 10x structure audited |
| GSE230808 | Exploratory discovery, advanced-degeneration AF/NP | Human | scRNA-seq | Ten male donors, AF n = 10 and NP n = 8; ages 37-73; eight paired AF/NP donors | `GSE230808_RAW.tar` | downloaded; archive and 10x structure audited |
| GSE165722 | NP score-level severity direction support | Human | scRNA-seq | Eight NP donors; source publication reports grades II-V, two donors per grade; GEO SOFT labels I-IV and is not used as the clinical-grade authority | `GSE165722_RAW.tar` | downloaded and audited |
| GSE244889 | NP severity direction check | Human | scRNA-seq plus bulk | Seven scRNA-seq NP donor/library keys, mild n = 4 and severe n = 3; a separate six-sample bulk FPKM matrix | `GSE244889_RAW.tar`, `GSE244889_FPKM.txt.gz` | downloaded; provenance, archive, 10x, and FPKM structure audited |
| GSE153066 | NP external normal-versus-degenerated support | Human | scRNA-seq | 16 NP samples, relatively normal n = 8 and degenerated n = 8; sample-prefixed barcodes map cells to samples | `GSE153066_AllSample.counts.tsv.gz` | downloaded; matrix structure, sample-prefix mapping, and metadata audited |
| GSE251686 | Isolated NP exploratory severity check | Human | scRNA-seq | Six GEO records; after stream-integrity audit, mild n = 2 and severe n = 3 are usable; GSM7986002 is excluded | `GSE251686_RAW.tar` | downloaded; nested archive/stream audit and separate score/effect package complete; not a default-summary input |
| GSE186542 | Post hoc external NP score-level support | Human | bulk RNA-seq / processed expression | Early Pfirrmann I-III n = 3 versus advanced IV-V n = 3; six GSM columns map exactly | `GSE186542_gene_expression.txt.gz` | audited; all four locked modules map at 100%; GEO text conflicts between FPKM/TMM and raw-count descriptions, so score-level only |
| GSE167931 | Post hoc external normalized-expression support | Human | bulk RNA-seq of isolated NP cells | Normal n = 4 versus degenerated n = 5; FPKM and TPM representations share the same nine sample columns | `GSE167931_AllSamplesFPKMValue.txt.gz`, `GSE167931_AllSamplesTPMValue.txt.gz` | audited; S8 uses FPKM once, TPM is a same-sample processing sensitivity; patient-level independence is unverified |
| GSE245147 | S9 source-family replacement sensitivity | Human | bulk RNA-seq / RPKM | Native Degenerated n = 3 versus No-degenerated n = 3 comparison | `GSE245147_Degenerated_NO_Degenerated_RPKM.txt.gz` | audited; P2/P8 passage and DMSO/H-151 treatment arms excluded; replaces GSE167931 and is never pooled with it |
| GSE56081 | Candidate conventional microarray | Human | Arraystar Human LncRNA microarray V2.0 | Control grade I n = 5 versus degenerated grade IV/V n = 5 | `GSE56081_series_matrix.txt.gz` | candidate only; a global Ensembl 113/GRCh38 specificity audit fails the 0.80 mapping gate for all four modules, and identical-RNA alternatives are not independent samples |
| GSE160756 | Healthy reference | Human | scRNA-seq | NP n = 3, AF n = 2, CEP n = 2 | `GSE160756_RAW.tar` | not yet queued |
| GSE211407 | Exploratory animal support | Rat | scRNA-seq | Puncture-induced degeneration at 2 and 8 weeks | reported processed files | not primary |

## Interpretation boundaries

GSE229711 and GSE230808 are child studies of the GSE230809 super-series and are a single exploratory discovery project, not independent studies. The shared low-grade group has only three donors per compartment, so this project cannot supply a confirmatory inferential endpoint. GSE160756 is a reference atlas and not a disease-support cohort. GSE211407 is non-human and will not enter the current human direction summary. AF lacks an independent human support cohort in the current design, so AF conclusions are exploratory.

The post hoc accessions are separated from the frozen default result by both
processing scale and identity certainty. GSE186542 is accession/BioProject-level
separate from the default set, but GEO exposes no patient ID, age, sex, or disc
level and its supplied-value description is internally inconsistent. GSE167931
has the same missing patient-level fields; the two linked source articles
(PMIDs 35304463 and 35340126) document GSE167931 data but do not resolve the
GSE186542 provenance conflict. GSE245147 is a source-family replacement for
GSE167931; public metadata cannot establish donor independence from that family.
GSE266883 is excluded because its human design and author chain are highly
similar to GSE245147. These limitations permit exploratory score-level
standardization only and do not establish replication.

## Known metadata risks before download

- The discovery project has complete age-degeneration confounding and limited sex diversity: low-grade donors are male ages 21-27, while advanced-degeneration donors are male ages 37-73.
- GSE165722 reports incomplete demographics in the public record. Its GEO SOFT grade labels (I-IV) conflict with Table 1 of Tu et al. (PMID 34825784), which reports grades II-V. The source-publication grade mapping is used for the protocol (mild II-III versus severe IV-V, n=4 per group), while the conflict is retained in the ledger.
- GSE165722's supplementary matrices are integer-valued but described by GEO as "normalized counts." They can support within-cohort descriptive QC and donor-level score direction checks, but not the project's raw-count negative-binomial pseudobulk model unless raw UMI values are independently recovered.
- GSE244889 has a small mild-versus-severe scRNA contrast (4 versus 3 presumed donor/library keys).  SOFT titles encode ages 17, 55, 24, 30 in MDD and 41, 62, 59 in SDD; this supports neither an age-adjusted nor a confirmatory severity claim.  Its bulk FPKM table is a separate processed asset and must not be merged with scRNA counts or treated as raw counts.
- GSE153066 compares relatively normal trauma-derived NP with degenerated herniation-derived NP. Age and clinical source are confounded; it is independent support, not an unconfounded causal estimate.
- GSE251686 has six GEO records but only five usable matrix payloads after stream audit (mild n=2 versus severe n=3); `GSM7986002` contains NUL bytes and a malformed Matrix Market line and is excluded without repair. It lacks public age/sex/individual-grade metadata, and its platform annotation is internally inconsistent.
- GSE186542 has six count-labelled columns and complete locked-module mapping, but the GEO processing description mentions FPKM/TMM as well as raw counts. It is not eligible for raw-count negative-binomial inference until that conflict is resolved.
- GSE167931 has nine matching normal/degenerated sample columns in FPKM and TPM matrices. Both are normalized-expression representations; S8 predeclares FPKM as the single study representation and treats TPM as a paired processing sensitivity.
- GSE245147 has 18 public BioSamples, but S9 selects only the six native clinical-comparison samples (Degenerated_1-3 versus NO_Degenerated_1-3). Passage P2/P8 and DMSO/H-151 treatment arms are excluded; the matrix is RPKM and the inflammatory module maps 20/21 genes because `TNF` is absent.
- GSE56081 has an exact ten-column sample mapping, but its candidate GPL15314 probes were searched in both orientations against Ensembl release-113 cDNA/ncRNA and the GRCh38 primary assembly. Twenty-eight probes failed the specificity gate, leaving globally specific gene fractions of 0.7083, 0.5714, 0.5000, and 0.6000 across the four modules; all are below 0.80. It remains candidate-only, is not a negative biological result, and does not enter S7-S9.
- Every sample-level field below must be verified from raw archive or GEO metadata before analysis.

## Required sample ledger columns

`dataset`, `accession`, `source_file`, `cell_barcode`, `donor_id`, `library_id`, `batch`, `compartment`, `degeneration_status`, `severity_grade`, `age`, `sex`, `disc_level`, `chemistry`, `original_label`, `analysis_label`, `qc_status`, `exclusion_reason`.

## Current verified discovery ledger

`data/derived/GSE230809_discovery_sample_ledger.csv` is generated directly
from child-series SOFT metadata. Its donor key is `patient_id`, because title
labels such as `donor1` repeat in GSE230808. The ledger has 24 libraries from
13 donors: 11 AF/NP-paired donors (3 healthy and 8 advanced) plus 2 advanced
AF-only donors. The raw archives and all 24 10x triples now pass structural
audit; the separate raw-data ledger is
`data/derived/GSE230809_discovery_raw_data_ledger.csv`.
