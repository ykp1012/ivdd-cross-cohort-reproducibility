# Data Audit

Audit date: 2026-08-13. Updated: 2026-08-14.

This audit began with GEO SOFT metadata, supplement file lists, and linked publication metadata. GSE165722, both GSE230809 child-series archives, and the post hoc candidate accessions have now been downloaded and format-audited. Raw archives are never altered.

| Dataset | Biological samples and parseable metadata | Files available | Planned role | Major constraint |
|---|---|---|---|---|
| GSE229711 | Human low-grade Thompson II; three male donors, ages 21, 25, and 27; paired AF and NP for each donor | Six audited 10X barcode, feature, and MTX triplets; raw archive 464,261,120 bytes | Low-grade portion of the exploratory discovery project | Three donors only; all young male; below confirmatory inference threshold |
| GSE230808 | Human advanced degeneration; ten male donors ages 37-73; AF from ten donors and NP from eight, with eight paired AF/NP donors | Eighteen audited 10X barcode, feature, and MTX triplets; raw archive 737,361,920 bytes | Advanced-degeneration portion of the same exploratory discovery project | Age and degeneration status are completely confounded; not independent of GSE229711 |
| GSE165722 | Human NP only; eight presumed donor-level Sample/GSM keys. GEO labels grades I-IV; the source article reports II-V, two keys per grade | Eight BD Rhapsody count-like TSVs plus cell-name files; raw archive about 67 MiB; exact CellIndex-to-matrix-header audit and locked score table complete | Independent NP score-level severity direction support | No AF/incomplete public age-sex metadata; GEO calls files normalized counts despite integer values, so they are not raw-count-model input or confirmatory evidence |
| GSE153066 | Human NP; 16 samples, relatively normal n = 8 and degenerated n = 8 | Audited 100,366,013-byte combined count matrix with 36,847 genes and 85,205 unique sample-prefixed cell barcodes; each of 16 prefixes exactly matches one GEO sample title | Independent NP normal-versus-degenerated support | Age and clinical source are confounded; sample title is a presumed donor/library key and must be disclosed |
| GSE244889 | Human NP; scRNA mild n = 4 and severe n = 3, plus a separate six-sample mild/severe bulk FPKM table | Seven audited 10X triples, all 37,487 features; 546,119,680-byte raw archive; six-column processed FPKM table | Small external NP directional reproducibility and separately normalized bulk direction check | 4 versus 3 presumed donor/library keys is below confirmatory threshold; title-encoded age/grade structure and bulk processed-value boundary must be reported |
| GSE251686 | Human NP; six GEO records, of which mild n = 2 and severe n = 3 pass stream integrity | Six per-GSM nested matrix archives, about 275 MiB total; `GSM7986002` fails text-integrity audit | Incomplete, non-balanced NP severity direction check only | `GSM7986002` has 5,471,800 NUL bytes and one malformed line; exclude without repair. No public age/sex or individual grade; platform metadata conflict remains documented |
| GSE186542 | Human NP bulk expression; early Pfirrmann I-III n = 3 versus advanced IV-V n = 3 | Six count-labelled matrix columns; exact GSM mapping; all four locked modules at 100% coverage | S8/S9 post hoc score-level external support | GEO processing text conflicts between FPKM/TMM and raw-count descriptions; patient ID, age, sex, and disc level are not exposed; not raw-count-model eligible |
| GSE167931 | Human isolated NP bulk expression; normal n = 4 versus degenerated n = 5 | Matching nine-column FPKM and TPM matrices; all four modules at 100% coverage | S8 post hoc score-level external support; FPKM used once | TPM is a same-sample processing sensitivity, not an extra cohort; patient-level independence and publication-level overlap are unverified |
| GSE245147 | Human NP bulk expression; native Degenerated n = 3 versus No-degenerated n = 3 | Audited 12-column RPKM matrix; six native comparison columns selected; all selected values finite; 20/21 inflammatory genes mapped | S9 source-family replacement for GSE167931 | P2/P8 passage and DMSO/H-151 treatment arms excluded; source-family overlap cannot be ruled out; score-level only |
| GSE56081 | Human NP microarray; control grade I n = 5 versus degenerated grade IV/V n = 5 | Ten-column log2 quantile-normalized matrix with exact GSM mapping | Candidate only; excluded from S7-S9 | Ensembl-113 transcriptome and GRCh38 primary-assembly specificity audit fails the 0.80 gate in every module; identical-RNA alternatives are not independent samples |
| GSE160756 | Human healthy atlas; NP n = 3, AF n = 2, CEP n = 2 | Seven loom archives; raw archive about 433 MiB | Reference for AF/NP label and marker checks | Not a disease-validation dataset; donor nesting needs careful confirmation |
| GSE211407 | Rat puncture model; eight animals at two or eight weeks | Sixteen processed, cluster-sorted count TXT files; raw archive about 432 MiB | Exploratory non-human context only | Not a standard cell-by-gene raw matrix and cannot support a fresh full single-cell workflow |

## Consequences for the manuscript claim

The frozen default comparison is an exploratory advanced-degeneration-associated
contrast, not a degeneration-only effect, because age is fully confounded and
the comparison group has only three donors. Its cross-cohort endpoint is
descriptive direction alignment or discordance of pre-specified NP program
scores across independently processed cohorts. GSE165722 contributes score-
level direction evidence only; raw-count pseudobulk inference requires
separately verified raw UMI matrices. S7 is an independent four-cohort
standardized synthesis, whereas S8 and S9 are post hoc six-cohort extensions;
all three are non-confirmatory and cannot establish patient-level replication.
The AF/NP shared-versus-compartment-specific findings are discovery-stage and
exploratory until an independent AF cohort is identified.

## Verified discovery metadata before matrix ingestion

GSE229711 and GSE230808 together yield 24 libraries from 13 biological donors:
three healthy and ten advanced-degeneration donors, all male. The healthy group
contains three paired AF/NP donors; the advanced group contains eight paired
AF/NP donors and two AF-only donors. Therefore AF has 13 donor-compartment
libraries and NP has 11. `patient id`, not the title's repeating `donor1`,
is the canonical donor key. Ages are 21, 25, and 27 in the healthy group and
37-73 in the advanced group, making age and disease state fully confounded.
Disc level is also imbalanced and will not be fitted as an adjusted covariate
at this sample size. The machine-readable ledger is
`data/derived/GSE230809_discovery_sample_ledger.csv`.

## Current disposition after audit

GSE165722, GSE229711, GSE230808, GSE153066, and GSE244889 contribute to the
frozen default descriptive package. GSE251686 remains an isolated, non-balanced
exploratory display. GSE186542 and GSE167931 contribute only to the separately
packaged S8 post hoc synthesis, with GSE167931 represented once by FPKM.
GSE245147 contributes only to S9 as a source-family replacement for GSE167931.
GSE56081 remains blocked at the candidate stage: its global probe-specificity
audit found 28 failed candidate probes, and no module reached 80% globally
specific-gene coverage. GSE266883 is excluded because of its similar human
design and author chain to GSE245147. GSE160756 remains a healthy reference
and GSE211407 remains non-human context; neither enters S7-S9.
