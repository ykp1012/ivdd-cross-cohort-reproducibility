# Analysis Protocol

## Scope

This is an observational, public-data, cross-cohort direction-and-heterogeneity study of human intervertebral disc degeneration (IVDD). The primary estimand is the donor-level or presumed donor/sample-key difference in pre-specified molecular-program scores between the recorded lower- and higher-severity groups in nucleus pulposus (NP), described separately by cohort. Annulus fibrosus (AF) and AF-versus-NP analyses in the discovery project are exploratory. The study is not designed to establish cellular lineage, molecular mechanism, clinical prediction, therapeutic efficacy, replication, or an age-independent disease effect.

## Primary question

Across independently processed human NP datasets, do pre-specified IVDD-related transcriptional programs show descriptively aligned or discordant donor-level directions between recorded severity groups, and how are those patterns limited by age, technical platform, processing scale, or cell composition? AF analyses in the GSE230809 discovery project are exploratory and do not have an independent AF support cohort.

## Pre-specified cohorts and roles

| Resource | Role | Unit of inference | Main use |
|---|---|---:|---|
| GSE229711 + GSE230808 (super-series GSE230809) | Exploratory discovery | donor x compartment | Descriptive AF/NP advanced-degeneration effect estimates and paired-compartment sensitivity analyses; not a primary inferential endpoint |
| GSE165722 | External score-level direction support | presumed donor-level Sample/GSM key | NP mild-versus-severe direction and severity-score trend; count model prohibited because GEO describes the supplied matrices as normalized |
| GSE244889 | External directional support | presumed donor/library key | NP mild-versus-severe direction only; 4 versus 3 scRNA groups are not confirmatory, and the separate bulk FPKM asset is not used in this score summary |
| GSE153066 | External count-level support only | presumed donor/library | NP relatively-normal-versus-degenerated direction; retained-cell count aggregation with prior GEO filtering disclosed |
| GSE251686 | Isolated incomplete exploratory NP severity check | presumed sample/library key | Separate audited mild n=2 versus severe n=3 descriptive score/effect package only; deliberately excluded from the default 20-effect summary and never called validation, replication, or decisive support |
| GSE160756 | Reference only | donor/sample | Healthy AF/NP/CEP label and marker reference |
| GSE211407 | Exploratory non-human context | animal | Not part of the current human direction summary |

The healthy and disease subseries under GSE230809 arise from the same parent study and will never be counted as independent discovery and validation cohorts. In that parent study, young low-grade donors and older advanced-degeneration donors are age-confounded. In addition, only three low-grade donors are available per compartment. Its estimates are therefore explicitly exploratory, labelled as advanced-degeneration-associated contrasts, and never interpreted as age-adjusted disease effects, confirmatory inferential results, or an externally replicated finding.

### Post hoc extension analyses

The frozen default 20-effect summary remains the primary descriptive product.
After that summary was locked, three separately packaged NP syntheses were run
as non-confirmatory extensions. They do not alter the default directory,
default effect count, or default sign-alignment display.

| Package | Cohorts and design | Role and boundary |
|---|---|---|
| S7 (`np_exploratory_meta_analysis`) | Four default NP cohorts, one effect per cohort and module (`k = 4`) | Independent exploratory random-effects synthesis; all contrasts remain `confirmatory_eligible=false` |
| S8 (`np_post_hoc_external_expansion_meta_analysis`) | S7 plus GSE186542 (Pfirrmann I-III `n = 3` vs IV-V `n = 3`) and GSE167931 FPKM (normal `n = 4` vs degenerated `n = 5`); `k = 6` | Post hoc score-level expansion; GSE167931 TPM is a paired processing sensitivity, not another cohort; patient-level independence is unverified |
| S9 (`np_source_family_replacement_meta_analysis`) | S8 design with GSE245147 native Degenerated `n = 3` vs No-degenerated `n = 3` replacing GSE167931; `k = 6` | Source-family sensitivity only; GSE167931 and GSE245147 are not pooled, and P2/P8 and DMSO/H-151 arms are excluded |

GSE186542 is used only at the score level because its GEO metadata mention
incompatible FPKM/TMM and raw-count descriptions. GSE167931 is represented by
its FPKM matrix in S8, with the TPM matrix retained only as a same-sample
processing sensitivity. GSE245147 is an RPKM score-level subset; its native
clinical comparison excludes passage and treatment arms. Distinct BioProjects
provide accession-level separation, but the public records do not establish
patient-level independence. GSE56081 remains outside all extensions: an exact
Ensembl release-113 transcriptome and GRCh38 primary-assembly audit found 28
failed candidate probes and fewer than 80% globally specific genes in every
locked module. GSE266883 is excluded because its human design and author chain
are highly similar to GSE245147.

## Primary outcomes

Four pathway families are evaluated as pre-specified, directional program scores:

1. extracellular-matrix organization and collagen remodeling;
2. inflammatory or NF-kB response;
3. hypoxia and oxidative-stress response;
4. disc matrix-homeostasis or matrix-synthesis program.

Before any expression result is examined, each program must receive a locked gene-list source, date, gene-symbol mapping, and hash in the program ledger. Any signature derived from the discovery cohort is exploratory and must be frozen before its first external test.

## Data and QC rules

- Retain raw UMI counts separately from transformed values. A matrix described by GEO as normalized or processed is never treated as raw UMI merely because values are integers.
- Maintain a complete sample ledger: accession, dataset, donor, library, batch, compartment, degeneration status or grade, age, sex, disc level, chemistry, and exclusion reason.
- Choose cell-level QC thresholds from pooled distributions before examining disease-specific results. Initial, non-binding starting checks are at least 200 detected genes and at least 500 UMIs for scRNA-seq; nucleus data are evaluated separately.
- Record mitochondrial percentage, total UMI, detected genes, doublet strategy, retained cells, and exclusion rate per donor.
- Main pseudobulk requires at least 30 high-quality, source-labelled,
  non-excluded compartment cells per donor x compartment. Thresholds of 20 and
  50 cells are sensitivity analyses. These are source-restricted compartment
  pseudobulks, not purified resident-cell populations.
- A group-level difference in QC exclusion rate or library yield triggers investigation and reporting; it is not silently corrected by integration.

## Annotation rules

Anatomical source label is primary. Fixed AF/NP support panels are reported as
concordance evidence but do not overwrite or exclude a source label. Cells with
strong multi-gene immune, endothelial, mural, erythroid, or mixed evidence are
retained in reporting but excluded from source-restricted compartment
pseudobulks. A complete protocol and amendment history are in
`docs/13_annotation_protocol.md`.

## Primary statistical analysis

For each gene within each compartment, raw counts are aggregated over cells:

`PB(g, donor, compartment) = sum(raw UMI count(g, cell))`.

Within AF and NP separately in the discovery project, raw counts may be aggregated to donor-level pseudobulks after QC. Because the low-grade group has only three donors per compartment, the discovery project is restricted to exploratory effect-size estimation, donor-level plots, leave-one-donor-out stability, and direction selection for the already locked modules. It is not used for confirmatory gene-level or module-level significance claims. A donor-level negative-binomial pseudobulk model with TMM library-size normalization and robust quasi-likelihood estimation may be used to obtain descriptive estimates conditional on verified raw UMI matrices, but its p-values are not a primary result. No covariate-adjusted causal interpretation is attempted because age and degeneration status are fully confounded.

Exploratory shared AF/NP patterns require direction agreement and retention of direction in at least 80% of leave-one-donor-out discovery analyses. No shared-effect conjunction p-value is used as a confirmatory claim in this dataset.

Degeneration-by-compartment interaction may be estimated descriptively only when AF and NP are demonstrably paired within donor and design balance permits it. Otherwise, AF and NP are analyzed as separate exploratory contrasts.

## External support and cross-cohort display

Every external dataset is processed independently. The locked gene sets and scoring rule are applied unchanged, and each cohort is displayed separately at the module-score level. GSE165722 is a score-level direction-support cohort because GEO describes its count-like values as normalized; it is not included in raw-count negative-binomial inference. The frozen default comparison is limited to descriptive sign alignment, cohort-specific intervals, and leave-one-key-out stability and has no pooled estimate, p-value, or replication adjudication. S7, S8, and S9 are separate exploratory standardization packages described below; they do not make any cohort confirmatory eligible. AF conclusions remain exploratory unless an independent AF cohort passes audit.

### Exploratory standardized synthesis specification

For S7-S9, one higher-minus-lower effect is supplied for each NP cohort and
module. The primary effect is the heteroscedastic standardized mean difference
(`metafor::escalc(measure = "SMDH")`) and the primary random-effects estimator
is REML with Knapp-Hartung (HKSJ) confidence intervals. Conventional pooled-SD
Hedges *g* and Paule-Mandel tau-squared are model sensitivities. S7-S9 use
`metafor` REML `maxiter = 10000`; all other `control` parameters retain package
defaults. Four module-level HKSJ p-values and their Benjamini-Hochberg (BH)
adjustments are written to the supplementary tables solely for transparent
description. With `k = 4` or `k = 6`, heterogeneous processing scales, and
unverified patient-level independence, these p-values, heterogeneity
statistics, prediction intervals, and leave-one-cohort-out results cannot be
interpreted as tests of confirmation, replication, biomarkers, mechanisms,
causal effects, or treatment response.

The primary SMDH estimates are recorded here to keep the analysis layers
traceable. S7 gives ECM 1.0481 (95% CI -0.9965 to 3.0926), inflammatory 0.2809
(-1.1432 to 1.7050), hypoxia 0.8184 (-0.1285 to 1.7654), and matrix homeostasis
0.1453 (-0.3679 to 0.6585). S8 gives ECM 0.7780 (-0.3532 to 1.9093),
inflammatory 0.4032 (-0.3991 to 1.2055), hypoxia 0.7694 (0.1706 to 1.3682),
and matrix homeostasis 0.3762 (-0.3017 to 1.0540); the HKSJ p-value for hypoxia
is 0.0214 and its four-module BH value is 0.0856. S9 gives ECM 1.0931
(-0.9186 to 3.1048), inflammatory 0.4056 (-1.0004 to 1.8115), hypoxia 0.5746
(-0.7231 to 1.8723), and matrix homeostasis 0.1046 (-1.0200 to 1.2292).
The corresponding S9 BH values are 0.6133, 0.6556, 0.6133, and 0.8205,
respectively; all four exceed 0.61. The S9 replacement therefore does not
provide a stable confirmatory signal, and the S8 hypoxia interval that was
above zero is not robust to this source-family replacement.

GSE165722 has now passed exact CellIndex-to-matrix-header and module mapping
audits. Its 8 presumed Sample/GSM keys (mild 4, severe 4) contribute a
score-level severe-minus-mild direction and leave-one-key-out display only.
GEO does not expose an independent patient identifier or demographic covariates,
and the supplied values remain normalized-count inputs rather than raw UMI
counts. It is not confirmatory and cannot support a raw-count model.

GSE153066 provides an additional independent NP normal-versus-degenerated
support contrast. Its compact count matrix, 16 sample-prefixed barcode groups,
and exact barcode-prefix-to-GEO-sample mapping have passed audit. Its age and
clinical source are confounded and it is not interpreted as an unconfounded
disease effect. GSE251686 now has a separately generated, audit-gated score
table and isolated descriptive mild n=2 versus severe n=3 effect package. It
is deliberately excluded from the default 20-effect cross-cohort summary:
`GSM7986002` is permanently excluded for malformed Matrix Market text and
must not be repaired or used, while the five remaining records have only
presumed sample/library identities and unavailable covariates. The separate
package is not validation, replication, or decisive external support.

GSE244889 has seven audited raw 10x NP libraries, but its MDD-versus-SDD
contrast has only four versus three title-derived presumed donor/library keys.
It is therefore restricted to directional support, effect-size, and
leave-one-donor-out stability displays.  Its six-sample FPKM table is a
separate processed bulk asset: it is not merged with scRNA counts or used for
raw-count inference.

## Sensitivity analyses

- pseudobulk cell thresholds of 20, 30, and 50;
- source labels versus marker-based versus reference-based annotation;
- leave-one-donor-out analysis and donor bootstrap;
- per-donor cell downsampling only if it can be defined before outcome review;
- no adjusted model where covariates are structurally confounded or unavailable;
- exclusion of low-quality or influential donors with reasons reported;
    - separate cohort effects with descriptive sign alignment; S7-S9 standardized syntheses remain isolated, exploratory, and non-confirmatory;
- ordinal and binary severity analyses where grades permit;
- gene mapping retention for every score, requiring at least 80% measured genes for comparable module scoring.

## Stop and downgrade rules

| Gate | Action if failed |
|---|---|
| Cells cannot be uniquely linked to donor or library | Stop the donor-level analysis for that cohort |
| A compartment has fewer than four disease and four comparison donors after QC | Do not make inferential claims for that compartment; retain only exploratory effect-size and stability displays if the cohort still passes all QC and provenance gates |
| Degeneration is fully confounded with batch, location, sex, or age | Do not attribute the effect to degeneration independently; restrict the cohort to descriptive advanced-degeneration association |
| No verified raw or count-like matrix exists | Do not run the count-model primary analysis; use only descriptive score-level support |
| AF/NP cannot be reliably separated | Do not make compartment-specific claims |
| Discovery direction stability is below 80% | Do not emphasize the module as a stable discovery pattern; retain it only in the complete exploratory display |
| External cohort directions are discordant | Do not call the program consistently aligned; report heterogeneity without pooling scales |
| Only one cohort supports a degeneration-by-compartment interaction | Do not use interaction language in title or abstract |

## Reporting

All results will report donor counts, cell counts, effect sizes, confidence intervals where estimable, exclusions, mapping losses, and all pre-specified results. The frozen default descriptive analysis reports no p-values or formal tests. S7-S9 transparently report HKSJ and four-module BH p-values, but these post hoc/non-confirmatory values cannot be used for confirmation, replication adjudication, biomarker, mechanistic, causal, or therapeutic claims. Any future confirmatory hypothesis-testing analysis must pre-specify its multiplicity family, adjustment method, cohort-independence criteria, and estimand. Conclusions will use association and exploratory language only.
