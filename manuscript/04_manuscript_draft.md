# Manuscript Draft

## Title
Cohort-aware direction and heterogeneity of locked nucleus pulposus transcriptional programs across human intervertebral disc degeneration datasets

## Abstract

Public single-cell studies of intervertebral disc degeneration (IVDD) often
contain many cells but few independently observed donors, and differences in
matrix processing and clinical sampling complicate cross-cohort comparison. We
used public functional-genomics records from the NCBI Gene Expression Omnibus
(GEO) [1] and treated cells as nested observations, consistent with warnings
that ignoring biological replicates can produce false discoveries in
single-cell analyses [2]. We performed an evidence-limited, cohort-specific
analysis of four project-locked nucleus pulposus (NP) transcriptional programs:
extracellular-matrix and collagen remodeling, inflammatory/NF-kB response,
hypoxia/oxidative stress, and disc-matrix homeostasis. The experimental unit
was the donor or, where a public identifier was incomplete, the explicitly
labelled presumed donor/sample-library key; cells were nested observations and
were not used as independent replicates. We analyzed the GSE230809 parent
project (GSE229711 plus GSE230808) as one exploratory AF/NP discovery project
[3-6] and displayed NP scores separately in GSE244889 [7,8], GSE153066 [9],
and GSE165722 [10,11]. Raw 10x matrices were aggregated only after archive,
QC, annotation, and sample-identity audits. GSE153066 used its contributor-
retained dense count matrix, and GSE165722 was restricted to normalized-count
score-level directional analysis. The frozen default summary reports cohort-
specific unweighted group differences, 95% Welch and donor/bootstrap intervals,
and leave-one-key-out stability; it contains no p-values, formal meta-analysis,
or replication adjudication. Separately, S7, S8, and S9 are non-confirmatory
SMDH random-effects syntheses with REML and Knapp-Hartung intervals; their HKSJ
and BH p-values are reported only for transparency.

Across the four cohorts in the default NP summary, the hypoxia/oxidative-stress
program had positive higher-recorded-severity point estimates, but all four
corresponding Welch 95% intervals included zero. ECM, inflammatory/NF-kB, and
disc-matrix-homeostasis directions were discordant between cohorts. The S8
post hoc six-cohort expansion yielded a hypoxia SMDH of 0.7694 (95% CI 0.1706
to 1.3682), but the S9 source-family replacement of GSE167931 with GSE245147
gave 0.5746 (-0.7231 to 1.8723). A separately audited GSE251686 analysis
remained isolated because one record failed stream-integrity audit and the
usable comparison was non-balanced. This study provides an auditable,
cohort-aware account of where program-level directions align descriptively and
where they remain heterogeneous; it does not establish a universal IVDD
program, mechanism, biomarker, or therapeutic target.

## Introduction

Intervertebral disc degeneration is a clinically important and biologically
heterogeneous condition. Public single-cell transcriptomic datasets provide an
opportunity to compare molecular programs across human samples [3,7,10,12,13], but
cell counts do not replace independent donor- or sample-key-level observations
[2,14].
Analyses that pool cells across cohorts can make a direction appear precise
while obscuring donor nesting, clinical-source differences, and processing
boundaries.

Recent human IVDD studies have also integrated multiple public datasets or
focused on specific NP states [12,13,15]. We therefore do not claim that the
present study is the first multi-cohort IVDD integration; its narrower aim is
to preserve the available observation structure and report where cohort-specific
NP program directions agree or disagree.

Here, we asked whether four project-locked NP expression programs showed
aligned or discordant directions across independently processed public
datasets. The study was designed as a provenance-aware descriptive evidence
audit. It was not designed to identify an age-independent disease effect,
establish a cellular lineage or mechanism, develop a clinical classifier, or
nominate a treatment target.

## Methods

### Cohorts and inference boundary

Table 1 and Supplementary Table S2 describe cohort roles, recorded group
structures, observation keys, and identity checks.

GSE229711 and GSE230808 were treated as the single GSE230809 parent project,
not as discovery and validation cohorts. This project contains three young
low-grade donors per compartment and older advanced-degeneration donors, so
age and recorded disease state are fully confounded. AF and NP contrasts from
this parent project were therefore exploratory advanced-degeneration-associated
displays only.

External NP support cohorts were analyzed independently. GSE244889 supplied
four mild and three severe title-derived presumed donor/library keys. GSE153066
supplied eight relatively normal and eight degenerated sample-prefixed keys,
but clinical source and age were confounded with disease status. GSE165722
supplied four mild and four severe presumed sample keys using the source
publication's severity grouping; GEO describes its integer-like matrices as
normalized counts, so it was restricted to score-level direction and stability
analysis [10,11]. GSE251686 underwent a separate audited 2-versus-3 exploratory
score analysis [16,17]. During this project's stream-integrity audit,
`GSM7986002` was permanently excluded. This non-balanced display was
deliberately omitted from the default 20-effect summary because of the
integrity failure and limitations in sample identity and covariate information;
it was neither used as validation nor treated as replication. Supplementary
Figure S2 depicts this cohort disposition and the separation of the default and
isolated analysis streams.

The cohort provenance was checked against the corresponding GEO series and,
where available, the linked source publication: GSE244889 [7,8], GSE153066
[9], and GSE251686 [16,17]. The GSE153066 GEO record currently lists no
associated citation, so that accession is cited as a database record rather
than attributed to an unverified article [9].

After the default analysis was frozen, GSE186542 was audited as a small
early-Pfirrmann-I--III versus advanced-IV--V NP score-level contrast (3 versus
3) [22]. Its GEO SOFT metadata do not list a linked PubMed record, so it is
cited as an accession only. GSE167931 supplied one normal-versus-degenerated
FPKM representation (4 versus 5); its paired TPM matrix was retained only as a
same-sample processing sensitivity [23-25]. GSE245147 supplied a native
Degenerated-versus-No-degenerated RPKM subset (3 versus 3) after passage and
treatment arms were excluded [26,27]. GSE167931 and GSE245147 were not pooled
together because their source family and patient-level independence could not
be resolved from public metadata.

### Locked module scoring and effect displays

The four module definitions were locked within the project before external
scoring. This was a timestamped, versioned analysis lock rather than a
prospective study registration. Supplementary Table S3 records the locking
time, source identifiers, gene lists, and SHA-256 hashes. At least 80% of the
locked genes had to map in a score. The pathway components were anchored to the
KEGG and Reactome knowledgebases [18,19] and to the cited IVDD source studies
[3,10]. A higher score means higher expression of the listed genes only; it
does not imply benefit, harm, causality, or therapeutic relevance.
For verified raw 10x matrices, source-restricted,
QC-passing cells were aggregated within donor or library and each module score
was the mean mapped-gene `log1p(CPM)` value. For the externally supplied dense
or normalized matrices, the same score definition was applied only within the
declared matrix-processing boundary. Original archives, input hashes, and
score-to-ledger identities were retained in the reproducibility contract
(Supplementary Table S6).

For every cohort, compartment, and module, we calculated the unweighted mean
difference between the recorded higher- and lower-severity groups. Welch and
within-group donor/library bootstrap intervals were displayed as descriptive
uncertainty summaries. Leave-one-key-out analyses assessed whether the sign
was retained; the complete default leave-one-key-out output is provided in
Supplementary Table S4. The frozen default cross-cohort sign-alignment display
was summarized without a pooled effect, formal meta-analysis, p-value, or
replication decision.

### Separate exploratory standardized syntheses

S7-S9 used one higher-minus-lower NP effect per cohort and module. The primary
effect was the heteroscedastic standardized mean difference (SMDH), fitted with
a REML random-effects model and Knapp-Hartung confidence intervals. Conventional
pooled-SD Hedges *g* and Paule-Mandel tau-squared estimates were model
sensitivities. All three packages recorded `metafor` REML `maxiter = 10000`,
while retaining all other `metafor` control defaults. Four module-level HKSJ
p-values and their Benjamini-Hochberg (BH) adjustments were included solely for
transparent reporting, not as confirmation or replication tests.

S7 independently standardized the four default NP cohorts (`k = 4`). S8 was a
post hoc six-cohort score-level expansion that added GSE186542 and GSE167931
FPKM (`k = 6`). S9 used the same design but replaced GSE167931 with the native
GSE245147 subset (`k = 6`). These packages did not modify the frozen default
20-effect summary, and small cohort counts, cross-platform score scales,
unverified patient-level independence, and source-family uncertainty precluded
confirmatory interpretation.

The donor-level treatment of cells follows the multi-sample principle used in
single-cell population-level analyses [14]. Welch intervals were used for the
unweighted difference because the two arms need not have equal variances [20],
and percentile intervals were generated by resampling donor/library keys
within each arm [21].

## Results

The default descriptive summary contains 20 cohort/compartment/module effects
from the four score cohorts listed in Table 1, with 55/55 exact
score-to-ledger sample-key matches (Supplementary Tables S2 and S6). All four
cohorts have `confirmatory_eligible=false`. GSE230809 contributes exploratory
AF and NP contrasts; GSE244889 contributes a 4-versus-3 NP directional-support
contrast; GSE153066 contributes an 8-versus-8 NP count-level-support contrast;
and GSE165722 contributes a 4-versus-4 NP normalized-count score-level
direction contrast. GSE251686 is not part of this default summary, sign
alignment, or either main figure (Supplementary Figure S2).

In NP, the hypoxia/oxidative-stress module had a positive
higher-recorded-severity point estimate in each of the four default contrasts:
0.1776 in GSE230809, 0.2381 in GSE244889, 0.0882 in GSE153066, and 0.4346 in
GSE165722 (Table 2 and Figures 1-2). The corresponding descriptive Welch 95%
intervals were `[-0.0245, 0.3798]`, `[-0.1314, 0.6076]`,
`[-0.2166, 0.3929]`, and `[-0.0619, 0.9310]`; every interval included zero.
Thus, the four positive point estimates are a descriptive sign pattern, not a
robust cross-cohort replication result. ECM/collagen remodeling was positive
in GSE230809, GSE244889, and GSE165722 but negative in GSE153066.
Inflammatory/NF-kB was positive in GSE244889 and GSE153066 but negative in
GSE230809 and GSE165722. Disc-matrix homeostasis was positive in GSE230809,
GSE244889, and GSE165722 but negative in GSE153066 (Table 2 and Figures 1-2).

The separately audited GSE251686 sensitivity analysis retained five presumed
sample/library keys after permanent exclusion of `GSM7986002` (mild n=2 and
severe n=3). Its hypoxia/oxidative-stress point estimate was negative
(`-0.0778`; Welch 95% interval `[-1.7214, 1.5658]`), and the four module
effects are reported in Supplementary Table S1 and Supplementary Figure S1.
This isolated result was not pooled, counted in default sign alignment, or used
to alter the four-cohort default result.

The separate S7 synthesis yielded SMDH values of 1.0481 (95% CI -0.9965 to
3.0926) for ECM/collagen remodeling, 0.2809 (-1.1432 to 1.7050) for
inflammatory/NF-kB, 0.8184 (-0.1285 to 1.7654) for hypoxia/oxidative stress,
and 0.1453 (-0.3679 to 0.6585) for disc-matrix homeostasis (Supplementary
Tables S7a-S7d and Supplementary Figure S3). These four exploratory intervals
all included zero.

The post hoc S8 expansion yielded SMDH values of 0.7780 (-0.3532 to 1.9093),
0.4032 (-0.3991 to 1.2055), 0.7694 (0.1706 to 1.3682), and 0.3762 (-0.3017 to
1.0540) for ECM, inflammatory, hypoxia, and homeostasis, respectively
(Supplementary Tables S8a-S8d and Supplementary Figure S4). The hypoxia HKSJ
p-value was 0.0214 and its four-module BH value was 0.0856; these values are
reported for transparency only. In S9, which replaced GSE167931 with GSE245147,
the corresponding SMDH values were 1.0931 (-0.9186 to 3.1048), 0.4056
(-1.0004 to 1.8115), 0.5746 (-0.7231 to 1.8723), and 0.1046 (-1.0200 to
1.2292) (Supplementary Tables S9a-S9d and Supplementary Figure S5). The four
S9 BH values were 0.6133, 0.6556, 0.6133, and 0.8205. The loss of the
non-zero-crossing S8 hypoxia interval after replacement is a source-family
sensitivity finding, not confirmation.

AF results were available only from the GSE230809 parent project. Because the
low-grade AF group contained three donors and was age-confounded with the
advanced-degeneration group, these estimates are retained as exploratory
context rather than external evidence; complete default leave-one-key-out
results are available in Supplementary Table S4.

The locked discovery retained-cell eligibility sensitivity analysis did
not change the GSE230809 input set. All 24 donor/libraries passed the 20-,
30-, and 50-cell thresholds; the minimum source-restricted, QC-passing count
was 471 cells. The 96 library-module scores and all eight AF/NP descriptive
effects were identical across these thresholds (Supplementary Tables S5a and
S5b). This result only shows that no library lay near the chosen eligibility
cutoffs; it does not test robustness to alternative annotation, cell
composition, or random cell downsampling.

## Discussion

This cohort-aware analysis separates a modest, audit-supported observation
from stronger claims the current data cannot sustain. The positive NP
hypoxia/oxidative-stress point estimates in the four default contrasts are a
descriptive cross-cohort pattern, but every associated Welch interval includes
zero. They therefore do not demonstrate replication, mechanism, a severity
biomarker, or a therapeutic target. The discordance of the ECM, inflammatory,
and matrix-homeostasis scores shows why pooled cells or a single cohort should
not be used to imply a universal degeneration program.

The S7-S9 syntheses add a standardized effect-size description but do not
remove this limitation. In particular, the S8 hypoxia interval excluded zero,
whereas the S9 source-family replacement interval included zero. That change,
along with post hoc cohort selection, `k = 4` or `k = 6`, processing-scale
differences, and unknown patient overlap, means that neither its p-values nor
its intervals can establish replication, a biomarker, mechanism, or treatment
effect.

The study is limited by small cohorts, presumed rather than fully verified
sample identities in several resources, complete age-disease confounding in
GSE230809, clinical-source and age confounding in GSE153066, normalized-count
restrictions in GSE165722, and the lack of an independent AF severity cohort.
GSE251686 was transparently reported as an isolated sensitivity analysis rather
than selected into the default summary after scoring. The four scores are locked
expression summaries rather than clinical indices, and their cross-platform
magnitudes are not pooled. Future work would require independently collected,
adequately powered human cohorts with verified sample nesting and orthogonal
biological measurements before any causal, mechanistic, prognostic, or
therapeutic conclusion could be made.

## Conclusions

Across four cohorts included in the default NP summary, only the direction of
the hypoxia/oxidative-stress point estimate aligned, and its four descriptive
Welch intervals all included zero. The remaining module directions were
discordant. The separate S7-S9 syntheses were also non-confirmatory and
source-family-sensitive. These results support transparent cohort-specific
reporting and do not establish a universal IVDD transcriptional program.

## Table and Figure Legends

**Table 1. Cohort roles and inference boundaries in the default descriptive
summary.** The table identifies the recorded group structure, observation key,
score-filtered keys, exact score-to-ledger identity check, and interpretation
boundary for each cohort contributing to the default 20-effect summary.
GSE229711 and GSE230808 are shown as the single GSE230809 parent project.
GSE251686 is intentionally excluded; its cohort disposition is documented in
Supplementary Table S2 and Supplementary Figure S2, while its isolated effects
are reported in Supplementary Table S1 and Supplementary Figure S1.

**Table 2. Cohort-specific NP module effects in the default descriptive
summary.** Each row reports the unweighted higher-minus-lower recorded-severity
module-score difference, descriptive Welch and bootstrap 95% intervals, and
leave-one-key-out direction-retention fraction. Values remain cohort specific
and are not pooled in this default layer. No p-values, formal meta-analysis, or
replication decision is reported for Table 2.

**Figure 1. Cohort-specific NP module-score differences and descriptive Welch
95% intervals.** Points are unweighted donor or presumed sample-key mean
differences between recorded higher- and lower-severity groups. Bars are
descriptive Welch 95% intervals. The four panels show the locked module scores;
the horizontal axis is unitless and cohort specific, and GSE165722 uses supplied
normalized-count values rather than raw CPM. No pooled estimate or replication
decision is shown.

**Figure 2. NP cohort-specific directions and descriptive sign alignment.**
Cells encode the sign of the cohort-specific higher-minus-lower module-score
difference, with the numerical difference printed in each cell. The alignment
display is descriptive only; it does not pool magnitudes, test a hypothesis, or
establish replication.

**Graphical abstract. Cohort-aware IVDD program audit.** Public human IVDD
cohorts were scored using donor or presumed sample/library keys as the
observation unit, with cells retained as nested observations. Program definitions
were locked before external scoring. Four default NP contrasts showed positive
hypoxia/oxidative-stress point estimates, although all corresponding Welch
intervals included zero; the other program directions were discordant. The
isolated GSE251686 sensitivity analysis did not enter the default summary.

**Supplementary Table S1. Isolated GSE251686 exploratory sensitivity effects.**
The table reports all four module effects for the mild n=2 versus severe n=3
comparison after permanent exclusion of `GSM7986002`. It is excluded from the
default summary, sign alignment, and main figures and is not validation or
replication evidence.

**Supplementary Table S2. Cohort disposition and identity audit.** The table
lists all default and isolated analysis streams, observation keys, group
structures, identity or eligibility checks, default-summary inclusion, and
interpretation boundaries.

**Supplementary Table S3. Locked program module definitions.** The table
provides the source classes, identifiers, complete locked gene lists,
timestamp, score-direction convention, and SHA-256 hash for each module.

**Supplementary Table S4. Complete default leave-one-key-out results.** The
table records each leave-one-key-out effect calculation for the 20-effect
default summary, including AF and NP displays from the GSE230809 parent project.

**Supplementary Table S5a. Discovery retained-cell threshold summary.** The
table records the 20-, 30-, and 50-cell eligibility runs, the number of passing
libraries, and the minimum observed source-restricted QC-passing cell count.

**Supplementary Table S5b. Discovery retained-cell threshold effect stability.**
The table compares all eight GSE230809 AF/NP descriptive effects at each
threshold with the 30-cell reference run.

**Supplementary Table S6. Reproducibility contract.** The table indexes the
default-summary manifest, identity crosswalk, locked program ledger,
threshold-sensitivity inputs, isolated GSE251686 manifest, S7-S9 manifests and
generators, artifact generators, and project-local Python environment lock by
SHA-256 hash.

**Supplementary Tables S7a-S7d. Exploratory four-cohort NP random-effects
synthesis.** The tables give S7 study-level effects, primary SMDH/REML/
Knapp-Hartung results, Hedges *g* and Paule-Mandel sensitivities, and
leave-one-cohort-out results. They are a separate non-confirmatory synthesis of
the four default NP cohorts and do not replace the default descriptive analysis.

**Supplementary Tables S8a-S8d. Post hoc six-cohort NP expansion.** The tables
give the separately packaged S8 study-level effects, primary results, model
sensitivities, and leave-one-cohort-out results after adding GSE186542 and the
GSE167931 FPKM representation. HKSJ and BH p-values are transparent descriptive
outputs only and not patient-level validation or confirmation.

**Supplementary Tables S9a-S9d. Source-family replacement sensitivity.** The
tables give the S9 results after substituting the native GSE245147 comparison
for GSE167931. The two source-family cohorts are never pooled together; S9 is a
non-confirmatory sensitivity analysis, not a seventh independent cohort.

**Supplementary Figure S1. Isolated GSE251686 exploratory sensitivity display.**
The figure presents cohort-specific effect estimates, descriptive Welch 95%
intervals, and leave-one-key-out direction retention for GSE251686. The display
remains isolated from the default summary.

**Supplementary Figure S2. Cohort disposition and analysis boundary.** The
flow diagram distinguishes the four default score cohorts from the separately
scored GSE251686 sensitivity package and states the inferential limits of the
analysis.

**Supplementary Figure S3. Exploratory four-cohort NP random-effects
synthesis.** The forest display gives S7 SMDH estimates and Knapp-Hartung
intervals. It is a separate exploratory quantitative summary and not a
confirmation or replication display.

**Supplementary Figure S4. Post hoc six-cohort NP expansion.** The forest
display gives S8 results after adding GSE186542 and GSE167931 FPKM. The plotted
SMDH estimates, intervals, and transparent p-values are non-confirmatory.

**Supplementary Figure S5. NP source-family replacement sensitivity.** The
forest display gives S9 results after replacing GSE167931 with the native
GSE245147 subset. It demonstrates source-family sensitivity and does not add
an independent validation cohort.

## Data Availability

All source data analyzed in this study are publicly available through NCBI GEO
under GSE230809, GSE229711, GSE230808, GSE244889, GSE153066, GSE165722, and
GSE251686, GSE186542, GSE167931, and GSE245147. GSE56081 was audited as a
candidate-only microarray extension and was not analyzed for module effects.
The project code, configurations, environment lock, derived tables, and
manifests have not yet been archived in a public versioned repository.
Before submission, the exact release used for the manuscript must be deposited
in a public repository and an archival service that provides a persistent DOI;
the repository URL and DOI should be inserted here: `[repository URL and DOI to
be added before submission]`.

## Code Availability

The versioned public code release must include the analysis scripts, locked
configuration files, environment lock, generated result tables, and manifests
referenced in Supplementary Table S6. This manuscript does not claim that code
has already been publicly archived.

## Ethics Statement

This work reanalyzes de-identified, publicly available data and involves no new
participant recruitment, intervention, or specimen collection. Before
submission, the corresponding author must confirm the applicable local
institutional determination and add any required approval, exemption, or waiver
identifier: `[ethics determination to be completed before submission]`.

## Funding

`[Funding sources, grant numbers, and the funders' role to be completed before
submission.]`

## Competing Interests

`[Competing-interest declaration to be completed by all authors before
submission.]`

## Author Contributions

`[Author names and CRediT contributions to be completed and approved by all
authors before submission.]`

## References
See `references.bib`.
