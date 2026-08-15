# Cohort-Aware Descriptive Analysis of Directionality and Heterogeneity of Nucleus Pulposus Transcriptional Programs across Public Human Intervertebral Disc Degeneration Datasets

## Title Page

**Article type:** Original Research

**Short title:** Cohort-aware IVDD program analysis

**Authors:** [Author names to be inserted after author approval]

**Affiliations:** [Affiliations to be inserted after author approval]

**Corresponding author:** [Name, postal address, email, and telephone number to be inserted]

**Word count:** approximately 2,456 words, excluding title page, references, figure legends, and tables

**Figures and tables:** graphical abstract; 2 main figures; 2 main tables; 5 supplementary figures; 9 supplementary table packages

## Abstract

Public single-cell intervertebral disc degeneration (IVDD) datasets often contain many cells but few independently observed donors, and clinical and processing differences complicate cross-cohort comparison. We conducted a cohort-aware descriptive analysis of four nucleus pulposus (NP) transcriptional programs: extracellular matrix (ECM)/collagen remodeling, inflammatory/nuclear factor-kappa B (NF-κB) response, hypoxia/oxidative stress, and disc-matrix homeostasis. The observation unit was the donor or, where needed, an explicitly labelled presumed donor, sample, or library key; cells were nested observations. Four audited NP cohorts contributed to the frozen default summary, which contained 20 effects: 16 NP effects and 4 exploratory annulus fibrosus (AF) effects from the GSE230809 parent project. We estimated cohort-specific differences between groups recorded as higher versus lower severity with descriptive Welch and donor or library bootstrap 95% intervals and leave-one-key-out stability. Exactly 55 of 55 score-to-ledger identity matches passed. Hypoxia/oxidative stress was the only NP program with positive point estimates in all four default cohorts (0.1776, 0.2381, 0.0882, and 0.4346), although every Welch interval included zero. The other three programs showed discordant directions. Separate exploratory standardized syntheses did not alter the default analysis. In a post hoc six-cohort expansion, the hypoxia standardized mean difference was 0.7694 (95% CI 0.1706 to 1.3682); after source-family replacement, it was 0.5746 (95% CI -0.7231 to 1.8723). These results suggest a descriptive hypoxia-related direction for future study amid substantial cohort-level heterogeneity. They do not establish a universal IVDD program, biological mechanism, biomarker, or therapeutic target.

## Keywords

intervertebral disc degeneration; nucleus pulposus; public transcriptomic data; single-cell transcriptomics; cohort heterogeneity; reproducibility

## Graphical Abstract

Public human IVDD cohorts were evaluated with donor or explicitly labelled presumed sample or library keys as the observation unit. Cells remained nested observations. Program definitions were locked before external scoring.

[[FIGURE:graphical_abstract]]

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

### Study design, public datasets, and inference boundary

Table 1 and Supplementary Table S2 describe cohort roles, recorded group
structures, observation keys, and identity checks. Annulus fibrosus (AF) denotes
the outer annular compartment; extracellular matrix (ECM) denotes the matrix
program label; and nuclear factor-kappa B (NF-κB) denotes the inflammatory
transcriptional response label.

GSE229711 and GSE230808 were treated as the single GSE230809 parent project,
not as discovery and validation cohorts. This project contains three young
low-grade donors per compartment and older advanced-degeneration donors, so
age and recorded disease state are fully confounded. AF and NP contrasts from
this parent project were therefore exploratory advanced-degeneration-associated
displays only.

External NP cohorts were analyzed independently. GSE244889 supplied
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
fragments per kilobase of transcript per million mapped reads (FPKM) representation (4 versus 5); its paired transcripts per million (TPM) matrix was retained only as a
same-sample processing sensitivity [23-25]. GSE245147 supplied a native
Degenerated-versus-No-degenerated reads per kilobase per million mapped reads (RPKM) subset (3 versus 3) after passage and
treatment arms were excluded [26,27]. GSE167931 and GSE245147 were not pooled
together because their source family and patient-level independence could not
be resolved from public metadata.

### Locked program definition, scoring, and descriptive effect displays

The four module definitions were locked within the project before external
scoring. This was a timestamped, versioned analysis lock rather than a
prospective study registration. Supplementary Table S3 records the locking
time, source identifiers, gene lists, and SHA-256 hashes. At least 80% of the
locked genes had to map in a score. The pathway components were anchored to the
KEGG and Reactome knowledgebases [18,19] and to the cited IVDD source studies
[3,10]. A higher score means higher expression of the listed genes only; it
does not imply benefit, harm, causality, or therapeutic relevance.
For verified raw 10x matrices, source-restricted,
QC-passing cells were aggregated within donor or library and counts per million (CPM) was the per-library normalization scale, and each module score
was the mean mapped-gene `log1p(CPM)` value. For the externally supplied dense
or normalized matrices, the same score definition was applied only within the
declared matrix-processing boundary. Original archives, input hashes, and
score-to-ledger identities were retained in the reproducibility contract
(Supplementary Table S6).

For every cohort, compartment, and module, we calculated the unweighted mean
difference between groups recorded as higher versus lower severity; positive values
therefore mean higher minus lower. Welch and
within-group donor/library bootstrap intervals were displayed as descriptive
uncertainty summaries. Leave-one-key-out (LOKO) analyses assessed whether the sign
was retained; the complete default leave-one-key-out output is provided in
Supplementary Table S4. The frozen default cross-cohort sign-alignment display
was summarized without a pooled effect, formal meta-analysis, p-value, or
replication decision.

### Exploratory supplementary standardized syntheses

S7-S9 used one higher-minus-lower NP effect per cohort and module. The primary
effect was the heteroscedastic standardized mean difference (SMDH), fitted with
a restricted maximum likelihood (REML) random-effects model and Knapp-Hartung confidence intervals. Conventional
pooled-SD Hedges *g* and Paule-Mandel tau-squared estimates were model
sensitivities. All three packages recorded `metafor` REML `maxiter = 10000`,
while retaining all other `metafor` control defaults. Four module-level Hartung-Knapp-Sidik-Jonkman (HKSJ) p-values and their Benjamini-Hochberg (BH) adjustments were included solely for
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

### Computational reproducibility

The default summary used 95% Welch intervals and 10,000 independent
donor or library bootstrap draws per contrast. Cohort-, compartment-, and
module-specific deterministic seeds were derived from root seed 20260814.
Scores with a mapped-gene fraction below 80% were excluded with an audit row.
The S7-S9 syntheses used R 4.4.1, `metafor` 4.8.0, and REML control
`maxiter = 10000`. Input, output, environment, and generator hashes are
indexed in Supplementary Table S6 and the submission-support manifest.

## Results

### Cohort disposition and reproducibility

The default descriptive summary contains 20 cohort/compartment/module effects
from the four score cohorts listed in Table 1, with 55/55 exact
score-to-ledger sample-key matches (Supplementary Tables S2 and S6). All four
cohorts have `confirmatory_eligible=false`. GSE230809 contributes exploratory
AF and NP contrasts; GSE244889 contributes a 4-versus-3 NP external score-level
contrast; GSE153066 contributes an 8-versus-8 NP external dense-count contrast;
and GSE165722 contributes a 4-versus-4 NP normalized-count score-level
descriptive contrast. GSE251686 is not part of this default summary, sign
alignment, or either main figure (Supplementary Figure S2).

### Default NP directional patterns

In NP, the hypoxia/oxidative-stress module had a positive
higher-recorded-severity point estimate in each of the four default contrasts:
0.1776 in GSE230809, 0.2381 in GSE244889, 0.0882 in GSE153066, and 0.4346 in
GSE165722 (Table 2 and Figures 1-2). The corresponding descriptive Welch 95%
intervals were `[-0.0245, 0.3798]`, `[-0.1314, 0.6076]`,
`[-0.2166, 0.3929]`, and `[-0.0619, 0.9310]`; every interval included zero.
Thus, the four positive point estimates are a descriptive sign pattern, not a
robust cross-cohort replication result. ECM/collagen remodeling was positive
in GSE230809, GSE244889, and GSE165722 but negative in GSE153066.
Inflammatory/NF-κB was positive in GSE244889 and GSE153066 but negative in
GSE230809 and GSE165722. Disc-matrix homeostasis was positive in GSE230809,
GSE244889, and GSE165722 but negative in GSE153066 (Table 2 and Figures 1-2).

### Isolated GSE251686 analysis

The separately audited GSE251686 sensitivity analysis retained five presumed
sample/library keys after permanent exclusion of `GSM7986002` (mild n=2 and
severe n=3). Its hypoxia/oxidative-stress point estimate was negative
(`-0.0778`; Welch 95% interval `[-1.7214, 1.5658]`), and the four module
effects are reported in Supplementary Table S1 and Supplementary Figure S1.
This isolated result was not pooled, counted in default sign alignment, or used
to alter the four-cohort default result.

### Exploratory supplementary syntheses

The separate S7 synthesis yielded SMDH values of 1.0481 (95% CI -0.9965 to
3.0926) for ECM/collagen remodeling, 0.2809 (-1.1432 to 1.7050) for
inflammatory/NF-κB, 0.8184 (-0.1285 to 1.7654) for hypoxia/oxidative stress,
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

### Exploratory AF context and retained-cell threshold sensitivity

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

## Informed Consent Statement

No new participant consent was obtained for this secondary analysis of de-identified publicly available data. The corresponding author must confirm whether any dataset-specific consent wording is required before submission.

## Funding

`[Funding sources, grant numbers, and the funders' role to be completed before
submission.]`

## Competing Interests

`[Competing-interest declaration to be completed by all authors before
submission.]`

## Author Contributions

`[Author names and CRediT contributions to be completed and approved by all
authors before submission.]`

## Acknowledgements

[To be completed after author approval. Do not list contributors or support not confirmed by the authors.]

## References

1. Barrett T, Wilhite SE, Ledoux P, Evangelista C, Kim IF, Tomashevsky M, Marshall KA, Phillippy KH, Sherman PM, Holko M, Yefanov A, Lee H, Zhang N, Robertson CL, Serova N, Davis S, Soboleva A. NCBI GEO: archive for functional genomics data sets-update. Nucleic Acids Research. 2013;41(Database issue):D991-D995. doi:10.1093/nar/gks1193.

2. Squair JW, Gautier M, Kathe C, Anderson MA, James ND, Hutson TH, Hudelle R, Qaiser T, Matson KJE, Barraud Q, Levine AJ, La Manno G, Skinnider MA, Courtine G. Confronting false discoveries in single-cell differential expression. Nature Communications. 2021;12(1):5692. doi:10.1038/s41467-021-25960-2.

3. Swahn H, Mertens J, Olmer M, Myers K, Mondala TS, Natarajan P, Head SR, Alvarez-Garcia O, Lotz MK. Shared and Compartment-Specific Processes in Nucleus Pulposus and Annulus Fibrosus During Intervertebral Disc Degeneration. Advanced Science. 2024;11(17):e2309032. doi:10.1002/advs.202309032.

4. NCBI Gene Expression Omnibus. GSE230809: Shared and compartment-specific processes in nucleus pulposus and annulus fibrosus during intervertebral disc degeneration [Internet]. GEO Series accession GSE230809. 2024. SuperSeries containing GSE229711 and GSE230808; accessed 2026-08-14. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230809.

5. NCBI Gene Expression Omnibus. GSE229711: The cellular landscape of the healthy human intervertebral disc [Internet]. GEO Series accession GSE229711. 2024. SubSeries of GSE230809; accessed 2026-08-14. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229711.

6. NCBI Gene Expression Omnibus. GSE230808: Shared and compartment-specific processes in nucleus pulposus and annulus fibrosus during intervertebral disc degeneration [Internet]. GEO Series accession GSE230808. 2024. SubSeries of GSE230809; accessed 2026-08-14. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230808.

7. Chen F, Lei L, Chen S, Zhao Z, Huang Y, Jiang G, Guo X, Li Z, Zheng Z, Wang J. Serglycin secreted by late-stage nucleus pulposus cells is a biomarker of intervertebral disc degeneration. Nature Communications. 2024;15(1):47. doi:10.1038/s41467-023-44313-9.

8. NCBI Gene Expression Omnibus. GSE244889: Gene expression profile at single cell level of nucleus pulposus cells from mild and severe degenerative intervertebral discs [Internet]. GEO Series accession GSE244889. 2023. Accessed 2026-08-14; GEO record links PMID 38167807 and later related publications. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244889.

9. NCBI Gene Expression Omnibus. GSE153066: Single cell sequencing of human nucleus pulposus [Internet]. GEO Series accession GSE153066. 2023. Accessed 2026-08-14; GEO record explicitly lists citation as missing. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153066.

10. Tu J, Li W, Li W, Yang S, Yang P, Yan Q, Wang S, Lai K, Bai X, Wu C, Ding W, Cooper-White J, Diwan A, Yang C, Yang H, Zou J. Single-Cell Transcriptome Profiling Reveals Multicellular Ecosystem of Nucleus Pulposus during Degeneration Progression. Advanced Science. 2022;9(3):e2103631. doi:10.1002/advs.202103631.

11. NCBI Gene Expression Omnibus. GSE165722: Single-cell transcriptome profiling reveals nucleus pulposus heterogeneity and immunity during degeneration progression [Internet]. GEO Series accession GSE165722. 2021. Accessed 2026-08-14; GEO supplementary matrices are described as normalized counts. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165722.

12. Wang D, Li Z, Huang W, Cao S, Xie L, Chen Y, Li H, Wang L, Chen X, Yang JR. Single-cell transcriptomics reveals heterogeneity and intercellular crosstalk in human intervertebral disc degeneration. iScience. 2023;26(5):106692. doi:10.1016/j.isci.2023.106692.

13. Sun Y, Peng Y, Su Z, So KKH, Lu Q, Lyu M, Zuo J, Huang Y, Guan Z, Cheung KMC, Zheng Z, Zhang X, Leung VYL. Fibrocyte enrichment and myofibroblastic adaptation causes nucleus pulposus fibrosis and associates with disc degeneration severity. Bone Research. 2025;13(1):10. doi:10.1038/s41413-024-00372-2.

14. Crowell HL, Soneson C, Germain PL, Calini D, Collin L, Raposo C, Malhotra D, Robinson MD. muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. Nature Communications. 2020;11(1):6077. doi:10.1038/s41467-020-19894-4.

15. Niu H, Qi H, Zhang P, Meng H, Liu N, Zhang D. Single-Cell Analysis Reveals Aspirin Restores Intervertebral Disc Integrity via Ferroptosis Regulation. Journal of Inflammation Research. 2025;18:6889-6905. doi:10.2147/JIR.S519218.

16. Jia S, Liu H, Yang T, Gao S, Li D, Zhang Z, Zhang Z, Gao X, Liang Y, Liang X, Wang Y, Meng C. Single-cell sequencing reveals cellular heterogeneity of nucleus pulposus in intervertebral disc degeneration. Scientific Reports. 2024;14(1):27245. doi:10.1038/s41598-024-78675-x.

17. NCBI Gene Expression Omnibus. GSE251686: Single-cell sequencing reveals cellular heterogeneity of nucleus pulposus in intervertebral disc degeneration [Internet]. GEO Series accession GSE251686. 2024. Accessed 2026-08-14; GSM7986002 is excluded in the present audit for failed stream integrity. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251686.

18. Kanehisa M, Furumichi M, Sato Y, Kawashima M, Ishiguro-Watanabe M. KEGG for taxonomy-based analysis of pathways and genomes. Nucleic Acids Research. 2023;51(D1):D587-D592. doi:10.1093/nar/gkac963.

19. Jassal B, Matthews L, Viteri G, Gong C, Lorente P, Fabregat A, Sidiropoulos K, Cook J, Gillespie M, Haw R, Loney F, May B, Milacic M, Rothfels K, Sevilla C, Shamovsky V, Shorser S, Varusai T, Weiser J, Wu G, Stein L, Hermjakob H, D'Eustachio P. The Reactome pathway knowledgebase. Nucleic Acids Research. 2020;48(D1):D498-D503. doi:10.1093/nar/gkz1031.

20. Welch BL. The generalization of Student's problem when several different population variances are involved. Biometrika. 1947;34(1-2):28-35. doi:10.1093/biomet/34.1-2.28.

21. Efron B. Bootstrap methods: another look at the jackknife. The Annals of Statistics. 1979;7(1):1-26. doi:10.1214/aos/1176344552.

22. NCBI Gene Expression Omnibus. GSE186542: Nucleus pulposus related lncRNA and mRNA expression profiles in intervertebral disc degeneration [Internet]. GEO Series accession GSE186542. 2021. Submitted 2021; accessed 2026-08-14; GEO SOFT lists no linked PubMed citation. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186542.

23. NCBI Gene Expression Omnibus. GSE167931: Next Generation Sequencing analysis at single-cell level of normal and degenerated nucleus pulposus cells transcriptomes [Internet]. GEO Series accession GSE167931. 2021. Accessed 2026-08-14; GEO SOFT lists PMIDs 35304463 and 35340126. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167931.

24. Li G, Ma L, He S, Luo R, Wang B, Zhang W, Song Y, Liao Z, Ke W, Xiang Q, Feng X, Wu X, Zhang Y, Wang K, Yang C. WTAP-mediated m6A modification of lncRNA NORAD promotes intervertebral disc degeneration. Nature Communications. 2022;13(1):1469. doi:10.1038/s41467-022-28990-6.

25. Li G, Luo R, Zhang W, He S, Wang B, Liang H, Song Y, Ke W, Shi Y, Feng X, Zhao K, Wu X, Zhang Y, Wang K, Yang C. m6A hypomethylation of DNMT3B regulated by ALKBH5 promotes intervertebral disc degeneration via E4F1 deficiency. Clinical and Translational Medicine. 2022;12(3):e765. doi:10.1002/ctm2.765.

26. NCBI Gene Expression Omnibus. GSE245147: CytoDNA triggered NP cell inflammatory senescence via cGAS-STING axis sensing but not AIM2 inflammasome activation [Internet]. GEO Series accession GSE245147. 2024. Accessed 2026-08-14; GEO SOFT lists PMID 38488012. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245147.

27. Zhang W, Li G, Zhou X, Liang H, Tong B, Wu D, Yang K, Song Y, Wang B, Liao Z, Ma L, Ke W, Zhang X, Lei J, Lei C, Feng X, Wang K, Zhao K, Yang C. Disassembly of the TRIM56-ATR complex promotes cytoDNA/cGAS/STING axis-dependent intervertebral disc inflammatory degeneration. The Journal of Clinical Investigation. 2024;134(6):e165140. doi:10.1172/JCI165140.

## Supplementary Table and Figure Legends

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
synthesis.** SMDH denotes the heteroscedastic standardized mean difference, and REML
denotes restricted maximum likelihood. The tables give S7 study-level effects,
primary SMDH/REML/Knapp-Hartung results, Hedges *g* and Paule-Mandel sensitivities, and
leave-one-cohort-out results. They are a separate non-confirmatory synthesis of
the four default NP cohorts and do not replace the default descriptive analysis.

**Supplementary Tables S8a-S8d. Post hoc six-cohort NP expansion.** SMDH denotes the heteroscedastic standardized mean difference, and REML
denotes restricted maximum likelihood. The tables give the separately packaged S8
study-level effects, primary results, model
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
synthesis.** The forest display gives S7 SMDH estimates and Knapp-Hartung intervals;
SMDH is a heteroscedastic standardized mean difference and the display is
non-confirmatory. It is a separate exploratory quantitative summary and not a
confirmation or replication display.

**Supplementary Figure S4. Post hoc six-cohort NP expansion.** The forest display gives S8 results after adding GSE186542 and GSE167931 FPKM. The plotted SMDH estimates and intervals are non-confirmatory; SMDH is a
heteroscedastic standardized mean difference, and HKSJ/BH p-values are reported
for transparency only.

**Supplementary Figure S5. NP source-family replacement sensitivity.** The
forest display gives S9 results after replacing GSE167931 with the native
GSE245147 subset. It demonstrates source-family sensitivity and does not add an independent
validation cohort. Axis ranges may differ among panels and effect magnitudes are
not directly comparable across processing scales.

## Main Tables

**Table 1. Cohort roles and inference boundaries in the default descriptive summary.**

[[TABLE:1]]

**Table 2. Cohort-specific NP module effects in the default descriptive summary.** The default layer contains 20 effects overall (16 NP effects and 4 exploratory AF effects); this table displays the 16 NP effects.

[[TABLE:2]]

## Main Figures

**Figure 1. Cohort-specific NP module-score differences and descriptive Welch 95% intervals.** Colors distinguish cohorts only; they do not encode effect magnitude.

[[FIGURE:figure_1]]

**Figure 2. NP cohort-specific directions and descriptive sign alignment.** Blue and orange encode positive and negative directions only; color does not encode effect magnitude.

[[FIGURE:figure_2]]

## Supplementary Material

Supplementary Tables S1-S4, S5a-S5b, S6, S7a-S7d, S8a-S8d, and S9a-S9d, together with Supplementary Figures S1-S5, are provided as separate submission files. Each S7-S9 artifact remains a separately packaged, exploratory, non-confirmatory synthesis. The graphical abstract should be uploaded separately if requested by the target journal.
