# Submission Outline

## Working Title

Cohort-aware direction and heterogeneity of locked nucleus pulposus
transcriptional programs across human intervertebral disc degeneration datasets

## Central Question

Across independently processed public NP datasets, do four project-locked
expression programs show descriptively aligned or discordant directions between
recorded lower- and higher-severity groups when the observation unit is a donor
or explicitly labelled presumed sample/library key rather than an individual
cell?

## Narrative Structure

### Abstract

The abstract should frame the analysis as a public-data, cohort-aware evidence
audit. It should state that cells are nested observations, that the four program
definitions were locked before external scoring, and that the frozen default
20-effect analysis has no p-values, pooled estimate, or replication
adjudication. It should then distinguish S7, S8, and S9 as separately packaged,
non-confirmatory random-effects syntheses with transparent HKSJ/BH p-values.
Its central default result is that all four NP hypoxia/oxidative-stress point
estimates were positive while all four corresponding Welch intervals included
zero. It should also state that S8's hypoxia interval did not cross zero but
that S9 lost this pattern after the GSE167931-to-GSE245147 replacement.

### Introduction

The introduction should distinguish cell abundance from independent biological
observations, position the work alongside prior public IVDD single-cell studies,
and explain that the contribution is an audit of cohort-specific directions and
their limitations rather than a claim of molecular novelty or a universal IVDD
signature.

### Methods

The methods should document cohort provenance and the experimental-unit
boundary, treating GSE229711 and GSE230808 as the one GSE230809 parent project.
It should explain which cohorts enter the default 20-effect summary, why
GSE251686 is isolated after stream-integrity failure of GSM7986002, the
versioned module lock, mapping threshold, score definitions, descriptive
intervals, and leave-one-key-out stability. It should then specify S7-S9:
SMDH, REML, Knapp-Hartung intervals, Hedges g and Paule-Mandel sensitivities,
four-module BH reporting, REML `maxiter = 10000` with all other controls at
default, and the GSE245147-for-GSE167931 source-family replacement. It should
state that the project lock is not a prospective registration and refer readers
to Supplementary Tables S2-S9 and Supplementary Figures S2-S5.

### Results

The results should first establish the four-cohort default summary, its 20
effects, 55/55 score-to-ledger matches, and non-confirmatory boundary. It should
then report the four hypoxia/oxidative-stress point estimates together with the
fact that every Welch interval includes zero. ECM, inflammatory/NF-kB, and
disc-matrix-homeostasis directions should be reported as discordant. The
separate GSE251686 mild n=2 versus severe n=3 sensitivity result should be
reported transparently but excluded from the main sign-alignment statement.
The results should then report S7, S8, and S9 in separate paragraphs, including
the S8 hypoxia SMDH 0.7694 (95% CI 0.1706 to 1.3682) and the S9 replacement
SMDH 0.5746 (-0.7231 to 1.8723). Exploratory AF context and retained-cell
threshold stability should be routed to Supplementary Tables S4-S5b rather
than promoted as independent main results.

### Discussion and Conclusion

The discussion should interpret aligned point-estimate signs as descriptive,
not as replication. It should foreground the wide intervals, cohort-specific
processing, age and clinical-source confounding, uncertain sample-key identity,
unverified patient-level independence, source-family sensitivity, and lack of
an independent AF severity cohort. The conclusion should state that the project
supports cohort-specific reporting rather than a causal, mechanistic,
prognostic, or therapeutic IVDD claim.

## Main Display Plan

| Item | Current file | Manuscript role |
|---|---|---|
| Graphical abstract | results/graphical_abstract/graphical_abstract_cohort_aware_ivdd.pdf and .png | Evidence-bound visual summary of the cohort-aware workflow and the limits of the conclusion. |
| Table 1 | results/tables/table_1_current_cohort_roles.csv | Cohort roles, observation keys, identity checks, and interpretation boundaries for the four default score cohorts. |
| Table 2 | results/tables/table_2_np_module_effects.csv | Exact cohort-specific NP effects, descriptive intervals, and leave-one-key-out retention. |
| Figure 1 | results/figures/figure_1_np_cohort_module_effects.pdf and .png | Four-panel display of default NP module effects and descriptive Welch intervals. |
| Figure 2 | results/figures/figure_2_np_direction_alignment.pdf and .png | Sign-only alignment display for the default NP contrasts. |

## Supplementary Display Plan

| Item | Current file | Purpose |
|---|---|---|
| Supplementary Table S1 | results/supplementary_tables/supplementary_table_s1_gse251686_exploratory_effects.csv | Isolated GSE251686 exploratory sensitivity effects. |
| Supplementary Table S2 | results/supplementary_tables/supplementary_table_s2_cohort_disposition_and_identity.csv | Cohort disposition, identity checks, and default-summary eligibility. |
| Supplementary Table S3 | results/supplementary_tables/supplementary_table_s3_locked_program_modules.csv | Locked genes, source identifiers, timestamps, and hashes. |
| Supplementary Table S4 | results/supplementary_tables/supplementary_table_s4_default_leave_one_key_out.csv | Full default leave-one-key-out output. |
| Supplementary Table S5a | results/supplementary_tables/supplementary_table_s5a_discovery_threshold_summary.csv | Retained-cell threshold eligibility summary. |
| Supplementary Table S5b | results/supplementary_tables/supplementary_table_s5b_discovery_threshold_effect_stability.csv | Effect stability across retained-cell thresholds. |
| Supplementary Table S6 | results/supplementary_tables/supplementary_table_s6_reproducibility_contract.csv | Reproducibility files and SHA-256 contract. |
| Supplementary Tables S7a-S7d | results/supplementary_tables/supplementary_table_s7*_np_meta_analysis_*.csv | Independent four-cohort exploratory SMDH synthesis, model sensitivity, and leave-one-cohort-out results. |
| Supplementary Tables S8a-S8d | results/supplementary_tables/supplementary_table_s8*_np_post_hoc_external_expansion_*.csv | Post hoc six-cohort expansion with GSE186542 and GSE167931 FPKM. |
| Supplementary Tables S9a-S9d | results/supplementary_tables/supplementary_table_s9*_np_source_family_replacement_sensitivity_*.csv | GSE245147-for-GSE167931 source-family replacement sensitivity. |
| Supplementary Figure S1 | results/supplementary_figures/supplementary_figure_s1_gse251686_exploratory_sensitivity.pdf and .png | Isolated effect and stability display for GSE251686. |
| Supplementary Figure S2 | results/supplementary_figures/supplementary_figure_s2_cohort_disposition_and_analysis_boundary.pdf and .png | Cohort disposition and hard analysis boundary. |
| Supplementary Figure S3 | results/supplementary_figures/supplementary_figure_s3_np_exploratory_random_effects_meta_analysis.pdf and .png | S7 four-cohort non-confirmatory synthesis. |
| Supplementary Figure S4 | results/supplementary_figures/supplementary_figure_s4_np_post_hoc_external_expansion_meta_analysis.pdf and .png | S8 post hoc six-cohort expansion. |
| Supplementary Figure S5 | results/supplementary_figures/supplementary_figure_s5_np_source_family_replacement_sensitivity.pdf and .png | S9 source-family replacement sensitivity. |

## Claims Excluded

No title, abstract, figure legend, cover letter, or discussion paragraph should
claim a replicated molecular signature, universal IVDD program, age-independent
disease effect, causal pathway, biomarker, diagnostic classifier, treatment
target, AF external validation, or biological confirmation. The default layer
must not claim a formal meta-analysis or hypothesis test. S7-S9 may be named as
non-confirmatory exploratory syntheses, but their HKSJ/BH p-values cannot be
presented as a test of confirmation. No main Figure 3 is planned; exploratory
AF and threshold-sensitivity details remain supplementary.
