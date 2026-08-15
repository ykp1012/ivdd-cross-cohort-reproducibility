# Topic Brief

## Working title
Cohort-aware direction and heterogeneity of locked nucleus pulposus transcriptional programs across human intervertebral disc degeneration datasets

## Slug
`manuscript`

## Core question
- How do pre-specified IVDD-related NP transcriptional programs align or differ in donor-level direction between recorded lower- and higher-severity groups across independently processed public cohorts?

## Primary hypothesis
- The analysis will describe sign alignment and discordance of each pre-specified NP program across independently processed cohorts. Any aligned sign pattern remains descriptive and will be reported with cohort-specific limitations rather than as a universal IVDD signature.

## Data assets
- Discovery: GSE230809 parent project (GSE229711 + GSE230808), human AF/NP raw 10x matrices; AF analyses are exploratory.
- External NP support: GSE165722 (score-level direction evidence only because GEO describes integer matrices as normalized), GSE153066 (audited count-level normal-versus-degenerated support), and GSE244889 (small directional support). GSE251686 has a separate audit-gated mild n=2 versus severe n=3 exploratory score/effect package after permanent exclusion of malformed `GSM7986002`; it is deliberately excluded from the default 20-effect summary and is not a confirmatory cohort. After the default package was frozen, S7 independently synthesized the four default NP cohorts; S8 added GSE186542 and GSE167931 FPKM; and S9 replaced GSE167931 with the native GSE245147 clinical-comparison subset because source-family overlap cannot be excluded. All three packages are non-confirmatory.

## Core design
- Donor or presumed donor/sample key is the unit of observation. Each cohort is processed separately. Raw-UMI pseudobulk scoring is used only where raw counts are verified. The module ledger and scoring rule are locked before external scoring. The frozen default layer displays cohort-specific effect estimates, intervals, and leave-one-key-out stability without a pooled estimate, p-value, or replication adjudication. S7-S9 separately standardize one NP effect per cohort with SMDH random-effects REML and Knapp-Hartung intervals; HKSJ and BH p-values are transparent descriptors only.

## Primary outcome
- Descriptive per-cohort direction, effect size, interval, and leave-one-key-out stability for each pre-specified NP module in the frozen default analysis; separate non-confirmatory S7-S9 standardized synthesis results are supplementary outcomes.

## Secondary outcomes
- Discovery-only AF/NP compartment contrasts, donor-level QC, gene-level concordance, and heterogeneity/sensitivity analyses.

## Must-have analyses
- Provenance and sample ledger; raw/processed status audit; donor-level QC; pseudobulk only with raw UMI; locked module ledger; independent cohort effect estimates; donor bootstrap, leave-one-donor-out, and retained-cell-threshold sensitivity analyses; descriptive cross-cohort sign alignment; isolated SMDH/REML/Knapp-Hartung S7-S9 analyses with Hedges g and Paule-Mandel sensitivities, four-module BH reporting, and source-family replacement.

## Optional extensions
- Healthy atlas annotation reference (GSE160756); descriptive one-donor direction check (GSE199866). Neither is confirmatory support.

## Journal or audience
- Orthopaedic and spine research audience; an observational transcriptomic reproducibility study, not a mechanistic or therapeutic paper.
