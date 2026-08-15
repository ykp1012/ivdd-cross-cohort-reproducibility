# Post Hoc External Expansion: Exploratory NP Random-Effects Meta-Analysis

This package is separate from the frozen four-cohort S7/S3 analysis.
It adds GSE186542 and the GSE167931 FPKM representation after candidate GEO auditing.

## Scope and boundary

- GSE230809 remains one parent project and is not split into GSE229711 and GSE230808 studies.
- GSE251686 and GSM7986002 remain outside this package.
- Six human NP cohort/module contrasts are included once each; k = 6 per module.
- GSE167931 TPM is a same-sample processing sensitivity, not an additional study.
- GSE186542 and GSE167931 are accession/BioProject-level additions; patient-level overlap cannot be excluded from public metadata.
- All included contrasts remain post hoc, score-level, and confirmatory_eligible=false.

## Methods

The primary standardized effect is metafor SMDH, which permits unequal group variances.
A conventional pooled-SD Hedges g analysis and Paule-Mandel tau-squared estimate are sensitivity analyses.
All models use random effects with Knapp-Hartung intervals.
REML Fisher scoring uses maxiter=10000 with all other metafor control values left at default; this changes only the iteration ceiling, not the estimator or interval method.
With k=6 per module and heterogeneous source material, prediction intervals, Q tests, I2, and leave-one-cohort-out results are descriptive and unstable.
No meta-regression, funnel plot, Egger test, causal claim, patient-level replication claim, biomarker claim, or therapeutic inference is permitted.

## Primary SMDH results
- ECM / collagen remodeling: pooled SMDH 0.7780, 95% CI [-0.3532, 1.9093], prediction interval [-1.5494, 3.1055], I2 48.8%.
- Inflammatory / NF-kB: pooled SMDH 0.4032, 95% CI [-0.3991, 1.2055], prediction interval [-0.5818, 1.3882], I2 7.0%.
- Hypoxia / oxidative stress: pooled SMDH 0.7694, 95% CI [0.1706, 1.3682], prediction interval [0.1706, 1.3682], I2 0.0%.
- Disc matrix homeostasis: pooled SMDH 0.3762, 95% CI [-0.3017, 1.0540], prediction interval [-0.3017, 1.0540], I2 0.0%.

The expansion quantifies how two newly audited public datasets shift the cross-cohort pattern. It does not establish a universal IVDD program or an independent patient-level validation.
