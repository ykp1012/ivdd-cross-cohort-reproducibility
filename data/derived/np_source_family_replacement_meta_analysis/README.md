# Source-Family Replacement Sensitivity: Exploratory NP Random-Effects Meta-Analysis

This package is separate from the frozen four-cohort S7/S3 analysis.
It replaces the GSE167931 FPKM representation with the native clinical-comparison subset of GSE245147.

## Scope and boundary

- GSE230809 remains one parent project and is not split into GSE229711 and GSE230808 studies.
- GSE251686 and GSM7986002 remain outside this package.
- Six human NP cohort/module contrasts are included once each; k = 6 per module.
- GSE245147 includes only Degenerated_1-3 versus NO_Degenerated_1-3; P2/P8 and DMSO/H-151 arms are excluded.
- GSE167931 and GSE245147 are not pooled together because their source lab/author family overlaps and patient-level reuse cannot be excluded.
- All included contrasts remain post hoc, score-level, and confirmatory_eligible=false.

## Methods

The primary standardized effect is metafor SMDH, which permits unequal group variances.
A conventional pooled-SD Hedges g analysis and Paule-Mandel tau-squared estimate are sensitivity analyses.
All models use random effects with Knapp-Hartung intervals.
REML Fisher scoring uses maxiter=10000 with all other metafor control values left at default; this changes only the iteration ceiling, not the estimator or interval method.
With k=6 per module and heterogeneous source material, prediction intervals, Q tests, I2, and leave-one-cohort-out results are descriptive and unstable.
No meta-regression, funnel plot, Egger test, causal claim, patient-level replication claim, biomarker claim, or therapeutic inference is permitted.

## Primary SMDH results
- ECM / collagen remodeling: pooled SMDH 1.0931, 95% CI [-0.9186, 3.1048], prediction interval [-2.2873, 4.4734], I2 56.8%.
- Inflammatory / NF-kB: pooled SMDH 0.4056, 95% CI [-1.0004, 1.8115], prediction interval [-1.4342, 2.2453], I2 20.6%.
- Hypoxia / oxidative stress: pooled SMDH 0.5746, 95% CI [-0.7231, 1.8723], prediction interval [-0.7231, 1.8723], I2 0.0%.
- Disc matrix homeostasis: pooled SMDH 0.1046, 95% CI [-1.0200, 1.2292], prediction interval [-1.0200, 1.2292], I2 0.0%.

This is a source-family replacement sensitivity, not an independent validation analysis. It tests whether the pooled pattern depends on the choice of one potentially related external cohort.
