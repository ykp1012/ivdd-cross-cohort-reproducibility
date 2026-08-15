# Exploratory NP Random-Effects Meta-Analysis

This package adds a cohort-level exploratory quantitative synthesis to the existing IVDD project.
The pre-existing donor/library-level Welch mean differences remain the primary analysis.

## Scope and boundary

- Four default NP cohorts are included: GSE230809 parent project, GSE244889, GSE153066, and GSE165722.
- GSE229711 and GSE230808 remain one GSE230809 parent project, never two independent studies.
- GSE251686 and GSM7986002 remain outside this primary meta-analysis package.
- One effect is used per cohort and module; AF is not pooled with NP.
- All included cohorts retain confirmatory_eligible=false.

## Methods

The primary standardized effect is metafor SMDH, which allows unequal group variances.
A conventional pooled-SD Hedges g analysis and Paule-Mandel tau-squared estimate are sensitivity analyses.
All models use REML or Paule-Mandel random effects with Knapp-Hartung intervals.
REML Fisher scoring uses maxiter=10000 with all other metafor control values left at default; this changes only the iteration ceiling, not the estimator or interval method.
With k=4 per module, heterogeneity estimates, Q tests, prediction intervals, and leave-one-cohort-out results are descriptive and unstable.
No meta-regression, funnel plot, Egger test, causal claim, validation claim, or therapeutic inference is permitted.

## Primary SMDH results
- ECM / collagen remodeling: pooled SMDH 1.0481, 95% CI [-0.9965, 3.0926], prediction interval [-2.9918, 5.0879], I2 65.1%.
- Inflammatory / NF-kB: pooled SMDH 0.2809, 95% CI [-1.1432, 1.7050], prediction interval [-1.9452, 2.5070], I2 31.7%.
- Hypoxia / oxidative stress: pooled SMDH 0.8184, 95% CI [-0.1285, 1.7654], prediction interval [-0.1285, 1.7654], I2 0.0%.
- Disc matrix homeostasis: pooled SMDH 0.1453, 95% CI [-0.3679, 0.6585], prediction interval [-0.3679, 0.6585], I2 0.0%.

These results quantify the current cross-cohort pattern; they do not establish a universal IVDD program or a confirmatory biological conclusion.
