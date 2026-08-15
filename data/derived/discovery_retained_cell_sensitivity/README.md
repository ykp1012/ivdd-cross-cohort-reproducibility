# Discovery Retained-Cell Threshold Sensitivity

This directory summarizes the pre-specified 20, 30, and 50 retained-cell
threshold analysis for the combined GSE230809 discovery parent project
(GSE229711 plus GSE230808). The threshold is a donor/library eligibility rule:
each source-restricted, QC-passing donor/library must contain at least the
specified number of cells before its full pseudobulk is scored. It is not a
random cell downsampling analysis and it does not make cells independent
replicates.

For each threshold, `score_module_pseudobulk.py` was rerun from the original
10x TAR archives using the same locked module configuration, QC ledger, and
annotation ledger. `summarize_donor_module_effects.py` then calculated only
unweighted donor/library target-minus-comparison score differences and
leave-one-donor/library-out displays for AF and NP separately.

Read `threshold_run_summary.csv` before interpreting the effect tables.
`threshold_score_identity.csv` compares every threshold-specific module score
to the primary 30-cell run. `threshold_effect_stability_vs_30.csv` compares
the resulting donor/library effects and directions. All effects remain
exploratory because GSE230809 is one parent project with three healthy donors
per compartment and complete age-disease confounding.
