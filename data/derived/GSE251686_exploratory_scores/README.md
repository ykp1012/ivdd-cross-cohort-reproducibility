# GSE251686 isolated exploratory effect summary

This directory is intentionally separate from the default `data/derived/donor_module_effect_summary/` result. It contains a descriptive mild n=2 versus severe n=3 comparison of the five stream-integrity-passing presumed sample/library keys. `GSM7986002` is permanently excluded because its Matrix Market payload failed the independent stream-integrity audit.

The summary uses unweighted severe-minus-mild differences, descriptive Welch 95% intervals, 10000 within-arm presumed-key bootstrap intervals, and leave-one-key-out direction checks. These are not hypothesis tests: no p-values, multiple-testing adjustment, formal meta-analysis, validation, replication adjudication, causal claim, biomarker claim, or therapeutic claim is produced.

Every score row must meet the locked mapped-gene fraction of 0.800. The scorer records a 30-cell primary eligibility gate; this effect summary additionally requires every selected key to pass the 20-, 30-, and 50-cell source-restricted gates. Input hashes are recorded in `GSE251686_exploratory_effect_parameters.csv`; the run manifest records exact generated-artifact hashes. `run_artifacts.csv` and `run_manifest.json` are excluded from generated-artifact hashing because each would otherwise be self-referential.
