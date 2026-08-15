# Discovery Retained-Cell Threshold Sensitivity

Updated: 2026-08-14.

## Question and unit of analysis

This pre-specified sensitivity analysis asks whether applying a minimum of 20,
30 (primary), or 50 source-restricted, QC-passing cells per discovery
donor/library changes the inclusion set or donor/library-level program-score
contrast. It uses the combined GSE230809 parent project (`GSE229711` healthy
and `GSE230808` diseased); the child series are not independent cohorts.

Cells are nested observations. No cell-level test, p-value, or cell-level
resampling was used. For every eligible donor/library, all source-restricted,
QC-passing cells were aggregated into the pseudobulk; a threshold is an
eligibility gate, not a request to retain exactly 20, 30, or 50 cells.

## Reproducible runs

For each value in `20`, `30`, and `50`, both original archives were rescored
with `scripts/score_module_pseudobulk.py`, the shared discovery raw-data
ledger, the corresponding fixed QC and annotation ledgers, the locked module
configuration, and `--min-retained-cells <value>`. The exact threshold-specific
score paths are in:

- `config/discovery_retained_cell_sensitivity_threshold_20.csv`
- `config/discovery_retained_cell_sensitivity_threshold_30.csv`
- `config/discovery_retained_cell_sensitivity_threshold_50.csv`

Each score run was summarized separately with
`scripts/summarize_donor_module_effects.py`. It reports unweighted
donor/library mean-difference and leave-one-key-out descriptive outputs for
AF and NP; it does not calculate a confirmatory test.

## Result

All 24 discovery donor/libraries passed every threshold. The smallest observed
number of source-restricted, QC-passing cells was 471. Therefore the 20-,
30-, and 50-cell thresholds selected the same 24 libraries. The 96 scored
library-module rows were exactly identical to the 30-cell reference at each
threshold (maximum absolute score delta = 0), and all eight AF/NP
module-effect rows had identical sample sizes, effects, and directions.

This shows that the conclusions in this limited threshold range are unchanged
because no library lies near a threshold boundary. It does **not** demonstrate
robustness to composition changes, annotation changes, or random cell
downsampling; those require distinct analyses.

## Outputs

The machine-readable final outputs are under
`data/derived/discovery_retained_cell_sensitivity/`:

- `library_threshold_eligibility.csv`: every threshold x donor/library gate.
- `threshold_run_summary.csv`: high-level inclusion accounting.
- `threshold_score_identity.csv`: exact module-score comparison against the
  primary 30-cell result.
- `discovery_effects_by_retained_cell_threshold.csv`: all descriptive AF/NP
  donor/library effects with 10,000 donor/library bootstrap draws per row.
- `threshold_effect_stability_vs_30.csv`: effect and direction comparison to
  the 30-cell threshold.
- `input_artifact_hashes.csv`: SHA-256 hashes of threshold-specific inputs.

The threshold-specific score and effect-summary directories remain under
`data/derived/module_scores_sensitivity/threshold_{20,30,50}/` and each
contains its own parameter and artifact audit files.

## Interpretation boundary

The GSE230809 contrast remains exploratory: healthy and diseased donors are
age-confounded, the healthy group contains only three donors per compartment,
and AF lacks an independent external endpoint. Threshold invariance does not
upgrade any result to replication, confirmation, causal interpretation, or an
age-independent degeneration association.
