# GSE165722 Score-Level Module Readiness

Audit date: 2026-08-14.

## Purpose and boundary

This document records the only permitted GSE165722 expression use in the
current project: pre-specified, sample-level module-score direction checks for
NP severe versus mild degeneration. It is not a raw-count analysis,
pseudobulk count analysis, or cell-level analysis.

GEO describes the supplementary matrices as normalized counts. Their
integer-like appearance does not change that status. The supplied values are
therefore never called UMI counts, and no edgeR, DESeq2, or other
negative-binomial model is fitted to this cohort.

## Inputs and integrity checks

| Item | Result |
|---|---|
| Archive | `data/raw/GSE165722_RAW.tar`, retained unmodified |
| Samples | 8 presumed donor-level Sample/GSM keys; mild = 4, severe = 4 |
| Tissue | NP for every sample |
| Matrix layout | Dense gene-by-cell TSV within TAR, streamed without extraction |
| Matrix-to-cell mapping | Every `CellIndex` vector exactly matched the corresponding count-matrix header |
| Supplied cells | 49,637 total; 5,073 to 7,066 per sample |
| Genes streamed | 16,226 to 18,933 per sample |
| Value boundary | Non-negative and integer-like across every streamed supplied value; still treated as normalized values per GEO |
| Locked modules | Four modules in `config/program_modules.json` |
| Mapping rule | At least 80% of configured module genes required |

All eight samples generated all four module scores (32 score rows). Every
score was available at the locked mapping threshold. Per-sample mapping and
duplicate-feature audits are in
`data/derived/module_scores_external/GSE165722/GSE165722_RAW_module_mapping_audit.csv`.

## Score definition

For sample `s` and module `M`, after summing each mapped gene over the
sample's supplied cell columns, the score is:

`mean_g in M [ log1p(1,000,000 * supplied_value_sum(g, s) / total_supplied_value(s)) ]`.

The score is a transformed within-sample program summary. It does not claim
that the denominator or numerator is a raw molecule count, and it is not
comparable as an absolute abundance measurement across platforms.

The analysis unit is one presumed donor-level Sample/GSM key. The cells used
to calculate a key's score remain nested measurements and are never resampled
or tested as biological replicates.

## Severity map and metadata limitation

GEO's ordered grade labels I-IV disagree with the source publication's Table
1 labels II-V. The locked severity map follows Tu et al. (PMID 34825784):

| Samples | Source-publication grade | Score-level group |
|---|---:|---|
| Sample1-Sample2 | II | mild |
| Sample3-Sample4 | III | mild |
| Sample5-Sample6 | IV | severe |
| Sample7-Sample8 | V | severe |

GEO does not expose an independent patient ID, age, sex, disc level, or batch
field. A Sample/GSM is therefore a presumed donor-level key, not a verified
patient identifier. Results must retain this limitation and cannot be read as
age-adjusted, sex-adjusted, causal, clinical, or therapeutic evidence.

## Generated artifacts

- `GSE165722_RAW_module_scores.csv`: compatible module-score table for a
  descriptive severe-minus-mild direction check.
- `GSE165722_RAW_sample_score_ledger.csv`: per-sample stream and identity
  checks.
- `GSE165722_RAW_module_mapping_audit.csv`: mapped-gene and duplicate-feature
  audit.
- `GSE165722_RAW_module_gene_scores.csv`: module-gene supplied-value sums,
  retained solely for audit.
- `GSE165722_RAW_module_score_parameters.csv`: fixed formula and boundaries.

## Reproduction

```powershell
& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\score_gse165722_processed_modules.py" `
  ".\data\raw\GSE165722_RAW.tar" `
  ".\data\derived\GSE165722_donor_ledger.csv" `
  ".\config\program_modules.json" `
  --output-dir ".\data\derived\module_scores_external\GSE165722"
```

## Interpretation rule

This cohort may contribute a descriptive severe-minus-mild score direction,
effect-size interval, and leave-one-presumed-key-out display. It is not
confirmatory validation and cannot supply raw-count inference. A sign match
with another cohort is descriptive only and does not establish replication or
a universal IVDD program.
