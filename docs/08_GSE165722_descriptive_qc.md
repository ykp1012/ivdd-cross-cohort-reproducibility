# GSE165722 Processed-Matrix Descriptive QC

Audit date: 2026-08-14. Source table: `data/derived/GSE165722_descriptive_qc.csv`.

## Scope

These are descriptive summaries of GEO-supplied matrices, which GEO labels
"normalized counts". The word "supplied value" is used deliberately: these
summaries do not establish raw UMI status and are not inputs to edgeR, DESeq2,
pseudobulk count models, cell filtering, or biological group inference.

## Matrix-level summary

- Eight presumed donor-level NP matrices contain 49,637 supplied cells.
- Per-matrix cell counts range from 5,073 to 7,066.
- Per-matrix median supplied values range from 874.0 to 4033.0.
- Per-matrix median detected features range from 348.0 to 1417.0.
- Per-matrix median mitochondrial fractions range from 1.9% to 7.9%.
- The fraction meeting source-like descriptive cutoffs (supplied value >=500,
  detected features >=200, mitochondrial fraction <=20%) ranges from
  64.8% to 94.5%.

## Interpretation limits

The range in supplied-library characteristics requires transparent reporting
and later sensitivity analyses. It is not evidence for a disease effect,
cell-type shift, or sample exclusion. The cohort remains eligible only for
score-level directional support after gene sets and scoring rules are locked.
