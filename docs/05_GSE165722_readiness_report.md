# GSE165722 Data-Quality and Analysis-Readiness Report

Audit date: 2026-08-13.

## Verified facts

- The GEO archive contains 8 presumed donor-level NP matrices and 49,637 supplied cells in total.
- The cell-name files map supplied matrix column IDs to cell barcodes. Each GSM/Sample file is a presumed donor-level sample, but GEO does not provide a patient identifier, age, sex, disc level, or batch field; this nesting remains an explicit limitation.
- The source publication, Tu et al. (PMID 34825784), reports grades II-V. GEO SOFT labels the same ordered samples I-IV. This project preserves both fields and uses the source-publication grouping: mild II-III = 4 donors, severe IV-V = 4 donors.
- GEO sample metadata states that the supplementary files contain "normalized counts", even though inspected values are integer-like.

## Permitted role

GSE165722 is eligible for descriptive per-donor QC and pre-specified score-level direction checks. It is **not** eligible for raw-count pseudobulk aggregation, edgeR/DESeq2 inference, or a negative-binomial effect estimate unless independently verified raw UMI matrices are recovered.

## Stop conditions already resolved

- Donor mapping: provisional pass at sample level; one GSM/Sample file is one supplied matrix, but its biological-donor identity is not independently exposed in GEO metadata.
- Tissue label: pass; all samples are NP.
- Severity grouping: usable only with the cited source-publication mapping, not the GEO labels alone.
- Raw-count status: fail for primary count-model inference because GEO explicitly describes the values as normalized.

## Remaining before biological interpretation

1. Calculate matrix-level and cell-level descriptive QC from the supplied values without calling them UMI counts.
2. Lock module gene lists and score transformations before looking at group effects.
3. Use GSE165722 only as one component of cross-cohort directional support; report its grade conflict and processing limitation in all methods and supplements.
