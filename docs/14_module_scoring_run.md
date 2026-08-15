# Discovery Module-Scoring Run Record

Updated: 2026-08-14.

## Scope

This run aggregates raw 10x UMI counts over source-restricted, QC-passing cells
and scores the four locked program modules. A donor/library is the unit of
inference; cells remain nested observations. The GSE229711 and GSE230808 child
series are one GSE230809 discovery project and are not independent validation
cohorts.

## Input contract

For each active archive, the scorer requires:

- the shared discovery raw-data ledger;
- the matching per-cell QC TSV emitted by `qc_10x_archives.py`;
- the matching per-cell annotation ledger emitted by `annotate_10x_qc_cells.py`;
- the locked `config/program_modules.json`.

The archive GSM set must be contained in the shared ledger. Every archive
barcode must have exactly one QC row, and every QC-passing barcode must have
exactly one annotation row. Duplicate GSMs, barcodes, or donor-compartment
libraries are rejected before scores are written.

## Primary command

From `ivdd_cross_cohort_reproducibility/`:

```powershell
& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\score_module_pseudobulk.py" `
  ".\data\raw\GSE229711_RAW.tar" `
  ".\data\derived\GSE230809_discovery_raw_data_ledger.csv" `
  ".\data\derived\annotation\GSE229711_RAW_cell_annotation.csv.gz" `
  ".\config\program_modules.json" `
  --qc-cells ".\data\derived\qc\GSE229711_RAW_cell_qc.tsv.gz" `
  --output-dir ".\data\derived\module_scores" `
  --min-retained-cells 30
```

Repeat with `GSE230808_RAW.tar`, its matching annotation/QC files, and the same
output directory. Thresholds 20 and 50 are sensitivity analyses, not alternate
primary analyses.

## Scoring and audit rules

For each module gene, counts are summed across all matching feature rows after
uppercase symbol normalization. Duplicate feature rows are retained in the
mapping audit and summed, rather than silently dropped. A module score is the
mean over mapped genes of `log1p(1e6 * pseudobulk_gene_count /
total_UMI_included_cells)`, and a module is scoreable only when at least 80% of
its configured genes are measured and included-cell UMI is positive.

The scorer validates Matrix Market dimensions, coordinate bounds, non-negative
counts, and declared versus observed coordinate-record counts. It writes module
scores, gene-level pseudobulk counts, mapping audits, a library ledger, and a
parameter record containing the QC path and retained-cell threshold.

## Interpretation boundary

The discovery comparison is an exploratory advanced-degeneration-associated
contrast because age and disease state are fully confounded and the comparison
group has only three donors per compartment. Do not report p-values from these
scores as confirmatory evidence, and do not treat AF/NP child-series results as
independent replication.
