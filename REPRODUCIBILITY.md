# Reproducibility Guide

This document distinguishes two supported activities:

1. Rebuild the released figures and tables from the included derived result
   records.
2. Recompute analysis stages from externally downloaded public source data.

The repository deliberately does not claim a one-command, raw-data-to-paper
pipeline. Cohort-specific acquisition, audit, scoring, and summary steps are
kept explicit so that their eligibility and interpretation boundaries remain
inspectable.

## Software

Create a Python environment with Python 3.10 or newer and install the exact
versions declared in tools/python/requirements-lock.txt:

~~~powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -r .\tools\python\requirements-lock.txt
~~~

For the optional S7-S9 random-effects analyses, use R 4.4.1 and run the
non-mutating version check:

~~~powershell
Rscript .\tools\r\check_required_packages.R
~~~

The required R package versions are recorded in tools/r/requirements.tsv.

## Restore public inputs

The external-source manifest contains URLs, expected byte counts, and SHA-256
digests for every raw GEO/Ensembl asset. It supplies a reproducible route to
the files excluded from the public code release:

~~~powershell
# Frozen default-analysis inputs.
.\.venv\Scripts\python.exe .\scripts\fetch_public_data.py --group default

# Isolated exploratory and S8/S9 sensitivity inputs.
.\.venv\Scripts\python.exe .\scripts\fetch_public_data.py --group exploratory --group s8 --group s9

# Every declared asset, including candidate-only and Ensembl audit inputs.
.\.venv\Scripts\python.exe .\scripts\fetch_public_data.py --all

# Offline checksum verification only.
.\.venv\Scripts\python.exe .\scripts\fetch_public_data.py --all --verify-only
~~~

See data/DATA_ACCESS.md for upstream data terms and scope labels.

## Frozen display regeneration

The following commands rebuild the release presentation from the included
derived records. They should reproduce figures and tables without changing the
scientific result contract.

~~~powershell
.\.venv\Scripts\python.exe .\scripts\make_current_summary_deliverables.py --summary-dir .\data\derived\donor_module_effect_summary --table-dir .\results\tables --figure-dir .\results\figures
.\.venv\Scripts\python.exe .\scripts\make_gse251686_supplement.py --package-dir .\data\derived\GSE251686_exploratory_scores --table-dir .\results\supplementary_tables --figure-dir .\results\supplementary_figures
.\.venv\Scripts\python.exe .\scripts\make_results_viewer.py
~~~

The support-deliverable generator is intentionally strict about its source
contracts:

~~~powershell
.\.venv\Scripts\python.exe .\scripts\make_submission_support_deliverables.py --default-summary-dir .\data\derived\donor_module_effect_summary --gse251686-package-dir .\data\derived\GSE251686_exploratory_scores --program-ledger .\data\derived\program_module_ledger.csv --discovery-sensitivity-dir .\data\derived\discovery_retained_cell_sensitivity --results-root .\results
~~~

## Cohort-specific recomputation map

| Stage | Primary scripts | Main records |
| --- | --- | --- |
| Download and verification | scripts/fetch_public_data.py | data/public_data_manifest.tsv |
| 10x archive structure and QC | audit_10x_tar.py, audit_nested_10x_tar.py, qc_10x_archives.py, qc_nested_10x_archives.py | data/derived/qc/ and audit ledgers |
| Discovery pseudobulk scoring | score_module_pseudobulk.py, score_nested_10x_modules.py | data/derived/module_scores_recomputed/ |
| GSE165722 score-level analysis | audit_gse165722.py, score_gse165722_processed_modules.py | data/derived/module_scores_external/GSE165722/ |
| External count and processed-expression scoring | audit_gse153066_count_matrix.py, audit_gse244889_assets.py, score_dense_matrix_modules.py, score_external_human_np_bulk_candidates.py | data/derived/module_scores_external/ |
| Default descriptive summary | summarize_donor_module_effects.py | data/derived/donor_module_effect_summary/ |
| Retained-cell sensitivity | summarize_discovery_retained_cell_sensitivity.py | data/derived/discovery_retained_cell_sensitivity/ |
| Isolated GSE251686 supplement | score_gse251686_exploratory.py, summarize_gse251686_exploratory.py, make_gse251686_supplement.py | data/derived/GSE251686_exploratory_scores/ |
| S7-S9 exploratory syntheses | run_np_exploratory_meta_analysis.R, run_np_post_hoc_external_expansion_meta_analysis.R | data/derived/np_*_meta_analysis/ |
| Manuscript sources | build_formal_manuscripts.py | manuscript/formal_submission/ |

The exact input and output hashes for frozen result products are retained in
their respective run-artifact and manifest records. The public release
snapshot normalizes text line endings and rewrites machine-root prefixes in
those records, leaving analytical values unchanged. RELEASE_FILE_MANIFEST.tsv
is the authoritative checksum inventory for the published snapshot; the
per-run hashes preserve the provenance of the original local execution
artifacts.

## Release verification

Build a clean candidate outside the project root, then audit it:

~~~powershell
.\.venv\Scripts\python.exe .\scripts\build_public_release_snapshot.py --output ..\ivdd_cross_cohort_reproducibility_release
.\.venv\Scripts\python.exe .\scripts\audit_public_release.py ..\ivdd_cross_cohort_reproducibility_release
~~~

After author-approved LICENSE, CITATION.cff, and .zenodo.json have been added,
repeat the audit with the metadata requirement:

~~~powershell
.\.venv\Scripts\python.exe .\scripts\audit_public_release.py ..\ivdd_cross_cohort_reproducibility_release --require-release-metadata
~~~
