# IVDD Cross-Cohort Reproducibility

Code, configuration, derived result records, and manuscript materials for a
cohort-aware descriptive analysis of pre-specified nucleus pulposus (NP)
transcriptional programs across public human intervertebral disc degeneration
(IVDD) datasets.

The unit of observation is the donor or clearly labeled presumed
donor/sample/library key, not the individual cell. The frozen default
analysis is descriptive: it reports effect sizes, intervals, donor bootstrap
stability, and sign alignment. It does not establish confirmatory,
mechanistic, causal, biomarker, or therapeutic conclusions.

## Release status

The first public code release is v1.0.0. The permanent Zenodo DOI will be
added here after the versioned GitHub release has been archived. Creator and
license metadata are recorded in CITATION.cff, .zenodo.json, LICENSE, and
RELEASE_METADATA.md.

## Research question

Which pre-specified NP programs have directionally aligned or discordant
donor-level associations with advanced degeneration or severity across
independently processed human datasets, and how are those patterns limited by
age, cohort, processing scale, or cell composition?

## Interpretation boundary

- Cells are not independent biological replicates.
- AF and NP analyses remain separate unless paired donor metadata support an
  interaction model.
- The AF/NP discovery project and external support cohorts remain separate;
  no current cohort is presented as confirmatory validation.
- The default 20-effect summary contains no p-values, formal meta-analysis, or
  replication adjudication.
- This observational public-data study does not support claims about migration,
  ligand-receptor mechanisms, therapeutic targets, causality, or molecular
  docking.

## Key finding

Across the four frozen default NP cohort roles, the hypoxia/oxidative-stress
program had directionally aligned point estimates, but all four Welch 95%
intervals crossed zero. The ECM, inflammatory/NF-kB, and
neurovascular/remodeling programs showed directionally discordant estimates.
These results describe cross-cohort heterogeneity and should not be treated as
confirmatory evidence.

## Quick start

Use Python 3.10 or later. The public release does not include an environment
directory; create one locally and install the version-pinned Python
dependencies:

~~~powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -r .\tools\python\requirements-lock.txt
~~~

The optional S7-S9 exploratory random-effects syntheses require R 4.4.1,
metafor 4.8.0, digest 0.6.37, and jsonlite 2.0.0. Verify an R environment with:

~~~powershell
Rscript .\tools\r\check_required_packages.R
~~~

## Data access

Original GEO and Ensembl files are not redistributed in this repository or
its Zenodo archive. The data are large, third-party assets and must be
retrieved from their source records. The checked manifest at
data/public_data_manifest.tsv records 36 upstream URLs, byte counts, and
SHA-256 digests. Full instructions and data-use boundaries are in
data/DATA_ACCESS.md.

~~~powershell
# Download and verify the frozen default-analysis inputs.
.\.venv\Scripts\python.exe .\scripts\fetch_public_data.py --group default

# Download every declared source asset.
.\.venv\Scripts\python.exe .\scripts\fetch_public_data.py --all

# Verify previously downloaded files without network access.
.\.venv\Scripts\python.exe .\scripts\fetch_public_data.py --all --verify-only
~~~

## Rebuild released displays

The release snapshot includes the derived result records needed to rebuild the
publication displays. These commands regenerate presentation artifacts; they
do not recompute cohort scores from raw archives.

~~~powershell
.\.venv\Scripts\python.exe .\scripts\make_current_summary_deliverables.py --summary-dir .\data\derived\donor_module_effect_summary --table-dir .\results\tables --figure-dir .\results\figures
.\.venv\Scripts\python.exe .\scripts\make_results_viewer.py
~~~

For the isolated GSE251686 exploratory supplement:

~~~powershell
.\.venv\Scripts\python.exe .\scripts\make_gse251686_supplement.py --package-dir .\data\derived\GSE251686_exploratory_scores --table-dir .\results\supplementary_tables --figure-dir .\results\supplementary_figures
~~~

Detailed environment, data-provenance, and analysis records are in docs/.
REPRODUCIBILITY.md maps the cohort-specific scoring and synthesis scripts to
their intended inputs and output classes.

## Result products

- data/derived/donor_module_effect_summary/: the only authoritative frozen
  default cross-cohort result directory.
- data/derived/GSE251686_exploratory_scores/: an intentionally separate,
  non-confirmatory mild n=2 versus severe n=3 exploratory package.
- data/derived/np_exploratory_meta_analysis/: S7, an isolated four-cohort
  exploratory random-effects synthesis.
- data/derived/np_post_hoc_external_expansion_meta_analysis/: S8, an isolated
  post hoc external-expansion synthesis.
- data/derived/np_source_family_replacement_meta_analysis/: S9, a
  source-family replacement sensitivity synthesis.
- results/: figures, tables, supplementary figures, and supplementary tables.
- ivdd_results_viewer.html: a local viewer for rendered figures and key tables.

S7-S9 use SMDH random-effects synthesis with REML and Knapp-Hartung intervals.
Their HKSJ and four-module BH values are reported only for transparency and do
not support confirmation, replication, biomarker, mechanism, causal, or
therapeutic claims. Patient-level independence is not established for the
public accessions used in those analyses.

## Analysis entry points

- scripts/fetch_public_data.py: downloads and verifies every declared external
  source asset.
- scripts/score_module_pseudobulk.py, score_gse165722_processed_modules.py,
  score_dense_matrix_modules.py, and score_external_human_np_bulk_candidates.py:
  cohort-specific score generation.
- scripts/summarize_donor_module_effects.py: default and sensitivity
  donor/library-key effect summaries.
- scripts/make_current_summary_deliverables.py,
  make_gse251686_supplement.py, and make_submission_support_deliverables.py:
  figures, tables, and supplementary deliverables.
- scripts/run_np_exploratory_meta_analysis.R and
  run_np_post_hoc_external_expansion_meta_analysis.R: isolated exploratory
  S7-S9 syntheses.
- scripts/build_formal_manuscripts.py: manuscript source generation.
- scripts/export_formal_pdfs.ps1: optional Windows/WPS PDF export only; it is
  not a core analysis dependency.

## Repository layout

- config/: locked program modules, marker panels, and contrast specifications.
- data/: public-data acquisition manifest, access notes, and generated records.
- docs/: protocol, provenance, audit, and environment documentation.
- manuscript/: manuscript source and formal submission materials.
- results/: generated figures and tables.
- scripts/: acquisition, audit, scoring, summary, and figure-generation code.
- tools/: Python dependency lock and R version specification.

## Public-release boundary

The public repository and Zenodo archive include original code, configuration,
documentation, manuscript materials, generated result records, and rendered
outputs. They exclude raw GEO/Ensembl data, local virtual environments,
machine-resource snapshots, caches, logs, and temporary files. The release
snapshot generator and audit are located at:

- scripts/build_public_release_snapshot.py
- scripts/audit_public_release.py

## Citation and license

Use the version DOI assigned by Zenodo when citing a specific release, and the
Zenodo concept DOI when citing the project across versions. CITATION.cff and
.zenodo.json provide the software citation metadata. The repository's original
code and documentation are released under the MIT License; third-party GEO and
Ensembl data remain subject to their respective source terms.
