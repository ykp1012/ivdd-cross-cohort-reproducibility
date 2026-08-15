# Candidate Human NP Cohort Audit

This directory records a reproducible, accession-level audit of candidate
human nucleus pulposus (NP) cohorts.  It is an extension screen only.  None of
these files changes the frozen default IVDD results under
`data/derived/donor_module_effect_summary/`.

## Re-run

Run from the project root with the project-local environment:

```powershell
& ".\tools\python\venv\Scripts\python.exe" `
  ".\scripts\audit_geo_candidate_cohorts.py" `
  --project-root "." `
  --output-dir ".\data\derived\geo_candidate_audit"
```

The audit deliberately refuses to overwrite an existing generated artifact.
This prevents a later run from silently changing its evidence record.  To
refresh it after a documented input change, archive the prior generated audit
files first and run again.

## Generated Evidence

- `candidate_input_artifact_hashes.csv`: input paths, sizes, and SHA-256 hashes.
- `candidate_sample_matrix_mapping.csv`: every ledger sample mapped to one
  processed-matrix column and its predeclared comparison group.
- `candidate_matrix_integrity_audit.csv`: processed-matrix dimensions, feature
  identifiers, numeric checks, and raw-count-model boundary.
- `candidate_module_mapping_audit.csv`: coverage of the locked program genes.
- `candidate_independence_ledger.csv`: BioProject, BioSample, parent-series,
  and stated identical-RNA relations.
- `candidate_cohort_assessment.csv`: the operational inclusion recommendation.

## Decision Boundary

`GSE186542` is conditionally suitable as an additional small external NP
score-level comparison: its early-stage (Pfirrmann I-III, n=3) and advanced
stage (IV-V, n=3) groups map exactly to the six count-labelled columns, and all
four locked modules map completely.  Its GEO processing text simultaneously
mentions FPKM/TMM and calls the supplementary file raw counts.  It must not be
used for a raw-count negative-binomial model until that provenance conflict is
resolved from primary count generation or SRA reprocessing.

`GSE167931` is conditionally suitable as an external normalized-expression NP
score cohort: normal n=4 and degenerated n=5 map exactly to both TPM and FPKM
files after the documented hyphen-to-underscore title normalization, and all
four locked modules map completely.  Its values are normalized expression, so
any extension must predeclare a log2(TPM + 1) or log2(FPKM + 1) score scale and
must not use DESeq2/edgeR raw-count inference.

`GSE56081` has a valid 5-versus-5 grade I versus IV/V human NP microarray
comparison. The initial accession-level audit recorded that the local GPL15314
table uses custom probe IDs and lacks an immediately usable gene identifier
map. A later fixed Ensembl release-113/GRCh38 exact-sequence audit searched 82
candidate probes in both orientations against the whole transcriptome and
primary assembly; 28 probes failed the global-specificity rule. The resulting
globally specific-gene fractions are 0.7083, 0.5714, 0.5000, and 0.6000 across
the four locked modules, so none meets the 0.80 gate. `GSE56081` is therefore
blocked as a candidate-only microarray extension and excluded from S7-S9; see
`GSE56081_probe_annotation/GLOBAL_SPECIFICITY_README.md` and
`GSE56081_probe_annotation/global_specificity_manifest.json`. It is a GSE67567
subseries and reports identical-RNA alternatives, so its biological sample set
must be counted only once.

The distinct BioProjects support accession-level independence from the current
default IVDD cohort set.  They do not prove donor-level independence: GEO does
not expose patient identifiers for GSE186542 or GSE167931, and publication-level
donor-overlap checks remain a separate requirement.
