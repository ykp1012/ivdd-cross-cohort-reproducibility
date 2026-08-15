# Journal-Neutral Submission Readiness Checklist

## Purpose and Current Boundary

This checklist prepares the cohort-aware IVDD manuscript for a journal
submission without widening its scientific claims. The current package is a
reproducible local project package, not a public code archive. It must not be
described as publicly archived until the exact version used for submission has a
public repository URL and an archival DOI.

## Current Package Status

| Area | Current status | Evidence or action required before submission |
|---|---|---|
| Default result contract | Complete locally | The default summary has 20 effects, 55/55 exact score-to-ledger matches, and all cohorts marked confirmatory_eligible=false; see Supplementary Tables S2 and S6. |
| GSE251686 boundary | Complete locally | GSM7986002 remains permanently excluded and the remaining 2-versus-3 result is isolated from the default summary; see Supplementary Table S1 and Supplementary Figure S1. |
| Main figures and tables | Complete locally | Check the target journal's image size, resolution, color, file-format, and table-template requirements. |
| Supplementary package | Complete locally | Submit Tables S1-S9 and Figures S1-S5 in the journal's requested formats; retain the graphical abstract separately if required. |
| Citations | Complete locally | citation_audit.md records 27 validated bibliography entries with no reported errors, warnings, or duplicates; Crossref `GET` metadata checks resolve the seven DOI entries whose publishers reject the generic validator's HTTP `HEAD` request. Recheck after any reference edit or journal-style conversion. |
| Public code and data-release archive | Not complete | Create a frozen public repository release and archive it with a persistent DOI before submission. Insert the final URL and DOI in the manuscript. |
| Author, funding, conflict, and ethics declarations | Not complete | Obtain and insert author-approved, institution-specific statements. Do not infer these details from the analysis package. |
| Target-journal formatting | Not complete | Select a journal, then apply its author guide, reporting requirements, word limits, citation style, title-page fields, and submission forms. |

## Scientific Claim Check

Before submission, confirm that every title, abstract, figure legend, response
letter, and cover letter preserves the following boundaries. Cells are nested
observations, not independent replicates. The default result refers only to the
four cohorts included in the 20-effect summary. Positive hypoxia/oxidative-
stress point estimates are descriptive because all four Welch intervals include
zero. GSE251686 is an isolated exploratory sensitivity result and must not be
pooled, counted in default sign alignment, or presented as validation or
replication. No text should claim a universal IVDD program, mechanism,
biomarker, treatment target, causal effect, or confirmatory result. The default
analysis must not be described as having a formal meta-analysis or p-values;
S7-S9 may be described only as non-confirmatory exploratory syntheses with
transparent HKSJ/BH p-values.

## Release and DOI Plan

1. Freeze the exact source tree used for submission, including scripts, config,
   docs, manuscript, tools/python/requirements-lock.txt, the derived tables
   needed to reproduce the displayed results, and the manifests referenced in
   Supplementary Table S6.

2. Review the release contents for public redistribution constraints. Public GEO
   accessions and retrieval instructions can be cited without redistributing
   source archives where a repository policy, file size, or license makes that
   inappropriate. Do not release identifiers or materials that should remain
   restricted.

3. Create a versioned public repository release, tag the exact commit, and
   archive that release with a DOI-providing service. Record the repository URL,
   release tag, archival DOI, release date, and the SHA-256 hashes from
   Supplementary Table S6 in the release notes.

4. Replace the manuscript placeholders in Data Availability and Code
   Availability only after the archive is live. Re-run the reproducibility and
   artifact checks against the frozen release rather than a subsequently edited
   working directory.

## Final File Checks

Use the project-local Python environment from the project root. The first
command checks the generators for syntax. The next three regenerate the main
tables/figures, isolated GSE251686 supplement, and submission-support artifacts
after any legitimate upstream result change. Regeneration must not be used to
change the analysis boundary.

~~~powershell
$ivddPython = '.\tools\python\venv\Scripts\python.exe'
& $ivddPython -m py_compile .\scripts\make_current_summary_deliverables.py .\scripts\make_gse251686_supplement.py .\scripts\make_submission_support_deliverables.py

& $ivddPython .\scripts\make_current_summary_deliverables.py --summary-dir .\data\derived\donor_module_effect_summary --table-dir .\results\tables --figure-dir .\results\figures

& $ivddPython .\scripts\make_gse251686_supplement.py --package-dir .\data\derived\GSE251686_exploratory_scores --table-dir .\results\supplementary_tables --figure-dir .\results\supplementary_figures

& $ivddPython .\scripts\make_submission_support_deliverables.py --default-summary-dir .\data\derived\donor_module_effect_summary --gse251686-package-dir .\data\derived\GSE251686_exploratory_scores --program-ledger .\data\derived\program_module_ledger.csv --discovery-sensitivity-dir .\data\derived\discovery_retained_cell_sensitivity --results-root .\results
~~~

After the generators finish, verify results/submission_support_manifest.json
against the generated files, excluding the manifest itself because it would be
self-referential. Confirm that the default effect count remains 20, the exact
identity contract remains 55/55, gse251686_default_summary_inclusion remains
false, and GSM7986002 remains the permanent exclusion. Verify the bibliography
after all author edits and re-inspect the PDF/PNG visual outputs at the target
journal's required dimensions.

## Submission Handoff

Before uploading, obtain author approval for the final manuscript and all
declarations, select the target journal, convert the citation and figure/table
formats, and submit the graphical abstract only if requested. Retain the
current quality-audit, citation-audit, and reproducibility-contract files with
the frozen release so that any reviewer question can be traced to an exact
input, script, and generated artifact.
