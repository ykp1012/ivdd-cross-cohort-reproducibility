# Citation Audit

Audit date: 2026-08-14.

## Verification procedure

Journal metadata were checked against PubMed/Europe PMC records using PMID or
DOI. GEO accession title, series relationships, sample counts, linked PMIDs,
and the `Citation missing` status were checked against the local GEO SOFT
archives and NCBI GEO accession pages. The bibliography was then checked with
the citation-management validator (`27` entries; `0` errors, `0` warnings,
`0` duplicates). Seventeen DOI-bearing entries were also checked against
Crossref metadata. Seven DOI resolvers rejected an HTTP `HEAD` request from the
generic validator, but all seven returned matching Crossref records via `GET`;
they are not citation errors.

## Numbered sources

| No. | Entry | Verification result and manuscript use |
|---:|---|---|
| 1 | Barrett et al., 2013; PMID 23193258; DOI 10.1093/nar/gks1193 | PubMed record verified. Supports the NCBI GEO repository description. |
| 2 | Squair et al., 2021; PMID 34584091; DOI 10.1038/s41467-021-25960-2 | PubMed/Europe PMC and Crossref metadata verified. Supports the warning that cell-level DE can ignore biological-replicate variation. It is not used to justify a specific IVDD effect. |
| 3 | Swahn et al., 2024; PMID 38403470; DOI 10.1002/advs.202309032 | PubMed/Europe PMC full record verified. Source publication for the GSE230809 parent project. |
| 4-6 | GSE230809, GSE229711, and GSE230808 | NCBI GEO SOFT and accession pages verified. The two child series are SubSeries of the one GSE230809 parent project; they are not independent validation cohorts. |
| 7 | Chen et al., 2024; PMID 38167807; DOI 10.1038/s41467-023-44313-9 | PubMed/Europe PMC full record and PMC data-availability statement verified. Linked publication for GSE244889. |
| 8 | GSE244889 | NCBI GEO SOFT/page verified. The record lists PMID 38167807 and later related publications; only the scRNA series records used by this project are cited here. |
| 9 | GSE153066 | NCBI GEO SOFT/page verified. The GEO page explicitly says `Citation missing`; no article is attributed. It is cited only as a database record. |
| 10 | Tu et al., 2022; PMID 34825784; DOI 10.1002/advs.202103631 | PubMed/Europe PMC full record and PMC data-availability statement verified. Source publication for GSE165722. |
| 11 | GSE165722 | NCBI GEO SOFT/page verified. GEO labels the supplied integer-like matrices as normalized counts; this supports the score-level boundary, not raw-count inference. |
| 12 | Wang et al., 2023; PMID 37216089; DOI 10.1016/j.isci.2023.106692 | PubMed/Europe PMC full record verified. Background/literature-boundary citation for prior human AF/NP single-cell work. |
| 13 | Sun et al., 2025; PMID 39828732; DOI 10.1038/s41413-024-00372-2 | PubMed/Europe PMC full record verified. Background/literature-boundary citation for prior multi-dataset NP/fibrocyte work. |
| 14 | Crowell et al., 2020; PMID 33257685; DOI 10.1038/s41467-020-19894-4 | PubMed/Europe PMC full record verified. Supports multi-sample/multi-condition single-cell analysis principles; no `muscat` result is reported here. |
| 15 | Niu et al., 2025; PMID 40453974; DOI 10.2147/JIR.S519218 | PubMed/Europe PMC full record verified. Background/literature-boundary citation; no novelty or therapeutic claim is adopted. |
| 16 | Jia et al., 2024; PMID 39516278; DOI 10.1038/s41598-024-78675-x | PubMed/Europe PMC full record and PMC data-availability statement verified. Source publication for GSE251686. |
| 17 | GSE251686 | NCBI GEO SOFT/page verified. Six records are listed; `GSM7986002` is excluded after this project's stream-integrity audit. The remaining records underwent a separate 2-versus-3 exploratory score analysis but were deliberately omitted from the default 20-effect summary and are not presented as confirmatory validation or replication. |
| 18 | Kanehisa et al., 2023; PMID 36300620; DOI 10.1093/nar/gkac963 | PubMed/Europe PMC full record verified. Supports KEGG as a pathway knowledgebase used in module definition. |
| 19 | Jassal et al., 2020; PMID 31691815; DOI 10.1093/nar/gkz1031 | PubMed/Europe PMC full record verified. Supports Reactome as a pathway knowledgebase used in module definition. |
| 20 | Welch, 1947; DOI 10.1093/biomet/34.1-2.28 | DOI resolves to Project Euclid/Crossref metadata; supports the unequal-variance interval convention. |
| 21 | Efron, 1979; DOI 10.1214/aos/1176344552 | DOI resolves through Crossref to the Project Euclid record; supports bootstrap resampling as a descriptive interval method. |
| 22 | GSE186542 | Local GEO SOFT verified. It has no linked PubMed ID in the GEO metadata and is cited only as an accession; it supplies the audited 3-versus-3 score-level S8/S9 input, not a raw-count or confirmatory result. |
| 23 | GSE167931 | Local GEO SOFT verified. It lists PMIDs 35304463 and 35340126; its FPKM representation is used once in S8 and its paired TPM matrix is not a second cohort. |
| 24 | Li et al., 2022; PMID 35304463; DOI 10.1038/s41467-022-28990-6 | PubMed and Crossref metadata verified. One GEO-linked publication for GSE167931; used only to establish accession provenance. |
| 25 | Li et al., 2022; PMID 35340126; DOI 10.1002/ctm2.765 | PubMed and Crossref metadata verified. The second GEO-linked publication for GSE167931; used only to establish accession provenance. |
| 26 | GSE245147 | Local GEO SOFT verified. It lists PMID 38488012 and supplies the native 3-versus-3 RPKM subset used only as the S9 source-family replacement for GSE167931. |
| 27 | Zhang et al., 2024; PMID 38488012; DOI 10.1172/JCI165140 | PubMed and Crossref metadata verified. Linked publication for GSE245147; used only to establish accession provenance. |

## Claims deliberately left unsupported

- No source supports calling any current cohort a confirmed replication or
  external validation. The manuscript uses `directional support`,
  `exploratory`, and `descriptive` language only.
- No source supports an age-independent GSE230809 disease effect; age and
  recorded disease state are fully confounded in that parent project.
- No source supports a causal mechanism, prognostic biomarker, therapeutic
  target, AF-to-NP migration claim, ROC classifier, WGCNA hub-gene claim,
  CellChat result, docking result, or clinical benefit/harm statement for this
  analysis. Those claims are excluded from the manuscript.
- The four study-defined modules are expression summaries. KEGG/Reactome and
  the cited IVDD papers document pathway/source context, but do not validate
  the exact locked gene lists as clinical signatures.
- GSE153066 has no verified linked publication in its GEO record as of the
  audit date; a manuscript sentence attributing it to a named paper would be
  unsupported and was not written.

## Citation-style note

The manuscript uses Vancouver-style numeric citations in order of first
appearance. The bibliography is maintained as BibTeX for reproducibility; a
target journal's final reference-style conversion remains a submission-stage
task.

## Audit outputs

- `manuscript/references.bib`: machine-readable source metadata in first-use
  numeric order.
- `manuscript/citation_validation.json` and
  `manuscript/citation_validation_final.json`: validator outputs for the final
  bibliography (27 entries; 0 errors, 0 warnings, 0 duplicates). DOI verification
  is documented above, including the Crossref `GET` fallback for seven resolvers
  that rejected the validator's HTTP `HEAD` request.
- `manuscript/04_manuscript_draft.md`: all external background, cohort
  provenance, pathway-resource, and named statistical-method statements have
  an in-text numeric citation; results generated by this project are cited to
  its tables, figures, and audit artifacts rather than to external papers.
