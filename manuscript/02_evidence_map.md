# Evidence Map

## Topic
Cohort-aware direction and heterogeneity of locked nucleus pulposus transcriptional programs across human intervertebral disc degeneration datasets

## Foundational studies
- Tu et al. characterized human NP degeneration with GSE165722 (PMID 34825784). The dataset supports a severity-ordered external cohort, but the available GEO supplementary matrices are described as normalized counts and have a GEO-versus-source grade-label conflict.
- Swahn et al. reported a human AF/NP study corresponding to the GSE230809 parent project (PMID 38403470). It is the sole AF/NP discovery project in this analysis, not an independent support pair.

## Recent studies
- Niu et al. (PMID 40453974) and Sun et al. (PMID 39828732) already performed multi-cohort human NP integrations, chiefly with pooled-cell workflows. This study does not claim the first multi-cohort integration.
- GSE186542 has no linked PubMed record in its GEO SOFT metadata and is cited as an accession only. Li et al. (PMIDs 35304463 and 35340126) are the two linked publications for GSE167931, and Zhang et al. (PMID 38488012) is the linked publication for GSE245147. These sources establish accession provenance only; they do not validate the present locked-module effects or any mechanistic interpretation.

## Conflicting or limiting evidence
- In GSE230809, age and advanced degeneration are completely confounded. Its result is an advanced-degeneration-associated discovery contrast, not an age-independent disease effect.
- Independent public direction support is currently NP-only. AF findings must remain exploratory.
- GSE165722's supplementary matrices are count-like but GEO labels them normalized; it is excluded from raw-count negative-binomial inference.
- GSE186542 has a conflicting GEO description of count, FPKM, and TMM provenance, and GSE167931/GSE245147 expose no patient identifiers. GSE245147 therefore replaces, rather than supplements, GSE167931 in the S9 source-family sensitivity analysis. GSE56081 remains outside the extension because its fixed Ensembl-113/GRCh38 audit leaves every locked module below the 0.80 global probe-specificity gate.

## Candidate gap statement
- Whether locked NP molecular programs show descriptively aligned or discordant directions after each public cohort is processed independently and the donor or presumed sample key, rather than thousands of cells, is treated as the unit of observation remains a narrow, testable evidence-audit question. This is a gap statement, not a priority claim.

## Keywords / search strings
- "intervertebral disc degeneration" AND "nucleus pulposus" AND "single-cell"
- "intervertebral disc degeneration" AND "cross-cohort" AND "pseudobulk"
- "GSE230809" OR "GSE165722" OR "GSE244889" OR "GSE251686" OR "GSE186542" OR "GSE167931" OR "GSE245147"
