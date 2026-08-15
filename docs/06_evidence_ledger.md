# Evidence Ledger

Audit date: 2026-08-13.

This ledger separates study evidence from implementation evidence. It does not
claim that a paper's code was reproduced unless a checked public repository is
recorded.

| Item | Role in this project | Paper evidence | Implementation evidence | Evidence and boundary |
|---|---|---:|---:|---|
| Swahn et al., 2024, PMID 38403470, DOI 10.1002/advs.202309032 | Source study for GSE230809 discovery project | B | F | Human AF/NP single project. Its child series GSE229711 and GSE230808 are one discovery cohort, never discovery plus validation. Existing AF/NP and fibrosis findings are background, not a novelty claim. |
| Tu et al., 2021/2022, PMID 34825784, DOI 10.1002/advs.202103631 | Source study for GSE165722 | B | F | Table 1 reports NP grades II-V, two donors per grade. GEO SOFT labels the ordered samples I-IV; this conflict is retained. GEO calls supplied values normalized counts, so this cohort is score-level direction support only. |
| Wang et al., 2023, PMID 37216089 | Literature boundary | B | T | Existing human AF/NP atlas and cell-state analyses preclude novelty claims based on an AF/NP atlas, trajectory, or CellChat alone. |
| Niu et al., 2025, PMID 40453974 | Literature boundary | C | T | Existing multi-cohort NP pooled-cell integration; do not claim the first multi-cohort human IVDD integration. |
| Sun et al., 2025, PMID 39828732 | Literature boundary | C | T | Existing multi-dataset NP integration and fibrosis/fibrocyte framing; do not repackage FibroNP/fibrosis as novelty. |
| Liang et al., 2025, PMID 40515444 | Literature boundary | B | T | Existing animal evidence relevant to AF-derived contributions; observational human data cannot support AF-to-NP migration claims. |
| Zhang et al., 2025, PMID 40570207 | Literature boundary | B | R | Do not recast LCN2/MDSC or IL1B-style staging hypotheses as new without independent evidence. |

## Claim rules

- The defensible contribution is a pre-specified, independently processed,
  donor-level NP reproducibility audit with heterogeneity reported explicitly.
- No "first" claim is permitted; absence of an exact matched publication is
  not proof of priority.
- AF results are discovery-stage only. No independent AF severity support
  cohort has yet passed audit.
- Evidence labels: A = original multi-centre/spatiotemporal with independent
  validation; B = human/animal omics plus mechanism or intervention; C =
  primarily public-data reanalysis. R = checked method/repository; F = full
  text/methods checked but repository unverified; T = title/abstract-level
  extraction; U = not checked.
