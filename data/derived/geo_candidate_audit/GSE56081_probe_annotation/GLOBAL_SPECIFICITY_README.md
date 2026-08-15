# GSE56081 Global Probe Specificity Audit

This is a candidate-only extension audit. It does not modify the frozen default IVDD result.

GPL15314 candidate probes were searched exactly (A/C/G/T, both orientations) against every
Ensembl release-113 human cDNA and ncRNA record and the GRCh38 primary assembly. The
dependency-free matcher retains overlapping hits. The ledger records all transcript gene IDs
and genomic loci; probes with cross-gene transcript hits, multiple genomic loci, or an
unexpected genomic overlap are excluded from any exploratory score.

Reference release: Ensembl 113 / GRCh38; cDNA records=410909; primary-assembly records=194.
Candidate probes searched: 82; failed specificity probes: 28.

## Module Gate

| Module | Candidate probes | Globally specific probes | Specific genes | Fraction | 0.80 gate |
|---|---:|---:|---:|---:|---|
| ecm_collagen_remodeling | 26 | 20 | 17 | 0.708333 | fail |
| inflammatory_nfkb | 18 | 12 | 12 | 0.571429 | fail |
| hypoxia_oxidative_stress | 23 | 12 | 9 | 0.500000 | fail |
| disc_matrix_homeostasis | 15 | 10 | 9 | 0.600000 | fail |

Decision: `candidate_only_blocked_global_probe_specificity_not_resolved`.
A pass here is still sequence-anchored evidence rather than manufacturer-certified
annotation; the 2011 Arraystar design/transcript-version and probe summarization boundaries
remain documented in the first-stage audit.
