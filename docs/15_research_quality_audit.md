# Research Quality Audit

Audit date: 2026-08-14. Scope: current IVDD cross-cohort module-score
pipeline, cohort roles, provenance, inference units, and interpretation.

## Executive decision

The pipeline is technically auditable for the analyses that have passed their
data gates, but the project does **not yet support a confirmatory biological
claim or a causal/mechanistic manuscript**. The defensible product is an
evidence-limited reproducibility report: pre-specified NP program scores are
compared at the presumed donor/sample level, with effect sizes, intervals,
leave-one-key-out stability, mapping audits, and explicit cohort limitations.

No cells were used as independent biological replicates. The frozen default
20-effect descriptive summary contains no p-values or formal meta-analysis.
Separately packaged S7, S8, and S9 are non-confirmatory random-effects
syntheses; their HKSJ and four-module BH p-values are reported for transparency
only and cannot support confirmation, replication, biomarker, mechanistic,
causal, or therapeutic claims. No negative-binomial model, ROC classifier,
WGCNA, hub-gene ranking, CellChat, molecular docking, migration claim,
treatment claim, or causal claim is supported by any current output.

## Five-dimension audit

| Dimension | Status | Main finding |
|---|---|---|
| Terminology | PASS with limits | Gene symbols and module names are explicit; use “presumed donor/sample key” where patient identity is not exposed. Use association/direction language only. |
| Statistical model | PASS with non-confirmatory extensions | Mean target-minus-comparison score, Welch interval, donor bootstrap, and LODO operate at the donor/library key in the default summary. S7-S9 add SMDH random-effects REML synthesis with Knapp-Hartung intervals and transparent HKSJ/BH p-values; they are not confirmatory tests. |
| Composite/module definition | PASS with limits | Four gene lists and the 80% mapping rule are locked before external scoring. The score is a transformed within-sample program summary, not a clinical index or absolute cross-platform abundance. |
| Provenance and reproducibility | PASS for audited inputs | Raw archives are retained unchanged; matrix dimensions, stream integrity, barcode/sample mappings, score parameters, and input hashes are recorded. |
| Data scope and generalizability | WARN/FAIL for strong claims | Small groups and `k = 4`/`k = 6` synthesis units, presumed identities, age/source/grade confounding, platform and processing-scale differences, unverified patient-level independence, source-family replacement, and one malformed excluded library limit inference and external validity. |

## Cohort-level decisions

| Cohort | Current role | Safe statement | Prohibited statement |
|---|---|---|---|
| GSE230809 parent (GSE229711 + GSE230808) | Exploratory discovery | Advanced-degeneration-associated donor-level effect estimates in AF/NP | Age-independent degeneration effect, independent validation, or confirmatory significance |
| GSE165722 | NP score-level direction support | Severe-minus-mild direction among 4 versus 4 presumed Sample/GSM keys | Raw-count differential expression, normalized values called UMIs, or confirmed replication |
| GSE153066 | NP external count-level support only | Degenerated-minus-relatively-normal donor/library score direction with contributor-retained cells disclosed | Age-independent or clinical-source-independent degeneration effect |
| GSE244889 | NP directional support | Severe-minus-mild direction in 4 versus 3 presumed keys | Confirmatory validation or age-independent severity effect |
| GSE251686 | Isolated incomplete exploratory check | Separate 2 versus 3 descriptive score/effect display, deliberately excluded from the default 20-effect summary | Repairing `GSM7986002`, decisive validation, replication, or a balanced severity result |
| GSE186542 | S8/S9 post hoc score-level support | Early Pfirrmann I-III versus advanced IV-V score contrast (3 versus 3) | Raw-count inference, patient-level replication, or confirmation while FPKM/TMM versus raw-count provenance remains unresolved |
| GSE167931 | S8 post hoc FPKM representation | Normal versus degenerated score contrast (4 versus 5); TPM is paired processing sensitivity | Counting TPM as another cohort, assuming patient-level independence, or treating normalized values as raw counts |
| GSE245147 | S9 source-family replacement sensitivity | Native Degenerated versus No-degenerated RPKM contrast (3 versus 3) after excluding passage/treatment arms | Independent validation, pooling with GSE167931, or causal interpretation; source-family independence is unverified |
| GSE56081 | Candidate only | Exact 5 versus 5 sample mapping is documented; the fixed Ensembl-113/GRCh38 audit found 28 failed candidate probes and no module passed the 0.80 global-specificity gate | Locked-module extension, a biological negative-result claim, or treating identical-RNA alternatives as independent samples |

## Current signal summary

The frozen default summary contains four NP score cohorts/roles, four modules,
and exploratory AF effects from the GSE230809 parent. The only default NP
module with the same positive higher-severity point-estimate direction in all
four cohorts is hypoxia/oxidative stress. This remains a **descriptive sign
pattern**, not a replicated or universal program: intervals, leave-one-key-out
stability, platform and processing differences, and cohort confounding must be
shown alongside it. ECM, inflammatory, and disc-matrix-homeostasis directions
are heterogeneous in the default layer.

S7 independently standardizes the four default NP contrasts and remains
non-confirmatory. S8 adds GSE186542 and GSE167931 FPKM, producing pooled SMDH
values of 0.7780 (ECM), 0.4032 (inflammatory), 0.7694 (hypoxia), and 0.3762
(homeostasis). The S8 hypoxia interval is 0.1706 to 1.3682 and its HKSJ p-value
is 0.0214 (four-module BH 0.0856), but this is a post hoc score-level result
with `k = 6` and unverified patient-level independence. S9 replaces GSE167931
with GSE245147 and yields SMDH 1.0931 (95% CI -0.9186 to 3.1048) for ECM,
0.4056 (-1.0004 to 1.8115) for inflammatory, 0.5746 (-0.7231 to 1.8723) for
hypoxia, and 0.1046 (-1.0200 to 1.2292) for homeostasis. The S9 BH values are
0.6133, 0.6556, 0.6133, and 0.8205, respectively. The loss of a positive,
non-zero-crossing hypoxia interval in S9 demonstrates source-family
sensitivity, not confirmation.

The wording “replicated,” “validated biomarker,” “driver,” “mechanism,”
“therapeutic target,” and “first discovered” is not allowed in the current
manuscript unless a future independent dataset and orthogonal evidence meet a
new, pre-specified gate.

## Required manuscript safeguards

1. Put the experimental unit in the abstract and Methods: donor or presumed
   donor/library key, never cell.
2. Report each cohort's group sizes, identity confidence, cell/library yields,
   mapping loss, excluded records, and processing boundary.
3. Keep GSE230809's two child series under one parent-project label.
4. Label GSE165722 normalized matrices as score-level only and preserve the
   GEO-versus-publication grade discrepancy.
5. State that GSE153066 normal/degenerate status is confounded with clinical
   source and age.
6. Display effect sizes and 95% intervals in the default layer; do not
   substitute p-values or cell-level tests. S7-S9 HKSJ/BH p-values are
   transparent, non-confirmatory descriptors only.
7. Treat default cross-cohort sign alignment as descriptive. S7-S9 may
   standardize one effect per cohort with SMDH, but their small `k`, different
   processing scales, and unverified patient-level independence prohibit
   confirmation, replication, biomarker, mechanism, causal, or treatment
   interpretation. Figure labels must describe default differences as
   cohort-specific, unitless score differences rather than a common CPM scale.
8. Keep `GSM7986002` out of every expression result and retain its integrity
   failure in the exclusion ledger.
9. Treat `data/derived/donor_module_effect_summary/` as the only authoritative
   default result path. `module_scores_recomputed/` is the canonical scored
   input location; `donor_module_effect_summary_canonical_audit/` is a legacy
   audit snapshot and must not be used as a default result.
10. Keep the separately scored GSE251686 package outside the default summary;
     it is intentionally not an additional cross-cohort effect or sign-alignment
     contributor.
11. Keep S7, S8, and S9 in their named supplementary directories. S8's
    GSE186542 and GSE167931 inputs are accession-level additions only, and S9's
    GSE245147 input replaces rather than supplements GSE167931 because source-
    family overlap cannot be excluded.

## Stopping criteria

Stop the current biological narrative if a score table has missing or
conflicting sample identities, mapping below 80%, non-finite values, an
unresolved malformed matrix, or a cohort direction that is opposite to the
pre-specified claim. A future confirmatory claim requires at least two
independent human validation contrasts with verified sample nesting,
pre-specified direction, interpretable intervals, and no unresolved opposite
direction; it also requires orthogonal biological evidence for any mechanism
or therapeutic wording. S7-S9 cannot satisfy this gate: they are post hoc,
score-level syntheses with small cohort counts and unverified patient-level
independence. Their p-values cannot be promoted to a stopping-rule exception.

## Reproducibility entry points

- Protocol: `docs/00_analysis_protocol.md`
- Locked modules: `config/program_modules.json`
- GSE165722 score audit: `docs/15_GSE165722_score_level_readiness.md`
- Default donor-level summary: `data/derived/donor_module_effect_summary/README.md`
- Default-summary run manifest: `data/derived/donor_module_effect_summary/run_manifest.json`
- GSE251686 readiness boundary: `docs/12_GSE251686_readiness.md`
- GSE251686 isolated exploratory scores/effects: `data/derived/GSE251686_exploratory_scores/`
- Candidate accession audit: `data/derived/geo_candidate_audit/`
- S7 exploratory synthesis: `data/derived/np_exploratory_meta_analysis/`
- S8 post hoc external expansion: `data/derived/np_post_hoc_external_expansion_meta_analysis/`
- S9 source-family replacement sensitivity: `data/derived/np_source_family_replacement_meta_analysis/`
