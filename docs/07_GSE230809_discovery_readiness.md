# GSE230809 Discovery Data-Readiness Report

Audit began: 2026-08-13. Updated: 2026-08-14.

## Design verified from GEO SOFT

GSE229711 and GSE230808 are child series of one discovery project. The sample
ledger contains 24 libraries from 13 human donors: 3 low-grade/healthy and 10
advanced-degeneration donors. All donors are male. There are 11 AF/NP-paired
donors (3 healthy and 8 advanced), plus 2 advanced AF-only donors. The
canonical donor key is `patient id`.

Age is fully confounded with disease state: healthy donors are 21, 25, and 27
years old; advanced donors are 37-73 years old. Thus the discovery comparison
is an advanced-degeneration-associated contrast only. It cannot establish an
age-independent degeneration effect.

## Raw-archive status

| Child series | Expected 10x libraries | Status | Result |
|---|---:|---|---|
| GSE229711 | 6 | complete and audited | Each expected GSM has one barcode, feature, and matrix member. All 6 feature/barcode dimensions agree with their Matrix Market headers. |
| GSE230808 | 18 | complete and audited | SHA-256 `cbff3048581362d69ef5eb33130d5dd97b52a501af1b0a1789c2484012c811ec`; each expected GSM has one barcode, feature, and matrix member, and all 18 feature/barcode dimensions agree with their Matrix Market headers. |

The combined raw-data ledger contains all 24 expected libraries. Both child
series contain verified Cell Ranger 10x-style MTX triples and are eligible for
pre-QC raw-count ingestion, not yet for inference. GSE229711 is not analyzed
separately as a study; it is the healthy component of the single GSE230809
discovery project.

## Pre-QC observations requiring later reporting

The healthy NP libraries differ markedly in supplied barcode counts (1,872,
16,302, and 8,185). This is a library-yield observation, not a biological
result. Per-donor QC distributions, filtering rates, and downsampling
sensitivity analyses must be reported before any group contrast.

## Inference boundary and next gate

The raw-archive stop condition is cleared. Technical cell-level QC and the
conservative resident/source annotation ledgers have been generated for both
child archives. The next computational gate is a successful module-scoring run
with the QC and annotation ledgers, followed by donor-level effect-size and
leave-one-donor-out summaries.

The discovery dataset has only three low-grade donors per compartment, below the
protocol's minimum of four comparison donors for inference. Therefore it can
generate exploratory donor-level effect estimates, module-direction checks, and
leave-one-donor-out stability displays only. It cannot generate a confirmatory
differential-expression or module-significance claim. See
`docs/14_module_scoring_run.md` for the exact command and audit contract.
