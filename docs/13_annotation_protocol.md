# Conservative Source-Restricted Cell Selection Protocol

## Purpose

This protocol defines the cell-level source-restricted selection used before
donor-level module aggregation. It is a filtering and provenance step, not a
claim that a cell's lineage or developmental origin has been established.

## Fixed evidence hierarchy

1. The GEO anatomical source (`AF` or `NP`) is the primary compartment label
   and cannot be overwritten by expression.
2. A fixed multi-gene support panel quantifies AF-like and NP-like
   transcriptional support. It is a concordance audit, not a cell-type gate.
   AF and NP matrix programs overlap in disc cells, especially with disease or
   tissue processing, so their coexpression does not establish a doublet.
3. Immune, endothelial, mural, and erythroid panels are exclusion evidence.
4. A source-marker discrepancy is retained as an audit flag. It does not by
   itself exclude a cell from the source-restricted compartment pseudobulk.
   A missing source label or strong nonresident evidence prevents inclusion.

## Locked thresholds

Scores are the mean log1p(CPM) over genes measured in the fixed panel, with CPM
computed from the cell's total UMI. Detection counts are the number of panel
genes with a positive count. The thresholds below are fixed before examining
condition-specific expression:

- resident support: at least 1 detected gene and score >= 0.2;
- nonresident exclusion: at least 3 detected genes and score >= 2.0;
- doublet flag: source-matched resident support together with a nonresident
  exclusion panel. AF/NP support coexpression alone is not a doublet criterion,
  because both are disc matrix programs;
- a source-marker discrepancy is flagged when the alternative source panel is
  more than 0.2 score units higher, but is not an exclusion criterion.

The exact gene lists and hashable JSON configuration are stored in
`config/cell_marker_panels.json`. Gene symbols absent from a dataset are
reported; they are not replaced with result-derived markers. This is protocol
amendment v1.1, made after pooled technical QC/marker inspection and before
any disease-group/module result was calculated. It replaces the earlier,
overly restrictive rule that treated AF/NP coexpression as incompatible.

## Labels

- `source_AF_nonexcluded`: source AF with no strong nonresident panel. AF/NP
  support and any source-marker discrepancy are reported separately.
- `source_NP_nonexcluded`: analogous rule for source NP.
- `nonresident_immune`, `nonresident_endothelial`, `nonresident_mural`, or
  `nonresident_erythroid`: the corresponding exclusion panel passes and no
  higher-confidence doublet rule supersedes it.
- `mixed_or_nonresident`: source support and a strong incompatible nonresident
  panel both pass.
- `ambiguous`: no source label is present.

Only source-labelled, non-excluded `source_AF_nonexcluded` or
`source_NP_nonexcluded` cells enter the source-restricted compartment
pseudobulk. A source-matched cell with absent support-panel genes is retained
under the source label and flagged as `support_insufficient`; it is not
silently relabelled or discarded. All exclusion and ambiguity labels are
retained for reporting and excluded. These are source-restricted compartment
pseudobulks, not claims of purified resident cell populations. `COL4A1` and
`COL4A2` are deliberately not used as mural exclusion markers because their
matrix expression is not specific enough in disc tissue.
Annotation proportions are descriptive and are not treated as independent
replicates.

## Sensitivity requirements

The primary resident pseudobulk requires at least 30 retained resident cells
per donor and compartment. Thresholds of 20 and 50 cells are reported as
sensitivity analyses. Any donor failing the threshold is excluded from that
specific pseudobulk, with no silent pooling across donors.
