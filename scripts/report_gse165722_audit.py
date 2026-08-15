"""Write a compact, evidence-first GSE165722 readiness report."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ledger", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.ledger.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    groups = Counter(row["source_publication_severity_group"] for row in rows)
    cells = sum(int(row["cells_in_supplied_matrix"]) for row in rows)

    report = f"""# GSE165722 Data-Quality and Analysis-Readiness Report

Audit date: 2026-08-13.

## Verified facts

- The GEO archive contains {len(rows)} presumed donor-level NP matrices and {cells:,} supplied cells in total.
- The cell-name files map supplied matrix column IDs to cell barcodes. Each GSM/Sample file is a presumed donor-level sample, but GEO does not provide a patient identifier, age, sex, disc level, or batch field; this nesting remains an explicit limitation.
- The source publication, Tu et al. (PMID 34825784), reports grades II-V. GEO SOFT labels the same ordered samples I-IV. This project preserves both fields and uses the source-publication grouping: mild II-III = {groups['mild']} donors, severe IV-V = {groups['severe']} donors.
- GEO sample metadata states that the supplementary files contain \"normalized counts\", even though inspected values are integer-like.

## Permitted role

GSE165722 is eligible for descriptive per-donor QC and pre-specified score-level direction checks. It is **not** eligible for raw-count pseudobulk aggregation, edgeR/DESeq2 inference, or a negative-binomial effect estimate unless independently verified raw UMI matrices are recovered.

## Stop conditions already resolved

- Donor mapping: provisional pass at sample level; one GSM/Sample file is one supplied matrix, but its biological-donor identity is not independently exposed in GEO metadata.
- Tissue label: pass; all samples are NP.
- Severity grouping: usable only with the cited source-publication mapping, not the GEO labels alone.
- Raw-count status: fail for primary count-model inference because GEO explicitly describes the values as normalized.

## Remaining before biological interpretation

1. Calculate matrix-level and cell-level descriptive QC from the supplied values without calling them UMI counts.
2. Lock module gene lists and score transformations before looking at group effects.
3. Use GSE165722 only as one component of cross-cohort directional validation; report its grade conflict and processing limitation in all methods and supplements.
"""
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(report, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
