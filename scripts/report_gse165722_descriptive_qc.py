"""Create a concise, non-inferential summary of processed-matrix QC."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


def fmt(value: float, digits: int = 1) -> str:
    return f"{value:.{digits}f}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("qc", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.qc.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 8:
        raise ValueError(f"Expected eight GSE165722 rows, found {len(rows)}")

    cells = np.array([int(row["cells"]) for row in rows])
    values = np.array([float(row["median_supplied_value"]) for row in rows])
    genes = np.array([float(row["median_detected_features"]) for row in rows])
    mt = np.array([float(row["median_mitochondrial_fraction"]) for row in rows])
    pass_fraction = np.array([
        float(row["fraction_meeting_source_like_descriptive_thresholds"])
        for row in rows
    ])

    report = f"""# GSE165722 Processed-Matrix Descriptive QC

Audit date: 2026-08-14. Source table: `data/derived/GSE165722_descriptive_qc.csv`.

## Scope

These are descriptive summaries of GEO-supplied matrices, which GEO labels
\"normalized counts\". The word \"supplied value\" is used deliberately: these
summaries do not establish raw UMI status and are not inputs to edgeR, DESeq2,
pseudobulk count models, cell filtering, or biological group inference.

## Matrix-level summary

- Eight presumed donor-level NP matrices contain {cells.sum():,} supplied cells.
- Per-matrix cell counts range from {cells.min():,} to {cells.max():,}.
- Per-matrix median supplied values range from {fmt(values.min())} to {fmt(values.max())}.
- Per-matrix median detected features range from {fmt(genes.min())} to {fmt(genes.max())}.
- Per-matrix median mitochondrial fractions range from {fmt(mt.min() * 100, 1)}% to {fmt(mt.max() * 100, 1)}%.
- The fraction meeting source-like descriptive cutoffs (supplied value >=500,
  detected features >=200, mitochondrial fraction <=20%) ranges from
  {fmt(pass_fraction.min() * 100, 1)}% to {fmt(pass_fraction.max() * 100, 1)}%.

## Interpretation limits

The range in supplied-library characteristics requires transparent reporting
and later sensitivity analyses. It is not evidence for a disease effect,
cell-type shift, or sample exclusion. The cohort remains eligible only for
score-level directional validation after gene sets and scoring rules are locked.
"""
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(report, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
