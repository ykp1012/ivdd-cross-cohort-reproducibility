"""Compute descriptive QC for GSE165722's GEO-described normalized matrices.

The matrix values are not labelled raw UMIs.  This script therefore reports
"supplied values" rather than counts and does not filter cells, normalize,
aggregate pseudobulk, cluster, or test groups.  It streams each compressed
matrix in the raw TAR without extraction.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import tarfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from audit_gse165722 import SAMPLES


def quantile(values: np.ndarray, probability: float) -> float:
    return float(np.quantile(values, probability))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    with tarfile.open(args.archive, "r") as tar:
        for sample in SAMPLES:
            count_name = f"{sample.gsm}_{sample.sample}.counts.tsv.gz"
            member = tar.extractfile(count_name)
            if member is None:
                raise FileNotFoundError(count_name)
            with gzip.open(member, "rt", encoding="utf-8", errors="replace") as handle:
                header = handle.readline().rstrip("\n\r").split("\t")
                cells = len(header) - 1
                if cells <= 0:
                    raise ValueError(f"No matrix columns in {count_name}")

                total_value = np.zeros(cells, dtype=np.float64)
                detected_features = np.zeros(cells, dtype=np.int32)
                mitochondrial_value = np.zeros(cells, dtype=np.float64)
                features = 0
                mitochondrial_features = 0
                nonzero_entries = 0

                for raw in handle:
                    line = raw.rstrip("\n\r")
                    if not line:
                        continue
                    gene, separator, value_text = line.partition("\t")
                    if not separator:
                        raise ValueError(f"Malformed matrix row in {count_name}: {gene!r}")
                    values = np.fromstring(value_text, sep="\t", dtype=np.float64)
                    if values.size != cells:
                        raise ValueError(
                            f"{count_name}: {gene} has {values.size} values, expected {cells}"
                        )
                    if np.any(values < 0) or not np.all(np.isfinite(values)):
                        raise ValueError(f"{count_name}: non-finite or negative supplied value at {gene}")
                    positive = values > 0
                    total_value += values
                    detected_features += positive
                    nonzero_entries += int(positive.sum())
                    features += 1
                    if gene.upper().startswith("MT-"):
                        mitochondrial_value += values
                        mitochondrial_features += 1

            mt_fraction = np.divide(
                mitochondrial_value,
                total_value,
                out=np.full(cells, np.nan, dtype=np.float64),
                where=total_value > 0,
            )
            source_like_qc_pass = (
                (detected_features >= 200)
                & (total_value >= 500)
                & (mt_fraction <= 0.20)
            )
            rows.append(
                {
                    "dataset": "GSE165722",
                    "gsm": sample.gsm,
                    "sample": sample.sample,
                    "presumed_donor_id": sample.sample,
                    "compartment": "NP",
                    "source_publication_severity_group": sample.severity_group,
                    "matrix_status": "GEO-described normalized counts; descriptive QC only",
                    "cells": cells,
                    "features": features,
                    "mitochondrial_features": mitochondrial_features,
                    "matrix_nonzero_fraction": nonzero_entries / (features * cells),
                    "matrix_zero_fraction": 1 - nonzero_entries / (features * cells),
                    "median_supplied_value": quantile(total_value, 0.50),
                    "q1_supplied_value": quantile(total_value, 0.25),
                    "q3_supplied_value": quantile(total_value, 0.75),
                    "median_detected_features": quantile(detected_features, 0.50),
                    "q1_detected_features": quantile(detected_features, 0.25),
                    "q3_detected_features": quantile(detected_features, 0.75),
                    "median_mitochondrial_fraction": quantile(mt_fraction, 0.50),
                    "q3_mitochondrial_fraction": quantile(mt_fraction, 0.75),
                    "fraction_meeting_source_like_descriptive_thresholds": float(source_like_qc_pass.mean()),
                    "threshold_note": "descriptive only: supplied value >=500; detected features >=200; mitochondrial fraction <=0.20",
                    "audited_at_utc": datetime.now(timezone.utc).isoformat(),
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["gsm"])
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
