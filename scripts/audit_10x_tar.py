"""Audit a GEO TAR containing prefixed 10x MTX triples without extraction.

This verifies that every ledger GSM has exactly one barcode, feature, and
matrix member, and streams Matrix Market headers for dimensions/NNZ. It writes
an audit table only; it does not normalize, filter, or analyze expression.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import tarfile
from pathlib import Path


SUFFIXES = {
    "barcodes": "_barcodes.tsv.gz",
    "features": "_features.tsv.gz",
    "matrix": "_matrix.mtx.gz",
}


def count_lines(member: object) -> int:
    with gzip.open(member, "rt", encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if line.strip())


def matrix_header(member: object) -> tuple[int, int, int]:
    with gzip.open(member, "rt", encoding="utf-8", errors="replace") as handle:
        first = handle.readline().strip()
        if not first.startswith("%%MatrixMarket matrix coordinate"):
            raise ValueError(f"Unexpected Matrix Market header: {first!r}")
        line = handle.readline().strip()
        while line.startswith("%"):
            line = handle.readline().strip()
        genes, cells, nnz = (int(value) for value in line.split())
        return genes, cells, nnz


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("ledger", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.ledger.open(newline="", encoding="utf-8") as handle:
        ledger = list(csv.DictReader(handle))
    result_rows: list[dict[str, object]] = []

    with tarfile.open(args.archive, "r") as tar:
        members = {member.name: member for member in tar.getmembers() if member.isfile()}
        archive_gsms = {name.split("_", maxsplit=1)[0] for name in members if re.match(r"^GSM\d+_", name)}
        archive_ledger = [row for row in ledger if row["gsm"] in archive_gsms]
        if not archive_ledger:
            raise ValueError(f"No ledger GSMs are represented in {args.archive.name}")
        for row in archive_ledger:
            gsm = row["gsm"]
            matches: dict[str, list[str]] = {}
            for label, suffix in SUFFIXES.items():
                matches[label] = sorted(name for name in members if name.startswith(f"{gsm}_") and name.endswith(suffix))
                if len(matches[label]) != 1:
                    raise ValueError(f"{gsm}: expected exactly one {label} member, found {matches[label]}")
            barcode_member = tar.extractfile(members[matches["barcodes"][0]])
            feature_member = tar.extractfile(members[matches["features"][0]])
            matrix_member = tar.extractfile(members[matches["matrix"][0]])
            if barcode_member is None or feature_member is None or matrix_member is None:
                raise FileNotFoundError(gsm)
            barcode_count = count_lines(barcode_member)
            feature_count = count_lines(feature_member)
            matrix_genes, matrix_cells, matrix_nnz = matrix_header(matrix_member)
            if barcode_count != matrix_cells or feature_count != matrix_genes:
                raise ValueError(
                    f"{gsm}: feature/barcode dimensions disagree with matrix "
                    f"({feature_count}/{barcode_count} vs {matrix_genes}/{matrix_cells})"
                )
            result_rows.append(
                {
                    "dataset": row["dataset"],
                    "gsm": gsm,
                    "donor_id": row["donor_id"],
                    "compartment": row["compartment"],
                    "disease_state": row["disease_state"],
                    "barcode_member": matches["barcodes"][0],
                    "feature_member": matches["features"][0],
                    "matrix_member": matches["matrix"][0],
                    "features": feature_count,
                    "barcodes": barcode_count,
                    "matrix_rows": matrix_genes,
                    "matrix_columns": matrix_cells,
                    "matrix_nnz": matrix_nnz,
                    "dimension_check_pass": True,
                    "matrix_status": "raw 10x Cell Ranger MTX triple reported by GEO",
                }
            )

        expected_gsms = {row["gsm"] for row in archive_ledger}
        unexpected = sorted(archive_gsms - expected_gsms)
        missing = sorted(expected_gsms - archive_gsms)
        if unexpected or missing:
            raise ValueError(f"Ledger/archive GSM mismatch: unexpected={unexpected}, missing={missing}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result_rows[0]) if result_rows else ["gsm"])
        writer.writeheader()
        writer.writerows(result_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
