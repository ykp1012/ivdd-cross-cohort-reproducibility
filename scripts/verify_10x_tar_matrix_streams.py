"""Verify complete Matrix Market gzip streams inside a GEO 10x TAR archive.

The companion ``audit_10x_tar.py`` checks feature/barcode dimensions against
Matrix Market headers.  This script additionally drains every matrix stream to
EOF (thereby exercising gzip CRC validation) and compares the declared NNZ
with the number of serialized coordinate lines.  It is an archive-integrity
check only; no expression values are analyzed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import tarfile
from pathlib import Path


def matrix_stream_audit(member: object) -> dict[str, object]:
    with gzip.open(member, "rb") as handle:
        first = handle.readline().strip()
        if not first.startswith(b"%%MatrixMarket matrix coordinate"):
            raise ValueError(f"Unexpected Matrix Market header: {first!r}")
        dimensions = handle.readline().strip()
        while dimensions.startswith(b"%"):
            dimensions = handle.readline().strip()
        try:
            rows, columns, nnz = (int(value) for value in dimensions.split())
        except ValueError as error:
            raise ValueError(f"Invalid Matrix Market dimensions: {dimensions!r}") from error

        coordinate_newlines = 0
        empty_coordinate_line = False
        previous_last = b""
        final_byte = b""
        while chunk := handle.read(1024 * 1024):
            coordinate_newlines += chunk.count(b"\n")
            empty_coordinate_line = empty_coordinate_line or b"\n\n" in chunk or (
                previous_last == b"\n" and chunk.startswith(b"\n")
            )
            previous_last = chunk[-1:]
            final_byte = chunk[-1:]

    coordinate_lines = coordinate_newlines + int(bool(final_byte) and final_byte != b"\n")
    return {
        "matrix_header": first.decode("ascii", errors="replace"),
        "matrix_rows": rows,
        "matrix_columns": columns,
        "matrix_nnz_header": nnz,
        "coordinate_lines_observed": coordinate_lines,
        "empty_coordinate_line_detected": empty_coordinate_line,
        "final_coordinate_line_terminated": final_byte == b"\n",
        "full_gzip_stream_read": True,
        "nnz_line_count_pass": coordinate_lines == nnz and not empty_coordinate_line,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("matrix_audit", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.matrix_audit.open(newline="", encoding="utf-8") as handle:
        expected_rows = list(csv.DictReader(handle))
    required = {"dataset", "gsm", "matrix_member", "matrix_rows", "matrix_columns", "matrix_nnz"}
    missing = required - set(expected_rows[0]) if expected_rows else required
    if missing:
        raise ValueError(f"Missing required audit columns: {sorted(missing)}")

    results: list[dict[str, object]] = []
    with tarfile.open(args.archive, "r") as tar:
        for expected in expected_rows:
            matrix_member = expected["matrix_member"]
            member = tar.extractfile(matrix_member)
            if member is None:
                raise FileNotFoundError(matrix_member)
            observed = matrix_stream_audit(member)
            prior_header_matches = (
                int(expected["matrix_rows"]) == observed["matrix_rows"]
                and int(expected["matrix_columns"]) == observed["matrix_columns"]
                and int(expected["matrix_nnz"]) == observed["matrix_nnz_header"]
            )
            results.append(
                {
                    "dataset": expected["dataset"],
                    "gsm": expected["gsm"],
                    "matrix_member": matrix_member,
                    **observed,
                    "matrix_header_matches_prior_audit": prior_header_matches,
                    "stream_integrity_pass": bool(observed["nnz_line_count_pass"]) and prior_header_matches,
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]) if results else ["gsm"])
        writer.writeheader()
        writer.writerows(results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
