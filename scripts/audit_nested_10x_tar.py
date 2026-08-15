"""Audit GEO TAR archives whose per-GSM payload is a nested 10x TAR.GZ.

The script reads the outer archive and each inner archive in memory only long
enough to inspect its three members. It validates structural dimensions and
writes a machine-readable inventory; it does not quantify expression.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import re
import tarfile
from pathlib import Path


GSM_RE = re.compile(r"^(GSM\d+)_")


def count_nonempty(stream: io.BufferedReader) -> int:
    with io.TextIOWrapper(stream, encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if line.strip())


def read_matrix_header(stream: io.BufferedReader) -> tuple[int, int, int]:
    with io.TextIOWrapper(stream, encoding="utf-8", errors="replace") as handle:
        first = handle.readline().strip()
        if not first.startswith("%%MatrixMarket matrix coordinate"):
            raise ValueError(f"Unexpected Matrix Market header: {first!r}")
        line = handle.readline().strip()
        while line.startswith("%"):
            line = handle.readline().strip()
        return tuple(int(value) for value in line.split())  # type: ignore[return-value]


def open_member(inner: tarfile.TarFile, matcher: re.Pattern[str]) -> tarfile.TarInfo:
    matches = [m for m in inner.getmembers() if m.isfile() and matcher.search(Path(m.name).name)]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one {matcher.pattern} in nested archive, found {[m.name for m in matches]}")
    return matches[0]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows: list[dict[str, object]] = []
    with tarfile.open(args.archive, "r") as outer:
        for outer_member in sorted((m for m in outer.getmembers() if m.isfile()), key=lambda m: m.name):
            match = GSM_RE.match(outer_member.name)
            if not match:
                raise ValueError(f"Outer member lacks GSM prefix: {outer_member.name}")
            stream = outer.extractfile(outer_member)
            if stream is None:
                raise OSError(outer_member.name)
            with stream:
                compressed = stream.read()
            with gzip.GzipFile(fileobj=io.BytesIO(compressed), mode="rb") as payload:
                nested_bytes = payload.read()
            with tarfile.open(fileobj=io.BytesIO(nested_bytes), mode="r:") as inner:
                barcode_member = open_member(inner, re.compile(r"(^|/)barcodes\.tsv$"))
                feature_member = open_member(inner, re.compile(r"(^|/)(genes|features)\.tsv$"))
                matrix_member = open_member(inner, re.compile(r"(^|/)matrix\.mtx$"))
                barcode_stream = inner.extractfile(barcode_member)
                feature_stream = inner.extractfile(feature_member)
                matrix_stream = inner.extractfile(matrix_member)
                if barcode_stream is None or feature_stream is None or matrix_stream is None:
                    raise OSError(outer_member.name)
                barcode_count = count_nonempty(barcode_stream)
                feature_count = count_nonempty(feature_stream)
                matrix_rows, matrix_columns, matrix_nnz = read_matrix_header(matrix_stream)
                dimension_ok = feature_count == matrix_rows and barcode_count == matrix_columns
                if not dimension_ok:
                    raise ValueError(
                        f"{outer_member.name}: feature/barcode {feature_count}/{barcode_count} "
                        f"do not match matrix {matrix_rows}/{matrix_columns}"
                    )
                rows.append(
                    {
                        "gsm": match.group(1), "outer_member": outer_member.name,
                        "outer_member_bytes": outer_member.size,
                        "nested_barcode_member": barcode_member.name,
                        "nested_feature_member": feature_member.name,
                        "nested_matrix_member": matrix_member.name,
                        "features": feature_count, "barcodes": barcode_count,
                        "matrix_rows": matrix_rows, "matrix_columns": matrix_columns,
                        "matrix_nnz": matrix_nnz, "dimension_check_pass": dimension_ok,
                        "matrix_status": "nested per-GSM matrix emitted by CeleScope; raw-count eligibility pending methodological review",
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
