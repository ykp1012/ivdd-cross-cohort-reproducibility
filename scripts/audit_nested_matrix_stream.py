"""Full-stream integrity audit for per-GSM nested Matrix Market archives.

This audit is deliberately separate from biological QC.  It checks the
decompressed matrix text byte-for-byte enough to detect embedded NUL bytes,
malformed coordinate records, dimension/range/count mismatches, and ordinary
versus sparse TAR metadata without rewriting the source archive.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
import shutil
import tarfile
import tempfile
from pathlib import Path


GSM_RE = re.compile(r"^(GSM\d+)_")
REQUIRED = {"genes.tsv", "barcodes.tsv", "matrix.mtx"}


def _read_text_member(stream: tarfile.ExFileObject) -> list[bytes]:
    with stream:
        return [line.rstrip(b"\r\n") for line in stream if line.strip()]


def _member_name_map(inner: tarfile.TarFile) -> dict[str, tarfile.TarInfo]:
    result: dict[str, tarfile.TarInfo] = {}
    for member in inner:
        if member.isfile():
            base = Path(member.name).name
            if base in result:
                raise ValueError(f"duplicate nested member basename: {base}")
            result[base] = member
    missing = REQUIRED - set(result)
    if missing:
        raise ValueError(f"nested archive missing members: {sorted(missing)}")
    return result


def _parse_int_triplet(raw: bytes) -> tuple[int, int, int] | None:
    fields = raw.split()
    if len(fields) < 3:
        return None
    try:
        return tuple(int(value) for value in fields[:3])  # type: ignore[return-value]
    except ValueError:
        return None


def audit_matrix(
    stream,
    n_features: int,
    n_barcodes: int,
    member_name: str,
) -> dict[str, object]:
    """Consume a matrix stream and return integrity counters."""
    byte_offset = 0
    first = stream.readline()
    byte_offset += len(first)
    if not first.startswith(b"%%MatrixMarket matrix coordinate"):
        raise ValueError(f"{member_name}: unexpected Matrix Market header {first[:80]!r}")

    dimensions = stream.readline()
    byte_offset += len(dimensions)
    while dimensions.startswith(b"%"):
        dimensions = stream.readline()
        byte_offset += len(dimensions)
    try:
        matrix_rows, matrix_columns, declared_nnz = (int(v) for v in dimensions.split())
    except ValueError as error:
        raise ValueError(f"{member_name}: invalid dimensions {dimensions[:100]!r}") from error

    counters = {
        "matrix_rows": matrix_rows,
        "matrix_columns": matrix_columns,
        "matrix_nnz_header": declared_nnz,
        "coordinate_lines_observed": 0,
        "valid_coordinate_lines": 0,
        "malformed_coordinate_lines": 0,
        "out_of_range_coordinates": 0,
        "negative_values": 0,
        "zero_values": 0,
        "nul_bytes": 0,
        "nul_lines": 0,
        "first_nul_offset": "",
        "last_nul_offset": "",
        "first_malformed_offset": "",
        "last_malformed_offset": "",
        "first_nul_context_before": "",
        "first_nul_context_after": "",
        "last_coordinate_byte_offset": byte_offset,
    }
    previous_tail = b""
    for raw in stream:
        start = byte_offset
        byte_offset += len(raw)
        counters["last_coordinate_byte_offset"] = byte_offset
        stripped = raw.strip()
        if not stripped or raw.startswith(b"%"):
            previous_tail = raw[-80:]
            continue
        counters["coordinate_lines_observed"] += 1
        nul_positions = [i for i, value in enumerate(raw) if value == 0]
        if nul_positions:
            counters["nul_bytes"] += len(nul_positions)
            counters["nul_lines"] += 1
            if not counters["first_nul_offset"]:
                counters["first_nul_offset"] = start + nul_positions[0]
                counters["first_nul_context_before"] = (previous_tail + raw[:nul_positions[0]])[-120:].decode("utf-8", "replace")
                counters["first_nul_context_after"] = raw[nul_positions[-1] + 1 : nul_positions[-1] + 121].decode("utf-8", "replace")
            counters["last_nul_offset"] = start + nul_positions[-1]
        triplet = _parse_int_triplet(stripped)
        if triplet is None or nul_positions:
            counters["malformed_coordinate_lines"] += 1
            if not counters["first_malformed_offset"]:
                counters["first_malformed_offset"] = start
            counters["last_malformed_offset"] = start
            previous_tail = raw[-80:]
            continue
        feature, barcode, value = triplet
        counters["valid_coordinate_lines"] += 1
        if not (1 <= feature <= n_features and 1 <= barcode <= n_barcodes):
            counters["out_of_range_coordinates"] += 1
        if value < 0:
            counters["negative_values"] += 1
        elif value == 0:
            counters["zero_values"] += 1
        previous_tail = raw[-80:]

    counters["dimension_match"] = matrix_rows == n_features and matrix_columns == n_barcodes
    counters["line_count_matches_header"] = counters["coordinate_lines_observed"] == declared_nnz
    counters["valid_count_matches_header"] = counters["valid_coordinate_lines"] == declared_nnz
    counters["text_integrity_pass"] = bool(
        counters["dimension_match"]
        and counters["line_count_matches_header"]
        and counters["valid_count_matches_header"]
        and counters["malformed_coordinate_lines"] == 0
        and counters["out_of_range_coordinates"] == 0
        and counters["negative_values"] == 0
        and counters["nul_bytes"] == 0
    )
    return counters


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--json", type=Path, default=None)
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    with tarfile.open(args.archive, mode="r|") as outer:
        for outer_member in outer:
            if not outer_member.isfile():
                continue
            match = GSM_RE.match(outer_member.name)
            if not match:
                raise ValueError(f"outer member lacks GSM prefix: {outer_member.name}")
            outer_stream = outer.extractfile(outer_member)
            if outer_stream is None:
                raise OSError(outer_member.name)
            with outer_stream, gzip.GzipFile(fileobj=outer_stream, mode="rb") as decompressed:
                with tempfile.TemporaryDirectory(prefix="ivdd_nested_audit_") as temp_dir:
                    temp_dir_path = Path(temp_dir)
                    members: dict[str, tarfile.TarInfo] = {}
                    paths: dict[str, Path] = {}
                    # A streaming tar reader cannot seek backwards.  Spool only
                    # the three required members, then audit them from disk.
                    with tarfile.open(fileobj=decompressed, mode="r|") as inner:
                        for nested in inner:
                            if not nested.isfile():
                                continue
                            base = Path(nested.name).name
                            if base in members:
                                raise ValueError(f"{outer_member.name}: duplicate nested basename {base}")
                            members[base] = nested
                            if base in REQUIRED:
                                extracted = inner.extractfile(nested)
                                if extracted is None:
                                    raise OSError(nested.name)
                                target = temp_dir_path / base
                                with extracted, target.open("wb") as output:
                                    shutil.copyfileobj(extracted, output, length=1024 * 1024)
                                paths[base] = target
                    missing = REQUIRED - set(members)
                    if missing:
                        raise ValueError(f"{outer_member.name}: missing {sorted(missing)}")
                    n_features = sum(1 for line in paths["genes.tsv"].open("rb") if line.strip())
                    n_barcodes = sum(1 for line in paths["barcodes.tsv"].open("rb") if line.strip())
                    with paths["matrix.mtx"].open("rb") as matrix_stream:
                        matrix_info = audit_matrix(matrix_stream, n_features, n_barcodes, outer_member.name)
                    row: dict[str, object] = {
                        "gsm": match.group(1),
                        "outer_member": outer_member.name,
                        "outer_member_bytes": outer_member.size,
                        "outer_typeflag": outer_member.type.decode("ascii", "replace"),
                        "outer_pax_headers": json.dumps(outer_member.pax_headers, sort_keys=True),
                        "nested_matrix_member": members["matrix.mtx"].name,
                        "nested_matrix_bytes": members["matrix.mtx"].size,
                        "nested_matrix_typeflag": members["matrix.mtx"].type.decode("ascii", "replace"),
                        "nested_matrix_pax_headers": json.dumps(members["matrix.mtx"].pax_headers, sort_keys=True),
                        "features_from_genes_tsv": n_features,
                        "barcodes_from_barcodes_tsv": n_barcodes,
                    }
                    row.update(matrix_info)
                    rows.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["gsm"]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(rows, indent=2, ensure_ascii=False), encoding="utf-8")
    failed = [row["gsm"] for row in rows if not row.get("text_integrity_pass")]
    print(json.dumps({"libraries": len(rows), "failed_libraries": failed}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
