"""Inspect the GSM7986002 nested TAR/GZIP Matrix Market payload.

This is a read-only archive-integrity audit for the one library whose
``matrix.mtx`` header was flagged for follow-up.  It parses the outer TAR
header and the decompressed inner TAR headers, records GNU/PAX sparse fields,
streams the matrix member without extraction, and counts NUL segments and
legal Matrix Market coordinate rows.  No expression result is calculated and
no raw file is modified.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path
from typing import BinaryIO


BLOCK = 512
TARGET_GSM = "GSM7986002"
TARGET_MEMBER_SUFFIX = "_NP3_EmptyDrops_CR_matrix.tar.gz"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def decode_field(raw: bytes) -> str:
    return raw.rstrip(b"\0 ").decode("ascii", errors="replace")


def parse_octal(raw: bytes) -> int | None:
    value = raw.strip(b"\0 ")
    if not value:
        return None
    try:
        return int(value, 8)
    except ValueError:
        return None


def checksum_info(header: bytes) -> dict[str, object]:
    stored = parse_octal(header[148:156])
    computed = sum(32 if 148 <= index < 156 else value for index, value in enumerate(header))
    return {
        "stored_checksum": stored,
        "computed_checksum": computed,
        "checksum_valid": stored == computed,
    }


def sparse_info(header: bytes) -> dict[str, object]:
    entries: list[dict[str, int | None]] = []
    for index in range(4):
        start = 386 + index * 24
        offset = parse_octal(header[start : start + 12])
        numbytes = parse_octal(header[start + 12 : start + 24])
        if offset not in (None, 0) or numbytes not in (None, 0):
            entries.append({"offset": offset, "numbytes": numbytes})
    return {
        "old_gnu_sparse_entries": entries,
        "old_gnu_sparse_isextended_raw": header[482],
        "old_gnu_sparse_isextended": bool(header[482] not in (0, ord("0"))),
        "old_gnu_sparse_realsize": parse_octal(header[483:495]),
    }


def parse_header(header: bytes) -> dict[str, object]:
    if len(header) != BLOCK:
        raise ValueError(f"TAR header must be {BLOCK} bytes, got {len(header)}")
    typeflag = header[156]
    result: dict[str, object] = {
        "name": decode_field(header[0:100]),
        "mode": decode_field(header[100:108]),
        "uid": decode_field(header[108:116]),
        "gid": decode_field(header[116:124]),
        "size_field": decode_field(header[124:136]),
        "size": parse_octal(header[124:136]),
        "mtime": decode_field(header[136:148]),
        "typeflag_byte": typeflag,
        "typeflag_ascii": chr(typeflag) if typeflag else "NUL",
        "linkname": decode_field(header[157:257]),
        "magic": decode_field(header[257:263]),
        "version": decode_field(header[263:265]),
        "uname": decode_field(header[265:297]),
        "gname": decode_field(header[297:329]),
        "devmajor": decode_field(header[329:337]),
        "devminor": decode_field(header[337:345]),
    }
    result.update(checksum_info(header))
    result.update(sparse_info(header))
    return result


def read_exact(stream: BinaryIO, count: int) -> bytes:
    chunks: list[bytes] = []
    remaining = count
    while remaining:
        chunk = stream.read(remaining)
        if not chunk:
            raise EOFError(f"Unexpected EOF while reading {count} bytes")
        chunks.append(chunk)
        remaining -= len(chunk)
    return b"".join(chunks)


def skip_exact(stream: BinaryIO, count: int) -> None:
    remaining = count
    while remaining:
        chunk = stream.read(min(1024 * 1024, remaining))
        if not chunk:
            raise EOFError(f"Unexpected EOF while skipping {count} bytes")
        remaining -= len(chunk)


def parse_pax(payload: bytes) -> dict[str, str]:
    result: dict[str, str] = {}
    cursor = 0
    while cursor < len(payload):
        separator = payload.find(b" ", cursor)
        if separator < 0:
            raise ValueError("Malformed PAX record without length")
        try:
            length = int(payload[cursor:separator])
        except ValueError as error:
            raise ValueError("Malformed PAX record length") from error
        record = payload[cursor : cursor + length]
        if not record.endswith(b"\n") or b"=" not in record:
            raise ValueError("Malformed PAX record payload")
        key, value = record[separator - cursor + 1 : -1].split(b"=", 1)
        result[key.decode("utf-8", errors="replace")] = value.decode("utf-8", errors="replace")
        cursor += length
    return result


class BoundedReader:
    """A binary reader limited to one TAR member's declared byte size."""

    def __init__(self, stream: BinaryIO, size: int):
        self.stream = stream
        self.remaining = size
        self.consumed = 0

    def read(self, count: int = -1) -> bytes:
        if self.remaining <= 0:
            return b""
        if count < 0 or count > self.remaining:
            count = self.remaining
        chunk = self.stream.read(count)
        if not chunk:
            raise EOFError("Unexpected EOF inside bounded TAR member")
        self.remaining -= len(chunk)
        self.consumed += len(chunk)
        return chunk

    def readline(self) -> bytes:
        if self.remaining <= 0:
            return b""
        chunks: list[bytes] = []
        while self.remaining:
            chunk = self.read(1)
            chunks.append(chunk)
            if chunk == b"\n":
                break
        return b"".join(chunks)

    def drain(self) -> None:
        while self.remaining:
            self.read(min(1024 * 1024, self.remaining))


def matrix_audit(stream: BoundedReader) -> dict[str, object]:
    header = stream.readline()
    dimensions = stream.readline()
    while dimensions.startswith(b"%"):
        dimensions = stream.readline()
    if not header.startswith(b"%%MatrixMarket matrix coordinate"):
        raise ValueError(f"Unexpected Matrix Market header: {header[:120]!r}")
    try:
        rows, columns, nnz = (int(value) for value in dimensions.split())
    except ValueError as error:
        raise ValueError(f"Invalid Matrix Market dimensions: {dimensions!r}") from error

    nul_count = 0
    nul_segments: list[dict[str, int]] = []
    current_nul_start: int | None = None
    last_nul_position: int | None = None
    non_ascii_count = 0
    control_byte_count = 0
    total_lines = 0
    blank_lines = 0
    comment_lines = 0
    legal_coordinate_lines = 0
    malformed_lines = 0
    out_of_bounds_lines = 0
    nonpositive_value_lines = 0
    line_with_nul_count = 0
    last_line_terminated = False
    coordinate_sample: list[str] = []

    def inspect_nuls(line: bytes, line_start: int) -> None:
        nonlocal nul_count, current_nul_start, last_nul_position
        for index, value in enumerate(line):
            if value == 0:
                position = line_start + index
                nul_count += 1
                last_nul_position = position
                if current_nul_start is None:
                    current_nul_start = position
            elif current_nul_start is not None:
                nul_segments.append(
                    {"start": current_nul_start, "end_exclusive": line_start + index}
                )
                current_nul_start = None

    while stream.remaining:
        line_start = stream.consumed
        line = stream.readline()
        if not line:
            break
        total_lines += 1
        last_line_terminated = line.endswith(b"\n")
        inspect_nuls(line, line_start)
        line_with_nul_count += int(b"\0" in line)
        non_ascii_count += sum(value >= 128 for value in line)
        control_byte_count += sum(value < 32 and value not in (9, 10, 13) for value in line)
        stripped = line.rstrip(b"\r\n")
        if not stripped:
            blank_lines += 1
            continue
        if stripped.startswith(b"%"):
            comment_lines += 1
            continue
        fields = stripped.split()
        try:
            if len(fields) != 3:
                raise ValueError
            row, column, value = (int(field) for field in fields)
        except ValueError:
            malformed_lines += 1
            continue
        if not (1 <= row <= rows and 1 <= column <= columns):
            out_of_bounds_lines += 1
        if value <= 0:
            nonpositive_value_lines += 1
        legal_coordinate_lines += 1
        if len(coordinate_sample) < 5:
            coordinate_sample.append(stripped.decode("ascii", errors="replace"))

    if current_nul_start is not None:
        nul_segments.append({"start": current_nul_start, "end_exclusive": stream.consumed})
    stream.drain()
    if len(nul_segments) > 100:
        nul_segments = nul_segments[:100]
    return {
        "matrix_header": header.rstrip(b"\r\n").decode("ascii", errors="replace"),
        "dimensions_line": dimensions.rstrip(b"\r\n").decode("ascii", errors="replace"),
        "matrix_rows": rows,
        "matrix_columns": columns,
        "matrix_nnz_header": nnz,
        "matrix_member_declared_bytes": stream.consumed,
        "matrix_member_bytes_consumed": stream.consumed,
        "nul_byte_count": nul_count,
        "nul_segment_count": len(nul_segments),
        "nul_segments_first_100": nul_segments,
        "first_nul_position": nul_segments[0]["start"] if nul_segments else None,
        "last_nul_position": last_nul_position,
        "non_ascii_byte_count": non_ascii_count,
        "unexpected_control_byte_count": control_byte_count,
        "total_lines_after_dimensions": total_lines,
        "blank_lines": blank_lines,
        "comment_lines_after_dimensions": comment_lines,
        "line_with_nul_count": line_with_nul_count,
        "legal_coordinate_lines": legal_coordinate_lines,
        "malformed_lines": malformed_lines,
        "out_of_bounds_lines": out_of_bounds_lines,
        "nonpositive_value_lines": nonpositive_value_lines,
        "coordinate_sample_first_5": coordinate_sample,
        "last_line_terminated": last_line_terminated,
        "line_count_matches_header": legal_coordinate_lines == nnz,
        "matrix_payload_legal": (
            legal_coordinate_lines == nnz
            and malformed_lines == 0
            and blank_lines == 0
            and out_of_bounds_lines == 0
            and nonpositive_value_lines == 0
            and nul_count == 0
            and non_ascii_count == 0
            and control_byte_count == 0
            and stream.remaining == 0
        ),
    }


def parse_inner_tar(stream: BinaryIO) -> dict[str, object]:
    members: list[dict[str, object]] = []
    pending_pax: dict[str, str] = {}
    pending_longname: str | None = None
    target_matrix: dict[str, object] | None = None
    inner_zero_blocks = 0
    inner_nonzero_trailing_bytes = 0
    inner_uncompressed_offset = 0
    while True:
        header_offset = inner_uncompressed_offset
        header = stream.read(BLOCK)
        if len(header) == 0:
            raise EOFError("Inner gzip stream ended before TAR end markers")
        if len(header) != BLOCK:
            raise EOFError("Truncated inner TAR header")
        inner_uncompressed_offset += BLOCK
        if header == b"\0" * BLOCK:
            inner_zero_blocks += 1
            second = read_exact(stream, BLOCK)
            inner_uncompressed_offset += BLOCK
            if second != b"\0" * BLOCK:
                raise ValueError("Inner TAR has one zero block followed by nonzero header")
            inner_zero_blocks += 1
            trailing = stream.read()
            inner_nonzero_trailing_bytes = sum(value != 0 for value in trailing)
            break
        parsed = parse_header(header)
        raw_name = str(parsed["name"])
        typeflag = int(parsed["typeflag_byte"])
        size = parsed["size"]
        if not isinstance(size, int):
            raise ValueError(f"Inner member has invalid size: {raw_name}")
        payload = BoundedReader(stream, size)
        if typeflag in (ord("x"), ord("g")):
            pax_payload = payload.read()
            payload.drain()
            pending_pax.update(parse_pax(pax_payload))
            inner_uncompressed_offset += size
            inner_uncompressed_offset += (-size) % BLOCK
            skip_exact(stream, (-size) % BLOCK)
            continue
        if typeflag == ord("L"):
            pending_longname = payload.read().rstrip(b"\0\n").decode("utf-8", errors="replace")
            payload.drain()
            inner_uncompressed_offset += size
            inner_uncompressed_offset += (-size) % BLOCK
            skip_exact(stream, (-size) % BLOCK)
            continue
        effective_name = pending_longname or str(pending_pax.get("path", raw_name))
        effective_pax = dict(pending_pax)
        parsed["effective_name"] = effective_name
        parsed["pax_headers"] = effective_pax
        parsed["header_offset"] = header_offset
        parsed["data_offset"] = header_offset + BLOCK
        parsed["data_size"] = size
        parsed["header_checksum_valid"] = bool(parsed["checksum_valid"])
        is_matrix = effective_name.endswith("matrix.mtx")
        if is_matrix:
            matrix_result = matrix_audit(payload)
            target_matrix = {**parsed, **matrix_result}
        else:
            payload.drain()
        inner_uncompressed_offset += size
        padding = (-size) % BLOCK
        skip_exact(stream, padding)
        inner_uncompressed_offset += padding
        members.append(parsed)
        pending_pax = {}
        pending_longname = None
    if target_matrix is None:
        raise ValueError("Target nested TAR did not contain matrix.mtx")
    return {
        "members": members,
        "target_matrix": target_matrix,
        "inner_zero_blocks": inner_zero_blocks,
        "inner_nonzero_trailing_bytes": inner_nonzero_trailing_bytes,
        "inner_tar_end_marker_valid": inner_zero_blocks == 2 and inner_nonzero_trailing_bytes == 0,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    archive = args.archive.resolve()
    outer_record: dict[str, object] | None = None
    inner_result: dict[str, object] | None = None
    outer_tar_read_pass = False
    with archive.open("rb") as raw:
        raw_header: bytes | None = None
        with tarfile.open(fileobj=raw, mode="r:") as outer:
            target = next((member for member in outer.getmembers() if member.name == f"{TARGET_GSM}{TARGET_MEMBER_SUFFIX}"), None)
            if target is None:
                raise FileNotFoundError(f"No outer member for {TARGET_GSM}")
            raw.seek(target.offset)
            raw_header = read_exact(raw, BLOCK)
            outer_record = parse_header(raw_header)
            outer_record["tarinfo_type"] = target.type.decode("ascii", errors="replace")
            outer_record["tarinfo_size"] = target.size
            outer_record["tarinfo_offset"] = target.offset
            outer_record["tarinfo_offset_data"] = target.offset_data
            outer_record["tarinfo_sparse"] = target.sparse
            outer_record["tarinfo_pax_headers"] = target.pax_headers
            stream = outer.extractfile(target)
            if stream is None:
                raise OSError(target.name)
            with stream:
                compressed = stream.read()
            gzip_header = compressed[:10]
            outer_record["member_bytes_read"] = len(compressed)
            outer_record["member_bytes_match_header"] = len(compressed) == target.size
            outer_record["member_bytes_sha256"] = hashlib.sha256(compressed).hexdigest()
            outer_record["gzip_header_hex_first_10"] = gzip_header.hex()
            outer_record["gzip_magic_valid"] = gzip_header[:2] == b"\x1f\x8b"
            outer_record["gzip_method"] = gzip_header[2] if len(gzip_header) > 2 else None
            outer_record["gzip_flags"] = gzip_header[3] if len(gzip_header) > 3 else None
            gz = gzip.GzipFile(fileobj=io.BytesIO(compressed), mode="rb")
            try:
                inner_result = parse_inner_tar(gz)
                gzip_trailing = gz.read()
                outer_record["gzip_trailing_decompressed_bytes"] = len(gzip_trailing)
                outer_record["gzip_crc_eof_read_pass"] = True
            finally:
                gz.close()
            outer_tar_read_pass = True

    if outer_record is None or inner_result is None:
        raise RuntimeError("Audit did not produce records")
    output = {
        "audit": {
            "dataset": "GSE251686",
            "gsm": TARGET_GSM,
            "archive": archive.name,
            "archive_bytes": archive.stat().st_size,
            "archive_sha256": sha256(archive),
            "outer_tar_read_pass": outer_tar_read_pass,
            "interpretation": "No GNU sparse representation or archive corruption detected if all pass flags are true; NULs are assessed within the declared matrix member only.",
        },
        "outer_member": outer_record,
        "inner_members": inner_result["members"],
        "inner_tar": {
            "member_count": len(inner_result["members"]),
            "inner_zero_blocks": inner_result["inner_zero_blocks"],
            "inner_nonzero_trailing_bytes": inner_result["inner_nonzero_trailing_bytes"],
            "inner_tar_end_marker_valid": inner_result["inner_tar_end_marker_valid"],
        },
        "matrix_member": inner_result["target_matrix"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, ensure_ascii=True), encoding="utf-8")

    # A compact CSV is useful for quick review without loading the JSON.
    matrix = inner_result["target_matrix"]
    assert isinstance(matrix, dict)
    csv_path = args.output.with_suffix(".csv")
    csv_fields = [
        "dataset", "gsm", "outer_member_typeflag", "outer_member_size", "outer_member_sparse",
        "outer_member_pax_headers", "outer_member_checksum_valid", "inner_member_typeflag",
        "inner_member_size", "inner_member_sparse", "inner_member_pax_headers",
        "inner_member_checksum_valid", "matrix_nnz_header", "coordinate_lines_observed",
        "nul_byte_count", "nul_segment_count", "first_nul_position", "last_nul_position",
        "malformed_lines", "out_of_bounds_lines", "nonpositive_value_lines", "line_count_matches_header",
        "matrix_payload_legal", "gzip_crc_eof_read_pass", "inner_tar_end_marker_valid",
    ]
    row = {
        "dataset": "GSE251686",
        "gsm": TARGET_GSM,
        "outer_member_typeflag": outer_record["typeflag_ascii"],
        "outer_member_size": outer_record["size"],
        "outer_member_sparse": json.dumps(outer_record["tarinfo_sparse"], ensure_ascii=True),
        "outer_member_pax_headers": json.dumps(outer_record["tarinfo_pax_headers"], ensure_ascii=True),
        "outer_member_checksum_valid": outer_record["checksum_valid"],
        "inner_member_typeflag": matrix["typeflag_ascii"],
        "inner_member_size": matrix["size"],
        "inner_member_sparse": json.dumps(matrix["old_gnu_sparse_entries"], ensure_ascii=True),
        "inner_member_pax_headers": json.dumps(matrix["pax_headers"], ensure_ascii=True),
        "inner_member_checksum_valid": matrix["checksum_valid"],
        "matrix_nnz_header": matrix["matrix_nnz_header"],
        "coordinate_lines_observed": matrix["legal_coordinate_lines"],
        "nul_byte_count": matrix["nul_byte_count"],
        "nul_segment_count": matrix["nul_segment_count"],
        "first_nul_position": matrix["first_nul_position"],
        "last_nul_position": matrix["last_nul_position"],
        "malformed_lines": matrix["malformed_lines"],
        "out_of_bounds_lines": matrix["out_of_bounds_lines"],
        "nonpositive_value_lines": matrix["nonpositive_value_lines"],
        "line_count_matches_header": matrix["line_count_matches_header"],
        "matrix_payload_legal": matrix["matrix_payload_legal"],
        "gzip_crc_eof_read_pass": outer_record["gzip_crc_eof_read_pass"],
        "inner_tar_end_marker_valid": inner_result["inner_tar_end_marker_valid"],
    }
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=csv_fields)
        writer.writeheader()
        writer.writerow(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
