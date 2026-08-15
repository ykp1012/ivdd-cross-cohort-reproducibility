"""Audit the GSE153066 combined count matrix by streaming its gzip payload.

GEO supplies one dense, sample-prefixed TSV rather than per-sample 10x
triples.  This script never extracts or materializes the expression matrix.
It verifies the tabular shape, gene identifiers, and cell-barcode prefixes,
then performs a deterministic sampled integer/non-negative value check and
links each prefix to a GEO SOFT sample.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


SAMPLE_START = re.compile(r"^\^SAMPLE\s*=\s*(?P<gsm>\S+)")
FIELD = re.compile(r"^!(?P<field>Sample_[^=]+)\s*=\s*(?P<value>.*)$")
PREFIX = re.compile(r"^(?P<prefix>(?:CTL|IDD)\d+)_")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_key(value: str) -> str:
    return re.sub(r"\s+", " ", value.strip().lower()).rstrip(":")


def split_characteristic(value: str) -> tuple[str, str] | None:
    if ":" not in value:
        return None
    key, content = value.split(":", maxsplit=1)
    return normalize_key(key), content.strip()


def first(values: dict[str, list[str]], key: str) -> str:
    return values.get(key, [""])[0]


def parse_soft(path: Path) -> list[dict[str, str]]:
    current: dict[str, object] | None = None
    samples: list[dict[str, object]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            start = SAMPLE_START.match(line)
            if start:
                if current is not None:
                    samples.append(current)
                current = {"gsm": start.group("gsm"), "fields": {}}
                continue
            if current is None:
                continue
            match = FIELD.match(line)
            if match:
                fields = current["fields"]
                assert isinstance(fields, dict)
                fields.setdefault(match.group("field").strip(), []).append(match.group("value"))
    if current is not None:
        samples.append(current)

    rows: list[dict[str, str]] = []
    for sample in samples:
        fields = sample["fields"]
        assert isinstance(fields, dict)
        characteristics: dict[str, list[str]] = {}
        for value in fields.get("Sample_characteristics_ch1", []):
            parsed = split_characteristic(value)
            if parsed is not None:
                characteristics.setdefault(parsed[0], []).append(parsed[1])
        rows.append(
            {
                "gsm": str(sample["gsm"]),
                "sample_title": first(fields, "Sample_title"),
                "tissue": first(characteristics, "tissue"),
                "status": first(characteristics, "status"),
                "age": first(characteristics, "age"),
                "sex": first(characteristics, "gender"),
                "raw_characteristics": " | ".join(fields.get("Sample_characteristics_ch1", [])),
                "data_processing": " | ".join(fields.get("Sample_data_processing", [])),
            }
        )
    return rows


def audited_value_positions(total_columns: int) -> set[int]:
    """Return deterministic, evenly distributed columns for sampled type checks."""
    if total_columns <= 256:
        return set(range(total_columns))
    spacing = max(1, total_columns // 256)
    positions = set(range(0, total_columns, spacing))
    positions.add(total_columns - 1)
    return positions


def audit_value_sample(values: list[str], positions: set[int], line_number: int) -> None:
    """Check a representative distributed sample without retaining numeric data."""
    for position in positions:
        value = values[position]
        try:
            integer = int(value)
        except ValueError as error:
            raise ValueError(f"row {line_number}: non-integer matrix value {value!r}") from error
        if integer < 0:
            raise ValueError(f"row {line_number}: negative matrix value {integer}")


def is_audited_gene_row(nonempty_gene_row_index: int) -> bool:
    """Sample the first row and every 500th non-empty gene row thereafter."""
    return nonempty_gene_row_index == 1 or nonempty_gene_row_index % 500 == 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path)
    parser.add_argument("soft", type=Path)
    parser.add_argument("--matrix-audit-output", type=Path, required=True)
    parser.add_argument("--sample-ledger-output", type=Path, required=True)
    args = parser.parse_args()

    matrix = args.matrix.resolve()
    soft = args.soft.resolve()
    if not matrix.is_file():
        raise FileNotFoundError(matrix)
    if not soft.is_file():
        raise FileNotFoundError(soft)

    soft_samples = parse_soft(soft)
    samples_by_title = {sample["sample_title"]: sample for sample in soft_samples}
    if len(samples_by_title) != len(soft_samples):
        raise ValueError("GEO sample titles are not unique")

    with gzip.open(matrix, "rt", encoding="utf-8", errors="strict") as handle:
        header = handle.readline().rstrip("\n\r").split("\t")
        if not header or header[0] != "gene":
            raise ValueError(f"Expected first header field 'gene', got {header[:1]!r}")
        barcodes = header[1:]
        if not barcodes:
            raise ValueError("No cell columns found")
        if len(set(barcodes)) != len(barcodes):
            raise ValueError("Cell-barcode headers are not unique")

        prefix_counts: Counter[str] = Counter()
        unparseable_barcodes: list[str] = []
        for barcode in barcodes:
            match = PREFIX.match(barcode)
            if match is None:
                unparseable_barcodes.append(barcode)
            else:
                prefix_counts[match.group("prefix")] += 1
        if unparseable_barcodes:
            raise ValueError(f"Unparseable sample prefix in barcode: {unparseable_barcodes[:5]!r}")

        gene_count = 0
        gene_ids: set[str] = set()
        duplicate_gene_ids: list[str] = []
        value_positions = audited_value_positions(len(barcodes))
        sampled_values_checked = 0
        sampled_gene_rows = 0
        all_rows_match_header = True
        for line_number, raw in enumerate(handle, start=2):
            if not raw.strip():
                continue
            gene, separator, value_text = raw.partition("\t")
            if not separator:
                raise ValueError(f"row {line_number}: missing tab after gene identifier")
            if not gene:
                raise ValueError(f"row {line_number}: empty gene identifier")
            if gene in gene_ids:
                duplicate_gene_ids.append(gene)
            else:
                gene_ids.add(gene)
            row_column_count = value_text.rstrip("\n\r").count("\t") + 1
            if row_column_count != len(barcodes):
                all_rows_match_header = False
                raise ValueError(
                    f"row {line_number}: {row_column_count} values, expected {len(barcodes)}"
                )
            gene_count += 1
            if is_audited_gene_row(gene_count):
                values = value_text.rstrip("\n\r").split("\t")
                audit_value_sample(values, value_positions, line_number)
                sampled_gene_rows += 1
                sampled_values_checked += len(value_positions)

    observed_prefixes = sorted(prefix_counts)
    missing_soft_titles = sorted(set(observed_prefixes) - set(samples_by_title))
    soft_titles_without_barcodes = sorted(set(samples_by_title) - set(observed_prefixes))
    if missing_soft_titles or soft_titles_without_barcodes:
        raise ValueError(
            "Barcode prefix/GEO title mismatch: "
            f"missing SOFT titles={missing_soft_titles}; "
            f"SOFT titles without barcodes={soft_titles_without_barcodes}"
        )

    matrix_hash = sha256(matrix)
    audit_time = datetime.now(timezone.utc).isoformat()
    matrix_status = (
        "GEO supplementary combined count TSV; GEO processing reports UMI-tools raw-count "
        "calculation followed by Seurat cell filtering"
    )
    analysis_note = (
        "Eligible for donor-level NP count aggregation only conditional on using the GEO-retained "
        "cells and reporting prior GEO cell filtering; not a reprocessed raw-read analysis"
    )

    audit_row = {
        "dataset": "GSE153066",
        "accession": "GSE153066",
        "matrix_asset": matrix.name,
        "matrix_sha256": matrix_hash,
        "matrix_bytes": matrix.stat().st_size,
        "matrix_format": "gzip-compressed dense TSV; genes in rows, cells in columns",
        "header_first_field": header[0],
        "genes": gene_count,
        "cells": len(barcodes),
        "unique_gene_identifiers": len(gene_ids),
        "duplicate_gene_identifier_count": len(duplicate_gene_ids),
        "duplicate_gene_identifier_examples": "|".join(duplicate_gene_ids[:10]),
        "unique_cell_barcodes": len(set(barcodes)),
        "barcode_prefixes": "|".join(observed_prefixes),
        "barcode_prefix_count": len(observed_prefixes),
        "prefixes_match_geo_sample_titles": True,
        "all_rows_match_header_columns": all_rows_match_header,
        "value_type_check": (
            "integer and non-negative in a deterministic sample: first and every 500th "
            f"non-empty gene row, with {len(value_positions)} evenly distributed columns per sampled row"
        ),
        "sampled_gene_rows": sampled_gene_rows,
        "sampled_columns_per_gene_row": len(value_positions),
        "sampled_matrix_values_checked": sampled_values_checked,
        "matrix_status": matrix_status,
        "analysis_note": analysis_note,
        "audited_at_utc": audit_time,
    }

    ledger_fields = [
        "dataset", "accession", "gsm", "donor_id", "library_id", "barcode_prefix",
        "cells_in_matrix", "compartment", "disease_state", "age", "sex", "sample_title",
        "tissue", "matrix_status", "barcode_to_sample_mapping", "eligible_for_count_pseudobulk",
        "eligible_for_score_level_validation", "metadata_note", "raw_characteristics",
    ]
    ledger_rows: list[dict[str, object]] = []
    for prefix in observed_prefixes:
        sample = samples_by_title[prefix]
        status = sample["status"].strip().lower()
        if status == "relatively normal":
            disease_state = "relatively normal"
        elif status == "degenerated":
            disease_state = "degenerated"
        else:
            raise ValueError(f"Unexpected GSE153066 status for {prefix}: {sample['status']!r}")
        ledger_rows.append(
            {
                "dataset": "GSE153066",
                "accession": "GSE153066",
                "gsm": sample["gsm"],
                "donor_id": prefix,
                "library_id": prefix,
                "barcode_prefix": prefix,
                "cells_in_matrix": prefix_counts[prefix],
                "compartment": "NP",
                "disease_state": disease_state,
                "age": sample["age"],
                "sex": sample["sex"],
                "sample_title": sample["sample_title"],
                "tissue": sample["tissue"],
                "matrix_status": matrix_status,
                "barcode_to_sample_mapping": "pass: exact barcode prefix equals unique GEO Sample_title",
                "eligible_for_count_pseudobulk": True,
                "eligible_for_score_level_validation": True,
                "metadata_note": (
                    "GEO provides one sample title per GSM and no separately named patient identifier; "
                    "sample title is used as the presumed donor/library key. Clinical source and age "
                    "are confounded with disease status."
                ),
                "raw_characteristics": sample["raw_characteristics"],
            }
        )

    for output in (args.matrix_audit_output, args.sample_ledger_output):
        output.parent.mkdir(parents=True, exist_ok=True)
    with args.matrix_audit_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audit_row))
        writer.writeheader()
        writer.writerow(audit_row)
    with args.sample_ledger_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ledger_fields)
        writer.writeheader()
        writer.writerows(ledger_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
