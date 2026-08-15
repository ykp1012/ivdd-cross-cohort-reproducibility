"""Audit GSE165722 matrices without extracting or altering the raw archive.

GEO calls these integer-valued supplementary matrices "normalized counts".
This script therefore produces descriptive, donor-level QC only.  It never
labels the values raw UMIs, creates pseudobulk counts, or fits count models.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import tarfile
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Sample:
    gsm: str
    sample: str
    geo_grade: int
    source_grade: int
    severity_group: str


SAMPLES = (
    # GEO labels I-IV conflict with Tu et al. PMID 34825784 Table 1 (II-V).
    # Preserve both and use source_grade only for the pre-specified grouping.
    Sample("GSM5048708", "Sample1", 1, 2, "mild"),
    Sample("GSM5048709", "Sample2", 1, 2, "mild"),
    Sample("GSM5048710", "Sample3", 2, 3, "mild"),
    Sample("GSM5048711", "Sample4", 2, 3, "mild"),
    Sample("GSM5048712", "Sample5", 3, 4, "severe"),
    Sample("GSM5048713", "Sample6", 3, 4, "severe"),
    Sample("GSM5048714", "Sample7", 4, 5, "severe"),
    Sample("GSM5048715", "Sample8", 4, 5, "severe"),
)


def read_head(handle, n: int = 6) -> list[str]:
    return [handle.readline().rstrip("\n") for _ in range(n)]


def looks_integer(values: list[str]) -> bool:
    try:
        return all(float(value).is_integer() for value in values if value)
    except ValueError:
        return False


def count_rows_and_cells(member: object) -> tuple[int, int, bool]:
    """Stream a matrix once to establish dimensions and sampled value type."""
    with gzip.open(member, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cells = max(0, len(header) - 1)
        genes = 0
        sampled_values: list[str] = []
        for line in handle:
            if not line.strip():
                continue
            genes += 1
            if genes <= 1_000:
                sampled_values.extend(line.rstrip("\n").split("\t")[1:11])
        return genes, cells, looks_integer(sampled_values)


def audit_cell_mapping(member: object, expected_ids: list[str]) -> tuple[int, bool, bool]:
    """Verify the complete CellName-to-CellIndex mapping against matrix IDs."""
    with gzip.open(member, "rt", encoding="utf-8", errors="replace") as handle:
        next(handle, None)
        pairs = [line.rstrip("\n").split("\t", maxsplit=1) for line in handle if line.strip()]
    indices = [pair[1] for pair in pairs if len(pair) == 2]
    return len(indices), indices == expected_ids, len(indices) == len(set(indices))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    archive = args.archive.resolve()
    if not archive.is_file():
        raise FileNotFoundError(archive)

    rows: list[dict[str, object]] = []
    with tarfile.open(archive, "r") as tar:
        for sample in SAMPLES:
            count_name = f"{sample.gsm}_{sample.sample}.counts.tsv.gz"
            cell_name = f"{sample.gsm}_{sample.sample}.cellname.txt.gz"
            count_member = tar.extractfile(count_name)
            cell_member = tar.extractfile(cell_name)
            if count_member is None or cell_member is None:
                raise FileNotFoundError(f"Missing archive member for {sample.gsm}")
            with gzip.open(count_member, "rt", encoding="utf-8", errors="replace") as counts:
                count_head = read_head(counts)
            cell_member = tar.extractfile(cell_name)
            if cell_member is None:
                raise FileNotFoundError(cell_name)
            with gzip.open(cell_member, "rt", encoding="utf-8", errors="replace") as cells:
                cell_file_header = cells.readline().rstrip("\n").split("\t")
                cell_head = read_head(cells)

            header = count_head[0].split("\t")
            first_data = count_head[1].split("\t") if len(count_head) > 1 else []
            cell_mapping = [line.split("\t", maxsplit=1) for line in cell_head if line.strip()]
            barcodes = [row[0] for row in cell_mapping if len(row) == 2]
            matrix_ids = [row[1] for row in cell_mapping if len(row) == 2]
            numeric_probe = first_data[1: min(len(first_data), 11)]
            count_member = tar.extractfile(count_name)
            cell_member = tar.extractfile(cell_name)
            if count_member is None or cell_member is None:
                raise FileNotFoundError(f"Missing archive member for {sample.gsm}")
            matrix_genes, matrix_cells, sampled_values_integer_like = count_rows_and_cells(count_member)
            mapping_rows, full_mapping_matches, full_cell_indices_unique = audit_cell_mapping(cell_member, header[1:])
            rows.append(
                {
                    "dataset": "GSE165722",
                    "gsm": sample.gsm,
                    "sample": sample.sample,
                    "donor_id": sample.sample,
                    "compartment": "NP",
                    "geo_reported_grade": sample.geo_grade,
                    "source_publication_grade": sample.source_grade,
                    "source_publication_severity_group": sample.severity_group,
                    "clinical_grade_authority": "Tu et al. PMID 34825784 Table 1",
                    "matrix_status": "integer count-like; GEO supplementary format states normalized counts",
                    "permitted_analysis": "descriptive QC and donor-level score direction only; no raw-count pseudobulk model",
                    "count_header_columns": len(header),
                    "first_data_columns": len(first_data),
                    "matrix_genes": matrix_genes,
                    "matrix_cells": matrix_cells,
                    "header_first_field": header[0] if header else "",
                    "first_feature": first_data[0] if first_data else "",
                    "matrix_cell_ids_in_header": max(0, len(header) - 1),
                    "cell_name_file_header": "|".join(cell_file_header),
                    "barcode_head": "|".join(barcodes),
                    "cell_index_head": "|".join(matrix_ids),
                    "matrix_cell_id_head": "|".join(header[1: 1 + len(matrix_ids)]),
                    "head_cell_indices_match": header[1: 1 + len(matrix_ids)] == matrix_ids,
                    "cell_mapping_rows_inspected": len(matrix_ids),
                    "cell_indices_unique_in_head": len(set(matrix_ids)) == len(matrix_ids),
                    "full_cell_mapping_rows": mapping_rows,
                    "full_cell_mapping_matches_matrix_header": full_mapping_matches,
                    "full_cell_indices_unique": full_cell_indices_unique,
                    "first_values_integer_like": looks_integer(numeric_probe),
                    "sampled_first_1000_rows_integer_like": sampled_values_integer_like,
                    "first_values": "|".join(numeric_probe),
                    "archive_member_bytes": tar.getmember(count_name).size,
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
