"""Build conservative GSE244889 provenance and processed-asset ledgers.

This script does not load single-cell expression values or perform biological
analysis.  It reads the GEO SOFT metadata and bulk FPKM table only to preserve
the sample/asset boundary required before later analysis.  The 10x Matrix
Market structure audit remains delegated to ``audit_10x_tar.py``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
import tarfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244889"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def title_prefix(title: str) -> str:
    return title.split(",", maxsplit=1)[0].strip()


def title_age_sex(title: str) -> tuple[str, str]:
    """Return only explicitly title-encoded age/sex values, if present."""
    suffix = re.search(r"(?<!\d)(?P<age>\d{1,3})(?P<sex>[FM])(?=\s*,|\s*$)", title)
    if suffix:
        return suffix.group("age"), suffix.group("sex")
    prefix = re.search(r"-(?P<sex>[FM])(?P<age>\d{1,3})(?=\s*,|\s*$)", title)
    return (prefix.group("age"), prefix.group("sex")) if prefix else ("", "")


def normalized_grade(value: str) -> str:
    match = re.search(r"(?:grade\s*)?(\d+)", value, flags=re.IGNORECASE)
    return match.group(1) if match else ""


def severity_group(status: str) -> str:
    normalized = status.strip().upper()
    if normalized == "MDD":
        return "mild"
    if normalized == "SDD":
        return "severe"
    return "unresolved"


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def audit_fpkm(path: Path) -> dict[str, object]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        genes: list[str] = []
        row_width_mismatches = 0
        non_numeric_values = 0
        empty_values = 0
        negative_values = 0
        minimum = float("inf")
        maximum = float("-inf")
        rows = 0
        for row in reader:
            rows += 1
            if not row:
                row_width_mismatches += 1
                continue
            genes.append(row[0])
            if len(row) != len(header):
                row_width_mismatches += 1
                continue
            for value in row[1:]:
                if value == "":
                    empty_values += 1
                    continue
                try:
                    numeric = float(value)
                except ValueError:
                    non_numeric_values += 1
                    continue
                negative_values += numeric < 0
                minimum = min(minimum, numeric)
                maximum = max(maximum, numeric)
    gene_counts = Counter(genes)
    duplicate_ids = sorted(gene for gene, count in gene_counts.items() if count > 1)
    return {
        "header_first_field": header[0] if header else "",
        "sample_columns": len(header) - 1,
        "bulk_columns": "|".join(header[1:]),
        "gene_rows": rows,
        "unique_gene_identifiers": len(set(genes)),
        "duplicate_gene_identifier_count": rows - len(set(genes)),
        "duplicate_gene_identifiers": "|".join(duplicate_ids),
        "row_width_mismatches": row_width_mismatches,
        "non_numeric_values": non_numeric_values,
        "empty_values": empty_values,
        "negative_values": negative_values,
        "minimum_numeric_value": minimum if minimum != float("inf") else "",
        "maximum_numeric_value": maximum if maximum != float("-inf") else "",
        "matrix_status": "GEO-supplied processed FPKM table; not a raw-count matrix",
        "eligible_for_raw_count_pseudobulk": False,
        "eligible_for_single_cell_aggregation": False,
    }


def audit_filelist(filelist: Path, raw_tar: Path) -> dict[str, object]:
    with filelist.open(newline="", encoding="utf-8", errors="replace") as handle:
        listed = list(csv.DictReader(handle, delimiter="\t"))
    archive_rows = [row for row in listed if row.get("#Archive/File") == "Archive"]
    member_rows = [row for row in listed if row.get("#Archive/File") == "File"]
    if len(archive_rows) != 1:
        raise ValueError(f"Expected exactly one archive row in {filelist.name}")
    archive_row = archive_rows[0]
    with tarfile.open(raw_tar, "r") as handle:
        actual_members = {member.name: member.size for member in handle.getmembers() if member.isfile()}
    listed_members = {row["Name"]: int(row["Size"]) for row in member_rows}
    missing_from_filelist = sorted(set(actual_members) - set(listed_members))
    missing_from_tar = sorted(set(listed_members) - set(actual_members))
    size_mismatches = sorted(
        name for name in set(actual_members) & set(listed_members)
        if actual_members[name] != listed_members[name]
    )
    return {
        "filelist_asset": filelist.name,
        "archive_name_matches": archive_row["Name"] == raw_tar.name,
        "archive_bytes_listed": int(archive_row["Size"]),
        "archive_bytes_actual": raw_tar.stat().st_size,
        "archive_size_matches": int(archive_row["Size"]) == raw_tar.stat().st_size,
        "filelist_member_count": len(listed_members),
        "tar_file_member_count": len(actual_members),
        "missing_from_filelist_count": len(missing_from_filelist),
        "missing_from_tar_count": len(missing_from_tar),
        "member_size_mismatch_count": len(size_mismatches),
        "missing_from_filelist": "|".join(missing_from_filelist),
        "missing_from_tar": "|".join(missing_from_tar),
        "member_size_mismatches": "|".join(size_mismatches),
        "filelist_tar_crosscheck_pass": not (missing_from_filelist or missing_from_tar or size_mismatches)
        and archive_row["Name"] == raw_tar.name
        and int(archive_row["Size"]) == raw_tar.stat().st_size,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("soft_metadata", type=Path)
    parser.add_argument("soft", type=Path)
    parser.add_argument("raw_tar", type=Path)
    parser.add_argument("fpkm", type=Path)
    parser.add_argument("filelist", type=Path)
    parser.add_argument("--asset-output", type=Path, required=True)
    parser.add_argument("--scrna-ledger-output", type=Path, required=True)
    parser.add_argument("--bulk-ledger-output", type=Path, required=True)
    parser.add_argument("--fpkm-audit-output", type=Path, required=True)
    parser.add_argument("--filelist-audit-output", type=Path, required=True)
    args = parser.parse_args()

    with args.soft_metadata.open(newline="", encoding="utf-8") as handle:
        metadata = list(csv.DictReader(handle))

    scrna_rows: list[dict[str, object]] = []
    bulk_rows: list[dict[str, object]] = []
    for row in metadata:
        title = row["sample_title"]
        title_lower = title.lower()
        age, sex = title_age_sex(title)
        common = {
            "dataset": "GSE244889",
            "gsm": row["gsm"],
            "sample_title": title,
            "tissue": row["tissue"],
            "source_name": row["source_name"],
            "age_reported_in_title": age,
            "sex_reported_in_title": sex,
            "degeneration_grade_reported": normalized_grade(row["degeneration_grade"] or title),
            "geo_status": row["disease_state"],
            "severity_group": severity_group(row["disease_state"]),
            "grouping_basis": "GEO Sample_status MDD= mild; SDD= severe",
        }
        if "scrnaseq" in title_lower:
            key = title_prefix(title)
            scrna_rows.append(
                {
                    **common,
                    "assay": "scRNA-seq",
                    "presumed_donor_library_key": key,
                    "donor_id": key,
                    "donor_key_basis": "unique scRNA sample-title prefix; GEO SOFT does not expose a separate patient ID",
                    "compartment": "NP",
                    "disease_state": row["disease_state"],
                    "matrix_status": "raw 10x Cell Ranger v6.1.1 MTX triple reported by GEO",
                    "eligible_for_raw_count_pseudobulk_after_qc": True,
                    "metadata_limit": "Sample-title key is a presumed donor/library key until independently verified; small mild/severe groups and title-encoded age require reporting.",
                }
            )
        elif "rnaseq" in title_lower:
            key = title_prefix(title).split("-", maxsplit=1)[0]
            bulk_rows.append(
                {
                    **common,
                    "assay": "bulk RNA-seq",
                    "presumed_donor_key": key,
                    "donor_key_basis": "S11-S16 token in GEO sample title; no separately named patient ID",
                    "fpkm_column_match_basis": "not yet assigned; resolved against FPKM header by shared S11-S16 token below",
                    "matrix_status": "GEO-supplied FPKM table, processed expression values",
                    "eligible_for_raw_count_pseudobulk": False,
                    "metadata_limit": "FPKM must remain separate from 10x data; it is not a raw-count or single-cell asset.",
                }
            )
        else:
            raise ValueError(f"Could not classify assay from sample title: {title!r}")

    if len(scrna_rows) != 7 or len(bulk_rows) != 6:
        raise ValueError(f"Expected 7 scRNA and 6 bulk SOFT samples, found {len(scrna_rows)} and {len(bulk_rows)}")

    fpkm = audit_fpkm(args.fpkm)
    bulk_columns = fpkm["bulk_columns"].split("|") if fpkm["bulk_columns"] else []
    bulk_by_key = {row["presumed_donor_key"]: row for row in bulk_rows}
    if len(bulk_columns) != len(bulk_rows):
        raise ValueError("FPKM columns do not match the expected six bulk SOFT samples")
    for column in bulk_columns:
        token = column.split("_", maxsplit=1)[0]
        if token not in bulk_by_key:
            raise ValueError(f"No bulk GEO title matches FPKM column {column!r}")
        bulk_by_key[token]["fpkm_column"] = column
        bulk_by_key[token]["fpkm_column_match_basis"] = (
            f"shared {token} token between FPKM header and GEO bulk sample title; "
            "no explicit bulk-column-to-GSM mapping is provided in the table header"
        )

    asset_specs = (
        (args.raw_tar, "GSE244889_RAW.tar", f"{GEO_BASE}/suppl/GSE244889_RAW.tar", "raw archive; 7 scRNA 10x triples"),
        (args.soft, "GSE244889_family.soft.gz", f"{GEO_BASE}/soft/GSE244889_family.soft.gz", "GEO SOFT metadata"),
        (args.fpkm, "GSE244889_FPKM.txt.gz", f"{GEO_BASE}/suppl/GSE244889_FPKM.txt.gz", "processed bulk FPKM table"),
        (args.filelist, "filelist.txt", f"{GEO_BASE}/suppl/filelist.txt", "GEO raw-archive member listing"),
    )
    inspected_at = datetime.now(timezone.utc).isoformat()
    assets = [
        {
            "dataset": "GSE244889",
            "asset": name,
            "url": url,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "asset_role": role,
            "audited_at_utc": inspected_at,
        }
        for path, name, url, role in asset_specs
    ]

    write_csv(
        args.asset_output,
        ["dataset", "asset", "url", "bytes", "sha256", "asset_role", "audited_at_utc"],
        assets,
    )
    write_csv(
        args.scrna_ledger_output,
        [
            "dataset", "gsm", "assay", "sample_title", "presumed_donor_library_key", "donor_id",
            "donor_key_basis", "tissue", "source_name", "age_reported_in_title", "sex_reported_in_title",
            "degeneration_grade_reported", "geo_status", "disease_state", "severity_group", "grouping_basis",
            "compartment", "matrix_status",
            "eligible_for_raw_count_pseudobulk_after_qc", "metadata_limit",
        ],
        scrna_rows,
    )
    write_csv(
        args.bulk_ledger_output,
        [
            "dataset", "gsm", "assay", "sample_title", "presumed_donor_key", "donor_key_basis", "fpkm_column",
            "fpkm_column_match_basis", "tissue", "source_name", "age_reported_in_title", "sex_reported_in_title",
            "degeneration_grade_reported", "geo_status", "severity_group", "grouping_basis", "matrix_status",
            "eligible_for_raw_count_pseudobulk", "metadata_limit",
        ],
        bulk_rows,
    )
    write_csv(
        args.fpkm_audit_output,
        list(fpkm),
        [fpkm],
    )
    filelist_audit = audit_filelist(args.filelist, args.raw_tar)
    write_csv(
        args.filelist_audit_output,
        list(filelist_audit),
        [filelist_audit],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
