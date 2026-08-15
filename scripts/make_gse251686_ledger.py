"""Build an auditable sample ledger for the GSE251686 nested archive.

The GEO SOFT record supplies the NP source label and mild/severe sample-title
groups, while the nested-archive audit supplies the matrix member and its
dimensions.  GEO exposes no patient identifier or individual Pfirrmann grade,
so a GSM is retained only as a presumed sample/library key.  The independent
stream audit is a hard downstream-integrity gate: a malformed Matrix Market
payload is retained in the ledger for traceability but excluded from every
analysis by default.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path


GROUP_RE = re.compile(r"^(Mildly|Severely)\s+degeneration", re.IGNORECASE)
REPLICATE_RE = re.compile(r"replicate\s+(\d+)", re.IGNORECASE)
LIBRARY_RE = re.compile(r"_(NP\d+)_")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def severity_from_title(title: str) -> str:
    match = GROUP_RE.search(title)
    if not match:
        raise ValueError(f"Cannot determine mild/severe group from title: {title!r}")
    return {"mildly": "mild", "severely": "severe"}[match.group(1).lower()]


def replicate_from_title(title: str) -> str:
    match = REPLICATE_RE.search(title)
    if not match:
        raise ValueError(f"Cannot determine replicate label from title: {title!r}")
    return match.group(1)


def library_from_member(member: str) -> str:
    match = LIBRARY_RE.search(member)
    if not match:
        raise ValueError(f"Cannot determine NP library label from member: {member!r}")
    return match.group(1)


def resolve_stream_audit_path(nested_audit: Path, requested: Path | None) -> Path:
    """Resolve the stream audit without allowing a silent integrity bypass."""
    if requested is not None:
        return requested
    inferred = nested_audit.with_name("GSE251686_nested_matrix_stream_audit.csv")
    if inferred.exists():
        return inferred
    raise FileNotFoundError(
        "A stream audit is required. Pass --stream-audit or place "
        f"{inferred.name} beside the nested audit."
    )


def stream_failures(audit: dict[str, str]) -> list[str]:
    """Return failed integrity gates using the machine-readable stream audit."""
    failures: list[str] = []
    if audit.get("dimension_match") != "True":
        failures.append("dimension mismatch")
    line_count = audit.get("coordinate_lines_observed")
    header_nnz = audit.get("matrix_nnz_header")
    line_count_matches = audit.get("line_count_matches_header")
    if line_count_matches is None:
        line_count_matches = str(line_count == header_nnz)
    if line_count_matches != "True":
        failures.append("Matrix Market line count mismatch")
    valid_lines = audit.get("valid_coordinate_lines")
    header_nnz = audit.get("matrix_nnz_header")
    if valid_lines is None or header_nnz is None or valid_lines != header_nnz:
        failures.append("valid coordinate count mismatch")
    malformed = int(audit.get("malformed_coordinate_lines", "0") or 0)
    out_of_range = int(audit.get("out_of_range_coordinates", "0") or 0)
    negative = int(audit.get("negative_values", "0") or 0)
    zero = int(audit.get("zero_values", "0") or 0)
    if malformed or out_of_range or negative or zero:
        failures.append("illegal Matrix Market payload")
    if audit.get("text_integrity_pass") != "True":
        failures.append("text integrity failure")
    # The CSV stream audit does not expose outer gzip/TAR checks. Those checks
    # are audited separately and must not become false failures merely because
    # their columns are absent here.
    return failures


def stream_derived_values(audit: dict[str, str]) -> dict[str, str]:
    """Materialize stable booleans for the ledger from the stream-audit CSV."""
    valid_lines = audit.get("valid_coordinate_lines")
    header_nnz = audit.get("matrix_nnz_header")
    line_count = audit.get("coordinate_lines_observed")
    dimension = audit.get("dimension_match", "False")
    line_match = audit.get("line_count_matches_header")
    if line_match is None:
        line_match = str(line_count == header_nnz)
    valid_match = str(valid_lines is not None and valid_lines == header_nnz)
    malformed = int(audit.get("malformed_coordinate_lines", "0") or 0)
    out_of_range = int(audit.get("out_of_range_coordinates", "0") or 0)
    negative = int(audit.get("negative_values", "0") or 0)
    zero = int(audit.get("zero_values", "0") or 0)
    payload_legal = str(not (malformed or out_of_range or negative or zero) and valid_match == "True")
    return {
        "text_integrity_pass": audit.get("text_integrity_pass", "False"),
        "matrix_payload_legal": payload_legal,
        "line_count_matches_header": line_match,
        "valid_count_matches_header": valid_match,
        "malformed_lines": str(malformed),
        "nul_bytes": audit.get("nul_bytes", "0"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("soft_metadata", type=Path)
    parser.add_argument("nested_audit", type=Path)
    parser.add_argument(
        "--stream-audit",
        type=Path,
        default=None,
        help=(
            "Independent nested matrix stream audit CSV. If omitted, infer "
            "GSE251686_nested_matrix_stream_audit.csv beside nested_audit."
        ),
    )
    parser.add_argument("--ledger-output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    args = parser.parse_args()

    metadata = read_csv(args.soft_metadata)
    audit_rows = read_csv(args.nested_audit)
    audit_by_gsm = {row["gsm"]: row for row in audit_rows}
    stream_audit_path = resolve_stream_audit_path(args.nested_audit, args.stream_audit)
    stream_rows = read_csv(stream_audit_path)
    stream_by_gsm = {row["gsm"]: row for row in stream_rows}
    if len(stream_rows) != len(stream_by_gsm):
        raise ValueError("Duplicate GSM in stream audit")
    metadata_gsms = {row["gsm"] for row in metadata}
    audit_gsms = set(audit_by_gsm)
    if metadata_gsms != audit_gsms:
        raise ValueError(
            "SOFT/nested-audit GSM mismatch: "
            f"missing_in_audit={sorted(metadata_gsms - audit_gsms)}, "
            f"unexpected_in_audit={sorted(audit_gsms - metadata_gsms)}"
        )
    if len(metadata) != len(metadata_gsms):
        raise ValueError("Duplicate GSM in SOFT metadata")
    if metadata_gsms != set(stream_by_gsm):
        raise ValueError(
            "SOFT/stream-audit GSM mismatch: "
            f"missing_in_stream={sorted(metadata_gsms - set(stream_by_gsm))}, "
            f"unexpected_in_stream={sorted(set(stream_by_gsm) - metadata_gsms)}"
        )

    fields = [
        "dataset",
        "gsm",
        "sample_title",
        "presumed_donor_or_library_key",
        "archive_library_label",
        "title_replicate_label",
        "compartment",
        "severity_group",
        "severity_group_source",
        "individual_pfirrmann_grade",
        "patient_id",
        "age",
        "sex",
        "disc_level",
        "sra_accession",
        "outer_matrix_member",
        "features",
        "pre_qc_barcodes_in_supplied_matrix",
        "matrix_nnz",
        "matrix_status",
        "stream_audit_text_integrity_pass",
        "stream_audit_matrix_payload_legal",
        "stream_audit_malformed_lines",
        "stream_audit_nul_bytes",
        "eligible_for_primary_count_model",
        "eligible_for_small_sample_score_direction_check",
        "analysis_role",
        "analysis_inclusion",
        "exclusion_reason",
        "metadata_note",
    ]
    rows: list[dict[str, object]] = []
    library_labels: set[str] = set()
    for meta in sorted(metadata, key=lambda row: row["gsm"]):
        audit = audit_by_gsm[meta["gsm"]]
        stream = stream_by_gsm[meta["gsm"]]
        library_label = library_from_member(audit["outer_member"])
        if library_label in library_labels:
            raise ValueError(f"Duplicate archive library label: {library_label}")
        library_labels.add(library_label)
        severity = severity_from_title(meta["sample_title"])
        failures = stream_failures(stream)
        derived = stream_derived_values(stream)
        included = not failures
        exclusion_reason = "; ".join(failures) if failures else ""
        if included:
            analysis_role = "exploratory NP mild-versus-severe direction check only"
            analysis_inclusion = "included after stream-integrity audit"
        else:
            analysis_role = "excluded from all downstream analysis"
            analysis_inclusion = "excluded by stream-integrity audit"
        rows.append(
            {
                "dataset": "GSE251686",
                "gsm": meta["gsm"],
                "sample_title": meta["sample_title"],
                "presumed_donor_or_library_key": meta["gsm"],
                "archive_library_label": library_label,
                "title_replicate_label": replicate_from_title(meta["sample_title"]),
                "compartment": "NP",
                "severity_group": severity,
                "severity_group_source": "GEO sample title",
                "individual_pfirrmann_grade": "not exposed at GSM level in GEO/SOFT",
                "patient_id": "not exposed",
                "age": "not exposed",
                "sex": "not exposed",
                "disc_level": "not exposed",
                "sra_accession": meta["sra"],
                "outer_matrix_member": audit["outer_member"],
                "features": audit["features"],
                "pre_qc_barcodes_in_supplied_matrix": audit["barcodes"],
                "matrix_nnz": audit["matrix_nnz"],
                "matrix_status": (
                    "integer UMI-derived matrix generated with CeleScope; "
                    "outer member names record EmptyDrops; not raw sequencing reads"
                ),
                "stream_audit_text_integrity_pass": derived["text_integrity_pass"],
                "stream_audit_matrix_payload_legal": derived["matrix_payload_legal"],
                "stream_audit_malformed_lines": derived["malformed_lines"],
                "stream_audit_nul_bytes": derived["nul_bytes"],
                "eligible_for_primary_count_model": False,
                "eligible_for_small_sample_score_direction_check": included,
                "analysis_role": analysis_role,
                "analysis_inclusion": analysis_inclusion,
                "exclusion_reason": exclusion_reason,
                "metadata_note": (
                    "GSM/title is a presumed sample/library key rather than independently "
                    "verified biological-donor identity. GEO/SOFT reports GEXSCOPE/Singleron "
                    "processing and CeleScope v1.10.0; the linked paper separately calls the "
                    "platform 10X Genomics, a documented conflict."
                ),
            }
        )

    groups_all = Counter(str(row["severity_group"]) for row in rows)
    included_rows = [row for row in rows if row["analysis_inclusion"] == "included after stream-integrity audit"]
    groups_included = Counter(str(row["severity_group"]) for row in included_rows)
    if len(rows) != 6:
        raise ValueError(f"Expected six GEO sample records, found {len(rows)}")
    if not included_rows:
        raise ValueError("No GSE251686 sample passed the stream-integrity gate")

    args.ledger_output.parent.mkdir(parents=True, exist_ok=True)
    with args.ledger_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    summary = [
        {"metric": "samples_total_geo_records", "value": len(rows)},
        {"metric": "samples_included_downstream", "value": len(included_rows)},
        {"metric": "samples_excluded_stream_integrity", "value": len(rows) - len(included_rows)},
        {"metric": "mild_presumed_sample_keys_total", "value": groups_all["mild"]},
        {"metric": "severe_presumed_sample_keys_total", "value": groups_all["severe"]},
        {"metric": "mild_presumed_sample_keys_included", "value": groups_included["mild"]},
        {"metric": "severe_presumed_sample_keys_included", "value": groups_included["severe"]},
        {
            "metric": "excluded_gsms",
            "value": ",".join(str(row["gsm"]) for row in rows if row["analysis_inclusion"] != "included after stream-integrity audit"),
        },
        {
            "metric": "pre_qc_barcodes_all_geo_records",
            "value": sum(int(row["pre_qc_barcodes_in_supplied_matrix"]) for row in rows),
        },
        {
            "metric": "pre_qc_barcodes_included_downstream",
            "value": sum(int(row["pre_qc_barcodes_in_supplied_matrix"]) for row in included_rows),
        },
        {
            "metric": "patient_identity_status",
            "value": "not exposed; GSM/title used as presumed sample/library key",
        },
        {
            "metric": "individual_grade_status",
            "value": "not exposed at GSM level in GEO/SOFT",
        },
        {
            "metric": "platform_status",
            "value": "GEO/SOFT: GEXSCOPE/Singleron and CeleScope; linked paper results: 10X Genomics (conflict documented)",
        },
        {
            "metric": "permitted_role",
            "value": "incomplete, non-balanced (mild n=2 versus severe n=3) exploratory NP severity direction check only",
        },
    ]
    with args.summary_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"])
        writer.writeheader()
        writer.writerows(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
