"""Aggregate source-restricted 10x cells and score locked gene modules.

The matrix is streamed from a TAR archive.  Cells are selected from the
annotation ledger, but inference remains at donor/library level.  Counts are
aggregated by gene symbol explicitly; duplicate feature rows are retained in
the audit and summed for module scoring.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import tarfile
from collections import defaultdict
from pathlib import Path

from qc_10x_archives import gsm_members, member_map, read_barcodes, read_features


GSM_RE = re.compile(r"^(GSM\d+)_")
REQUIRED_LEDGER_COLUMNS = {"gsm", "dataset", "compartment", "disease_state"}


def read_module_config(path: Path) -> tuple[dict[str, list[str]], str]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    modules: dict[str, list[str]] = {}
    for module in payload["modules"]:
        module_id = str(module["module_id"])
        genes = [str(gene).strip().upper() for gene in module["genes"]]
        if not genes or any(not gene for gene in genes):
            raise ValueError(f"Module {module_id} has an empty gene symbol")
        if len(genes) != len(set(genes)):
            raise ValueError(f"Module {module_id} contains duplicate gene symbols")
        if module_id in modules:
            raise ValueError(f"Duplicate module_id in module config: {module_id}")
        modules[module_id] = genes
    return modules, str(payload.get("score_direction", ""))


def _parse_bool(value: str, field: str, key: tuple[str, str]) -> bool:
    normalized = value.strip().lower()
    if normalized not in {"true", "false"}:
        raise ValueError(f"{field} for {key} must be true or false, got {value!r}")
    return normalized == "true"


def read_qc_passes(path: Path) -> dict[tuple[str, str], bool]:
    """Read and validate the complete per-cell QC ledger."""
    result: dict[tuple[str, str], bool] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gsm", "barcode", "qc_pass"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"QC file missing required columns: {sorted(missing)}")
        for row in reader:
            key = (row["gsm"].strip(), row["barcode"].strip())
            if not key[0] or not key[1]:
                raise ValueError("QC file contains an empty GSM or barcode")
            if key in result:
                raise ValueError(f"Duplicate QC row for {key}")
            result[key] = _parse_bool(row["qc_pass"], "qc_pass", key)
    if not result:
        raise ValueError(f"QC file is empty: {path}")
    return result


def read_inclusion(path: Path) -> dict[tuple[str, str], bool]:
    """Read and validate the annotation rows (one row per QC-passing cell)."""
    result: dict[tuple[str, str], bool] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"gsm", "barcode", "compartment_pseudobulk_include"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Annotation file missing required columns: {sorted(missing)}")
        for row in reader:
            key = (row["gsm"].strip(), row["barcode"].strip())
            if not key[0] or not key[1]:
                raise ValueError("Annotation file contains an empty GSM or barcode")
            if key in result:
                raise ValueError(f"Duplicate annotation row for {key}")
            result[key] = _parse_bool(row["compartment_pseudobulk_include"], "compartment_pseudobulk_include", key)
    if not result:
        raise ValueError(f"Annotation file is empty: {path}")
    return result


def archive_gsms(members: dict[str, tarfile.TarInfo]) -> set[str]:
    return {match.group(1) for name in members if (match := GSM_RE.match(name))}


def validate_ledger(ledger: dict[str, dict[str, str]], path: Path) -> None:
    """Reject ambiguous or incomplete sample metadata before matrix ingestion."""
    if not ledger:
        raise ValueError(f"Empty ledger: {path}")
    for gsm, row in ledger.items():
        missing = [column for column in REQUIRED_LEDGER_COLUMNS if not row.get(column, "").strip()]
        if missing:
            raise ValueError(f"{gsm}: ledger fields missing or empty: {missing}")
        if not (row.get("donor_id", "").strip() or row.get("patient_id", "").strip()):
            raise ValueError(f"{gsm}: ledger has no donor_id or patient_id")


def validate_cell_ledgers(
    archive_members: dict[str, tarfile.TarInfo],
    tar: tarfile.TarFile,
    ledger: dict[str, dict[str, str]],
    qc_passes: dict[tuple[str, str], bool],
    inclusion: dict[tuple[str, str], bool],
) -> dict[str, list[str]]:
    """Check one-to-one GSM/barcode coverage across archive, QC, and annotation."""
    archive_set = archive_gsms(archive_members)
    if not archive_set:
        raise ValueError("Archive contains no GSM-prefixed files")
    ledger_set = set(ledger)
    # A shared discovery ledger may intentionally contain both child-series;
    # the active archive must still be fully represented.  Extra ledger rows
    # are permitted only because they are not part of this archive and are
    # never traversed below.
    if not archive_set <= ledger_set:
        raise ValueError(
            "Archive/ledger GSM mismatch: "
            f"archive_only={sorted(archive_set - ledger_set)}, "
            f"ledger_only={sorted(ledger_set - archive_set)}"
        )
    if ledger_set - archive_set:
        # Keep this distinction explicit in the error-free path: a superset
        # ledger is an accepted shared ledger, not evidence that extra samples
        # were analyzed.
        pass

    expected_by_gsm: dict[str, list[str]] = {}
    for gsm in sorted(archive_set):
        parts = gsm_members(archive_members, gsm)
        barcodes = read_barcodes(parts["barcodes"], tar)
        if len(barcodes) != len(set(barcodes)):
            raise ValueError(f"{gsm}: archive barcodes are not unique")
        expected_by_gsm[gsm] = barcodes

    expected_qc = {(gsm, barcode) for gsm, barcodes in expected_by_gsm.items() for barcode in barcodes}
    qc_keys = set(qc_passes)
    if qc_keys != expected_qc:
        raise ValueError(
            "QC/archive barcode mismatch: "
            f"qc_only={len(qc_keys - expected_qc)}, archive_only={len(expected_qc - qc_keys)}"
        )
    expected_annotation = {key for key, passed in qc_passes.items() if passed}
    annotation_keys = set(inclusion)
    if annotation_keys != expected_annotation:
        raise ValueError(
            "Annotation/QC-pass barcode mismatch: "
            f"annotation_only={len(annotation_keys - expected_annotation)}, "
            f"qc_pass_only={len(expected_annotation - annotation_keys)}"
        )
    return expected_by_gsm


def strict_load_ledger(path: Path) -> dict[str, dict[str, str]]:
    """Load the sample ledger without silently overwriting duplicate GSM rows."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = set(reader.fieldnames or [])
        missing = REQUIRED_LEDGER_COLUMNS - fieldnames
        if missing:
            raise ValueError(f"Ledger missing required columns: {sorted(missing)}")
        rows: dict[str, dict[str, str]] = {}
        for row in reader:
            gsm = row.get("gsm", "").strip()
            if not gsm:
                raise ValueError("Ledger contains an empty GSM")
            if gsm in rows:
                raise ValueError(f"Ledger contains duplicate GSM row: {gsm}")
            rows[gsm] = row
    validate_ledger(rows, path)
    return rows


def donor_group(row: dict[str, str], gsm: str) -> tuple[str, str, str, str]:
    donor = row.get("donor_id", "").strip() or row.get("patient_id", "").strip() or gsm
    return (row.get("dataset", "").strip(), donor, row.get("compartment", "").strip(), row.get("disease_state", "").strip())


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("ledger", type=Path)
    parser.add_argument("annotation", type=Path)
    parser.add_argument("module_config", type=Path)
    parser.add_argument("--qc-cells", type=Path, required=True, help="gzipped TSV emitted by qc_10x_archives.py")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-mapped-fraction", type=float, default=0.80)
    parser.add_argument(
        "--min-retained-cells", type=int, choices=(20, 30, 50), default=30,
        help="minimum source-restricted cells per donor/library (primary=30; sensitivity=20 or 50)",
    )
    args = parser.parse_args()
    if not 0 < args.min_mapped_fraction <= 1:
        raise ValueError("--min-mapped-fraction must be in (0, 1]")
    modules, score_direction = read_module_config(args.module_config)
    ledger = strict_load_ledger(args.ledger)
    qc_passes = read_qc_passes(args.qc_cells)
    inclusion = read_inclusion(args.annotation)

    score_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    mapping_rows: list[dict[str, object]] = []
    library_rows: list[dict[str, object]] = []

    with tarfile.open(args.archive, "r") as tar:
        members = member_map(tar)
        barcodes_by_gsm = validate_cell_ledgers(members, tar, ledger, qc_passes, inclusion)
        gsms = sorted(barcodes_by_gsm)
        seen_donor_compartments: dict[tuple[str, str, str], str] = {}
        for gsm in gsms:
            metadata = ledger[gsm]
            parts = gsm_members(members, gsm)
            barcodes = barcodes_by_gsm[gsm]
            feature_names, _ = read_features(parts["features"], tar)
            include_indices = {i for i, barcode in enumerate(barcodes) if inclusion.get((gsm, barcode), False)}
            if not include_indices:
                raise ValueError(f"{gsm}: no source-restricted cells in annotation ledger")
            if len(include_indices) < args.min_retained_cells:
                raise ValueError(
                    f"{gsm}: {len(include_indices)} source-restricted cells is below the "
                    f"required {args.min_retained_cells}-cell threshold"
                )
            group = donor_group(metadata, gsm)
            donor_compartment = group[:3]
            if donor_compartment in seen_donor_compartments:
                previous = seen_donor_compartments[donor_compartment]
                raise ValueError(
                    f"{gsm}: donor/compartment {donor_compartment} already represented by {previous}; "
                    "multiple libraries per donor/compartment must be aggregated before scoring"
                )
            seen_donor_compartments[donor_compartment] = gsm
            normalized_features = [name.strip().upper() for name in feature_names]
            # Reset all feature-to-module state for every library.  Keeping
            # this inside the GSM loop prevents duplicate-row audit counts
            # from leaking across libraries.
            module_gene_to_ids: defaultdict[str, dict[str, list[int]]] = defaultdict(
                lambda: defaultdict(list)
            )
            for module_id, genes in modules.items():
                for gene in genes:
                    module_gene_to_ids[module_id][gene] = []
            feature_to_modules: defaultdict[int, list[tuple[str, str]]] = defaultdict(list)
            measured_by_module: defaultdict[str, set[str]] = defaultdict(set)
            for feature_index, symbol in enumerate(normalized_features):
                for module_id, gene_map in module_gene_to_ids.items():
                    if symbol in gene_map:
                        gene_map[symbol].append(feature_index)
                        feature_to_modules[feature_index].append((module_id, symbol))
                        measured_by_module[module_id].add(symbol)
            module_counts: dict[str, dict[str, int]] = {
                module_id: {gene: 0 for gene in genes} for module_id, genes in modules.items()
            }
            total_umi = 0
            stream = tar.extractfile(parts["matrix"])
            if stream is None:
                raise OSError(parts["matrix"].name)
            with stream, gzip.open(stream, "rt", encoding="utf-8", errors="replace") as handle:
                first = handle.readline().strip()
                if not first.startswith("%%MatrixMarket matrix coordinate"):
                    raise ValueError(f"Unexpected Matrix Market header: {parts['matrix'].name}")
                line = handle.readline().strip()
                while line.startswith("%"):
                    line = handle.readline().strip()
                rows, columns, _ = (int(value) for value in line.split())
                if (rows, columns) != (len(feature_names), len(barcodes)):
                    raise ValueError(f"Dimension mismatch in {parts['matrix'].name}")
                declared_nnz: int | None = None
                entries_read = 0
                for raw in handle:
                    if not raw.strip() or raw.startswith("%"):
                        continue
                    fields = raw.split()
                    if len(fields) < 3:
                        raise ValueError(f"Malformed Matrix Market entry in {parts['matrix'].name}: {raw!r}")
                    feature_index, cell_index, value = (int(value) for value in fields[:3])
                    cell_index -= 1
                    feature_index -= 1
                    if not (0 <= feature_index < len(feature_names) and 0 <= cell_index < len(barcodes)):
                        raise ValueError(
                            f"Out-of-range coordinate in {parts['matrix'].name}: "
                            f"{feature_index + 1}, {cell_index + 1}"
                        )
                    if value < 0:
                        raise ValueError(f"Negative count in {parts['matrix'].name}: {value}")
                    entries_read += 1
                    if cell_index not in include_indices or value <= 0:
                        continue
                    total_umi += value
                    for module_id, symbol in feature_to_modules.get(feature_index, []):
                        module_counts[module_id][symbol] += value
                # Matrix Market's third dimension is the number of coordinate
                # records, including zero-valued records.  Validate it so a
                # truncated or silently malformed stream cannot produce a
                # plausible-looking pseudobulk.
                declared_nnz = int(line.split()[2])
                if entries_read != declared_nnz:
                    raise ValueError(
                        f"Matrix entry count mismatch in {parts['matrix'].name}: "
                        f"declared {declared_nnz}, read {entries_read}"
                    )

            _, donor_id, compartment, disease_state = group
            library_rows.append(
                {
                    "dataset": metadata.get("dataset", ""), "gsm": gsm, "donor_id": donor_id,
                    "compartment": compartment, "disease_state": disease_state,
                    "included_cells": len(include_indices), "total_umi_included_cells": total_umi,
                    "annotation_source": str(args.annotation),
                }
            )
            for module_id, genes in modules.items():
                mapped = sorted(measured_by_module[module_id])
                fraction = len(mapped) / len(genes) if genes else 0.0
                mapping_rows.append(
                    {
                        "dataset": metadata.get("dataset", ""), "gsm": gsm, "donor_id": donor_id,
                        "compartment": compartment, "module_id": module_id,
                        "configured_genes": len(genes), "mapped_genes": len(mapped),
                        "mapped_fraction": f"{fraction:.6f}",
                        "mapping_pass": fraction >= args.min_mapped_fraction,
                        "mapped_gene_symbols": ";".join(mapped),
                        "duplicate_feature_rows": sum(max(0, len(module_gene_to_ids[module_id][gene]) - 1) for gene in mapped),
                    }
                )
                for gene in genes:
                    gene_rows.append(
                        {
                            "dataset": metadata.get("dataset", ""), "gsm": gsm, "donor_id": donor_id,
                            "compartment": compartment, "disease_state": disease_state,
                            "module_id": module_id, "gene_symbol": gene,
                            "feature_rows_mapped": len(module_gene_to_ids[module_id][gene]),
                            "pseudobulk_count": module_counts[module_id][gene],
                            "total_umi_included_cells": total_umi,
                        }
                    )
                if fraction >= args.min_mapped_fraction and total_umi > 0:
                    gene_scores = [math.log1p(1_000_000.0 * module_counts[module_id][gene] / total_umi) for gene in mapped]
                    score = sum(gene_scores) / len(gene_scores) if gene_scores else float("nan")
                    status = "score_available"
                else:
                    score = float("nan")
                    status = "mapping_below_minimum_or_zero_library_umi"
                score_rows.append(
                    {
                        "dataset": metadata.get("dataset", ""), "gsm": gsm, "donor_id": donor_id,
                        "compartment": compartment, "disease_state": disease_state,
                        "module_id": module_id, "module_score_log1p_cpm": "" if math.isnan(score) else f"{score:.8f}",
                        "score_status": status, "mapped_fraction": f"{fraction:.6f}",
                        "included_cells": len(include_indices), "total_umi_included_cells": total_umi,
                    }
                )

    prefix = args.archive.stem
    write_csv(args.output_dir / f"{prefix}_module_scores.csv", score_rows, list(score_rows[0]) if score_rows else ["module_id"])
    write_csv(args.output_dir / f"{prefix}_module_gene_pseudobulk.csv", gene_rows, list(gene_rows[0]) if gene_rows else ["gene_symbol"])
    write_csv(args.output_dir / f"{prefix}_module_mapping_audit.csv", mapping_rows, list(mapping_rows[0]) if mapping_rows else ["module_id"])
    write_csv(args.output_dir / f"{prefix}_library_pseudobulk_ledger.csv", library_rows, list(library_rows[0]) if library_rows else ["gsm"])
    write_csv(
        args.output_dir / f"{prefix}_module_score_parameters.csv",
        [{
            "module_config": str(args.module_config), "annotation": str(args.annotation),
            "qc_cells": str(args.qc_cells), "min_retained_cells": args.min_retained_cells,
            "min_mapped_fraction": args.min_mapped_fraction,
            "formula": "mean over mapped genes of log1p(1e6 * pseudobulk_gene_count / total_UMI_included_cells)",
            "score_direction": score_direction,
            "duplicate_feature_rule": "sum all feature rows sharing the normalized uppercase gene symbol; record row count",
            "inference_unit": "donor/library; cells are nested observations",
        }],
        ["module_config", "annotation", "qc_cells", "min_retained_cells", "min_mapped_fraction", "formula", "score_direction", "duplicate_feature_rule", "inference_unit"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
