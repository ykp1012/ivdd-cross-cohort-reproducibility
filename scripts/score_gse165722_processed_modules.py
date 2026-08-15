"""Score locked modules from GSE165722 GEO-supplied processed matrices.

GEO describes these integer-like matrices as normalized counts.  This script
therefore creates per-sample program scores only: it never calls values UMIs,
does not create count pseudobulks, and does not treat individual cells as
independent observations.  Each matrix is streamed directly from the original
TAR archive without extraction or modification.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np


REQUIRED_LEDGER_COLUMNS = {
    "dataset",
    "gsm",
    "sample",
    "donor_id",
    "compartment",
    "source_publication_severity_group",
    "cells_in_supplied_matrix",
    "matrix_status",
    "eligible_for_score_level_validation",
}


def read_modules(path: Path) -> tuple[dict[str, list[str]], str]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    modules: dict[str, list[str]] = {}
    for item in payload["modules"]:
        module_id = str(item["module_id"]).strip()
        genes = [str(gene).strip().upper() for gene in item["genes"]]
        if not module_id or not genes or any(not gene for gene in genes):
            raise ValueError(f"Invalid module definition: {item!r}")
        if len(genes) != len(set(genes)):
            raise ValueError(f"Module {module_id} contains duplicate genes")
        if module_id in modules:
            raise ValueError(f"Duplicate module ID: {module_id}")
        modules[module_id] = genes
    if not modules:
        raise ValueError(f"No modules found in {path}")
    return modules, str(payload.get("score_direction", ""))


def parse_bool(value: str, field: str) -> bool:
    normalized = value.strip().lower()
    if normalized not in {"true", "false"}:
        raise ValueError(f"{field} must be true or false, got {value!r}")
    return normalized == "true"


def read_ledger(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        raise ValueError(f"Empty ledger: {path}")
    missing = REQUIRED_LEDGER_COLUMNS - set(rows[0])
    if missing:
        raise ValueError(f"Ledger missing required columns: {sorted(missing)}")

    ledger: dict[str, dict[str, str]] = {}
    for row in rows:
        gsm = row["gsm"].strip()
        if not gsm or gsm in ledger:
            raise ValueError(f"Ledger contains duplicate or empty GSM: {gsm!r}")
        if not parse_bool(row["eligible_for_score_level_validation"], "eligible_for_score_level_validation"):
            raise ValueError(f"{gsm} is not eligible for score-level validation")
        if not row["source_publication_severity_group"].strip():
            raise ValueError(f"{gsm} has no source-publication severity group")
        try:
            cells = int(row["cells_in_supplied_matrix"])
        except ValueError as exc:
            raise ValueError(f"{gsm} has an invalid cells_in_supplied_matrix value") from exc
        if cells < 1:
            raise ValueError(f"{gsm} has no supplied cells")
        ledger[gsm] = row
    return ledger


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def read_cell_mapping(
    archive: tarfile.TarFile, member_name: str, expected_ids: list[str], gsm: str
) -> None:
    member = archive.extractfile(member_name)
    if member is None:
        raise FileNotFoundError(member_name)
    with member, gzip.open(member, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        fields = handle.readline().rstrip("\r\n").split("\t")
        if fields != ["CellName", "CellIndex"]:
            raise ValueError(f"{gsm}: unexpected cell-name header {fields!r}")
        barcodes: list[str] = []
        cell_indices: list[str] = []
        for line_number, raw in enumerate(handle, start=2):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            barcode, separator, cell_index = line.partition("\t")
            if not separator or not barcode or not cell_index or "\t" in cell_index:
                raise ValueError(f"{gsm}: malformed cell-name row {line_number}")
            barcodes.append(barcode)
            cell_indices.append(cell_index)
    if len(barcodes) != len(expected_ids):
        raise ValueError(
            f"{gsm}: cell-name rows {len(barcodes)} do not match matrix columns {len(expected_ids)}"
        )
    if len(barcodes) != len(set(barcodes)) or len(cell_indices) != len(set(cell_indices)):
        raise ValueError(f"{gsm}: duplicate barcode or CellIndex in cell-name mapping")
    if cell_indices != expected_ids:
        raise ValueError(f"{gsm}: CellIndex order does not exactly match the count-matrix header")


def stream_sample(
    archive: tarfile.TarFile,
    count_member_name: str,
    cell_member_name: str,
    metadata: dict[str, str],
    modules: dict[str, list[str]],
    min_mapped_fraction: float,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    """Stream one dense matrix and return score, mapping, gene, and sample rows."""
    gsm = metadata["gsm"].strip()
    source_member = archive.extractfile(count_member_name)
    if source_member is None:
        raise FileNotFoundError(count_member_name)

    module_sets = {module_id: set(genes) for module_id, genes in modules.items()}
    gene_to_modules: dict[str, list[str]] = defaultdict(list)
    for module_id, genes in module_sets.items():
        for gene in genes:
            gene_to_modules[gene].append(module_id)
    mapped: dict[str, set[str]] = defaultdict(set)
    module_sums = {
        module_id: {gene: 0.0 for gene in genes}
        for module_id, genes in modules.items()
    }
    feature_rows: dict[str, dict[str, int]] = {
        module_id: {gene: 0 for gene in genes}
        for module_id, genes in modules.items()
    }
    seen_gene_rows: dict[str, int] = defaultdict(int)
    total_supplied_value = 0.0
    genes_streamed = 0
    values_integer_like = True

    with source_member, gzip.open(source_member, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if len(header) < 2 or header[0] != "gene":
            raise ValueError(f"{gsm}: expected a gene-by-cell TSV header, got {header[:3]!r}")
        cell_ids = header[1:]
        if len(cell_ids) != len(set(cell_ids)):
            raise ValueError(f"{gsm}: duplicate matrix cell IDs")
        expected_cells = int(metadata["cells_in_supplied_matrix"])
        if len(cell_ids) != expected_cells:
            raise ValueError(
                f"{gsm}: matrix columns {len(cell_ids)} do not match ledger cells {expected_cells}"
            )
        read_cell_mapping(archive, cell_member_name, cell_ids, gsm)

        for line_number, raw in enumerate(handle, start=2):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            gene, separator, value_text = line.partition("\t")
            if not separator or not gene.strip():
                raise ValueError(f"{gsm}: malformed count-matrix row {line_number}")
            gene_key = gene.strip().upper()
            values = np.fromstring(value_text, sep="\t", dtype=np.float64)
            if values.size != len(cell_ids):
                raise ValueError(
                    f"{gsm}: row {line_number} has {values.size} values; expected {len(cell_ids)}"
                )
            if not np.all(np.isfinite(values)) or np.any(values < 0):
                raise ValueError(f"{gsm}: non-finite or negative supplied value at row {line_number}")
            if values_integer_like and not np.all(values == np.floor(values)):
                values_integer_like = False
            row_sum = float(values.sum(dtype=np.float64))
            total_supplied_value += row_sum
            seen_gene_rows[gene_key] += 1
            for module_id in gene_to_modules.get(gene_key, []):
                mapped[module_id].add(gene_key)
                feature_rows[module_id][gene_key] += 1
                module_sums[module_id][gene_key] += row_sum
            genes_streamed += 1

    if not math.isfinite(total_supplied_value) or total_supplied_value <= 0:
        raise ValueError(f"{gsm}: nonpositive total supplied matrix value")

    disease_state = metadata["source_publication_severity_group"].strip()
    common = {
        "dataset": metadata["dataset"].strip(),
        "gsm": gsm,
        "donor_id": metadata["donor_id"].strip(),
        "compartment": metadata["compartment"].strip(),
        "disease_state": disease_state,
    }
    mapping_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    score_rows: list[dict[str, object]] = []
    for module_id, genes in modules.items():
        module_mapped = sorted(mapped[module_id])
        mapped_fraction = len(module_mapped) / len(genes)
        duplicate_feature_rows = sum(max(0, feature_rows[module_id][gene] - 1) for gene in genes)
        mapping_rows.append(
            {
                **common,
                "module_id": module_id,
                "configured_genes": len(genes),
                "mapped_genes": len(module_mapped),
                "mapped_fraction": f"{mapped_fraction:.6f}",
                "mapping_pass": mapped_fraction >= min_mapped_fraction,
                "mapped_gene_symbols": ";".join(module_mapped),
                "duplicate_feature_rows": duplicate_feature_rows,
            }
        )
        for gene in genes:
            gene_rows.append(
                {
                    **common,
                    "module_id": module_id,
                    "gene_symbol": gene,
                    "feature_rows_mapped": feature_rows[module_id][gene],
                    "supplied_value_sum": f"{module_sums[module_id][gene]:.8f}",
                    "total_supplied_value": f"{total_supplied_value:.8f}",
                }
            )
        if mapped_fraction >= min_mapped_fraction:
            score = sum(
                math.log1p(1_000_000.0 * module_sums[module_id][gene] / total_supplied_value)
                for gene in module_mapped
            ) / len(module_mapped)
            score_status = "score_available"
        else:
            score = float("nan")
            score_status = "mapping_below_minimum"
        score_rows.append(
            {
                **common,
                "module_id": module_id,
                "module_score_log1p_cpm": "" if math.isnan(score) else f"{score:.8f}",
                "score_status": score_status,
                "mapped_fraction": f"{mapped_fraction:.6f}",
                "included_cells": len(cell_ids),
                "total_supplied_value": f"{total_supplied_value:.8f}",
                "matrix_status": metadata["matrix_status"].strip(),
                "analysis_role": (
                    "processed-matrix donor-level score direction only; all supplied cells nested within "
                    "the presumed sample key; no raw-count or cell-level inference"
                ),
            }
        )
    sample_row = {
        **common,
        "sample": metadata["sample"].strip(),
        "included_cells": len(cell_ids),
        "total_supplied_value": f"{total_supplied_value:.8f}",
        "genes_streamed": genes_streamed,
        "unique_gene_symbols": len(seen_gene_rows),
        "duplicate_gene_symbol_rows": sum(max(0, count - 1) for count in seen_gene_rows.values()),
        "all_supplied_values_integer_like": values_integer_like,
        "matrix_status": metadata["matrix_status"].strip(),
        "barcode_mapping_status": "pass_exact_CellIndex_to_matrix_header",
        "score_scope": "all GEO-supplied cells; no re-QC or raw-count pseudobulk",
    }
    return score_rows, mapping_rows, gene_rows, sample_row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("ledger", type=Path)
    parser.add_argument("module_config", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-mapped-fraction", type=float, default=0.80)
    args = parser.parse_args()
    if not 0 < args.min_mapped_fraction <= 1:
        raise ValueError("--min-mapped-fraction must be in (0, 1]")

    modules, score_direction = read_modules(args.module_config)
    ledger = read_ledger(args.ledger)
    expected_members = {
        f"{row['gsm']}_{row['sample']}.counts.tsv.gz"
        for row in ledger.values()
    } | {
        f"{row['gsm']}_{row['sample']}.cellname.txt.gz"
        for row in ledger.values()
    }
    with tarfile.open(args.archive, "r") as archive:
        actual_members = {member.name for member in archive.getmembers() if member.isfile()}
        if actual_members != expected_members:
            raise ValueError(
                "Archive membership differs from the validated ledger: "
                f"archive_only={sorted(actual_members - expected_members)}, "
                f"ledger_only={sorted(expected_members - actual_members)}"
            )
        score_rows: list[dict[str, object]] = []
        mapping_rows: list[dict[str, object]] = []
        gene_rows: list[dict[str, object]] = []
        sample_rows: list[dict[str, object]] = []
        for gsm, metadata in sorted(ledger.items()):
            count_name = f"{gsm}_{metadata['sample'].strip()}.counts.tsv.gz"
            cell_name = f"{gsm}_{metadata['sample'].strip()}.cellname.txt.gz"
            scores, mappings, genes, sample = stream_sample(
                archive,
                count_name,
                cell_name,
                metadata,
                modules,
                args.min_mapped_fraction,
            )
            score_rows.extend(scores)
            mapping_rows.extend(mappings)
            gene_rows.extend(genes)
            sample_rows.append(sample)

    if len(sample_rows) != len(ledger):
        raise AssertionError("Not every ledger sample produced a score record")
    prefix = args.archive.stem
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / f"{prefix}_module_scores.csv", score_rows, list(score_rows[0]))
    write_csv(args.output_dir / f"{prefix}_module_mapping_audit.csv", mapping_rows, list(mapping_rows[0]))
    write_csv(args.output_dir / f"{prefix}_module_gene_scores.csv", gene_rows, list(gene_rows[0]))
    write_csv(args.output_dir / f"{prefix}_sample_score_ledger.csv", sample_rows, list(sample_rows[0]))
    write_csv(
        args.output_dir / f"{prefix}_module_score_parameters.csv",
        [
            {
                "archive": str(args.archive),
                "ledger": str(args.ledger),
                "module_config": str(args.module_config),
                "min_mapped_fraction": f"{args.min_mapped_fraction:.6f}",
                "samples": len(sample_rows),
                "formula": (
                    "mean over mapped genes of log1p(1e6 * sample-level gene supplied-value sum / "
                    "sample-level total supplied value)"
                ),
                "score_direction": score_direction,
                "inference_unit": "presumed donor-level sample; cells nested",
                "matrix_boundary": (
                    "GEO describes supplied matrices as normalized counts; score-level direction only; "
                    "no raw-count model, pseudobulk count inference, or cell-level test"
                ),
            }
        ],
        [
            "archive",
            "ledger",
            "module_config",
            "min_mapped_fraction",
            "samples",
            "formula",
            "score_direction",
            "inference_unit",
            "matrix_boundary",
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
