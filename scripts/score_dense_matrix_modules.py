"""Score locked modules from a sample-prefixed dense gene-by-cell TSV.

The matrix is streamed one gene row at a time.  Sample/library is the unit of
inference; cells are only columns nested within that key.  This is intended for
GSE153066's GEO-retained combined count matrix, whose barcode prefixes map
one-to-one to the presumed sample/library keys in its ledger.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def read_modules(path: Path) -> tuple[dict[str, list[str]], str]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    modules: dict[str, list[str]] = {}
    for item in payload["modules"]:
        module_id = str(item["module_id"])
        genes = [str(gene).strip().upper() for gene in item["genes"]]
        if not genes or len(genes) != len(set(genes)):
            raise ValueError(f"Invalid or duplicate genes in module {module_id}")
        if module_id in modules:
            raise ValueError(f"Duplicate module id: {module_id}")
        modules[module_id] = genes
    return modules, str(payload.get("score_direction", ""))


def read_ledger(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Empty ledger: {path}")
    required = {"barcode_prefix", "gsm", "dataset", "compartment", "disease_state"}
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"Ledger missing columns: {sorted(missing)}")
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        key = row["barcode_prefix"].strip()
        if not key or key in result:
            raise ValueError(f"Duplicate or empty barcode_prefix: {key!r}")
        result[key] = row
    return result


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path)
    parser.add_argument("ledger", type=Path)
    parser.add_argument("module_config", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-mapped-fraction", type=float, default=0.80)
    args = parser.parse_args()
    if not 0 < args.min_mapped_fraction <= 1:
        raise ValueError("--min-mapped-fraction must be in (0, 1]")

    modules, score_direction = read_modules(args.module_config)
    ledger = read_ledger(args.ledger)
    module_sets = {module: set(genes) for module, genes in modules.items()}
    mapped: dict[str, set[str]] = defaultdict(set)
    counts: dict[str, dict[str, dict[str, int]]] = {
        prefix: {module: {gene: 0 for gene in genes} for module, genes in modules.items()}
        for prefix in ledger
    }
    total_umi = {prefix: 0 for prefix in ledger}
    cells_by_prefix: dict[str, int] = {}
    genes_seen = 0
    duplicate_genes: set[str] = set()

    with gzip.open(args.matrix, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if not header or header[0] != "gene":
            raise ValueError(f"Expected first header field 'gene', got {header[:1]!r}")
        barcodes = header[1:]
        prefixes = [barcode.split("_", 1)[0] for barcode in barcodes]
        unknown = sorted(set(prefixes) - set(ledger))
        missing = sorted(set(ledger) - set(prefixes))
        if unknown or missing:
            raise ValueError(f"Matrix/ledger prefix mismatch: unknown={unknown}, missing={missing}")
        for prefix in ledger:
            cells_by_prefix[prefix] = prefixes.count(prefix)
        prefix_indices: dict[str, list[int]] = defaultdict(list)
        for index, prefix in enumerate(prefixes):
            prefix_indices[prefix].append(index)
        # GEO's combined matrix is block-ordered by sample.  Keep contiguous
        # runs so each row can be reduced with vectorized NumPy slices rather
        # than a Python loop over every cell value.
        runs: dict[str, tuple[int, int]] = {}
        cursor = 0
        for prefix in prefixes:
            if prefix not in runs:
                start = cursor
                end = start
                while end < len(prefixes) and prefixes[end] == prefix:
                    end += 1
                runs[prefix] = (start, end)
                cursor = end
            else:
                # Non-contiguous prefixes are valid but use the indexed path.
                cursor = 0
                break
        contiguous = bool(runs) and cursor == len(prefixes)

        seen_genes: set[str] = set()
        for line_number, raw in enumerate(handle, start=2):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            gene, separator, value_text = line.partition("\t")
            if not separator or not gene:
                raise ValueError(f"Malformed matrix row {line_number}")
            gene_key = gene.strip().upper()
            if gene_key in seen_genes:
                duplicate_genes.add(gene_key)
            seen_genes.add(gene_key)
            # Parse a complete row in optimized C code.  The dense matrix is
            # large, so retaining Python strings or iterating cell-by-cell is
            # prohibitively expensive.
            values = np.fromstring(value_text, sep="\t", dtype=np.float64)
            if values.size != len(barcodes):
                raise ValueError(f"Row {line_number}: {values.size} values, expected {len(barcodes)}")
            if not np.all(np.isfinite(values)) or np.any(values < 0):
                raise ValueError(f"Row {line_number}: non-finite or negative matrix value")
            if np.any(values != np.floor(values)):
                raise ValueError(f"Row {line_number}: non-integer matrix value")
            target_modules = [module for module, genes in module_sets.items() if gene_key in genes]
            if contiguous:
                subtotals = {prefix: float(values[start:end].sum(dtype=np.float64)) for prefix, (start, end) in runs.items()}
            else:
                subtotals = {prefix: float(values[indices].sum(dtype=np.float64)) for prefix, indices in prefix_indices.items()}
            for prefix, subtotal in subtotals.items():
                total_umi[prefix] += subtotal
                for module in target_modules:
                    mapped[module].add(gene_key)
                    counts[prefix][module][gene_key] += subtotal
            genes_seen += 1

    if duplicate_genes:
        raise ValueError(f"Duplicate gene identifiers are not supported: {sorted(duplicate_genes)[:5]}")

    score_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    mapping_rows: list[dict[str, object]] = []
    for prefix, metadata in sorted(ledger.items()):
        for module, genes in modules.items():
            module_mapped = sorted(mapped[module])
            fraction = len(module_mapped) / len(genes)
            mapping_rows.append({
                "dataset": metadata["dataset"], "gsm": metadata["gsm"], "donor_id": metadata.get("donor_id", prefix),
                "library_id": prefix, "compartment": metadata["compartment"], "disease_state": metadata["disease_state"],
                "module_id": module, "configured_genes": len(genes), "mapped_genes": len(module_mapped),
                "mapped_fraction": f"{fraction:.6f}", "mapping_pass": fraction >= args.min_mapped_fraction,
                "mapped_gene_symbols": ";".join(module_mapped), "duplicate_feature_rows": 0,
            })
            for gene in genes:
                gene_rows.append({
                    "dataset": metadata["dataset"], "gsm": metadata["gsm"], "donor_id": metadata.get("donor_id", prefix),
                    "library_id": prefix, "compartment": metadata["compartment"], "disease_state": metadata["disease_state"],
                    "module_id": module, "gene_symbol": gene, "feature_rows_mapped": int(gene in mapped[module]),
                    "pseudobulk_count": counts[prefix][module][gene], "total_umi_included_cells": total_umi[prefix],
                })
            if fraction >= args.min_mapped_fraction and total_umi[prefix] > 0:
                score = sum(math.log1p(1_000_000.0 * counts[prefix][module][gene] / total_umi[prefix]) for gene in module_mapped) / len(module_mapped)
                status = "score_available"
            else:
                score = float("nan")
                status = "mapping_below_minimum_or_zero_library_umi"
            score_rows.append({
                "dataset": metadata["dataset"], "gsm": metadata["gsm"], "donor_id": metadata.get("donor_id", prefix),
                "library_id": prefix, "compartment": metadata["compartment"], "disease_state": metadata["disease_state"],
                "module_id": module, "module_score_log1p_cpm": "" if math.isnan(score) else f"{score:.8f}",
                "score_status": status, "mapped_fraction": f"{fraction:.6f}", "included_cells": cells_by_prefix[prefix],
                "total_umi_included_cells": total_umi[prefix], "matrix_status": metadata.get("matrix_status", ""),
                "analysis_role": "donor-level score support; cells nested within presumed sample/library key",
            })

    library_rows = [{
        "dataset": row["dataset"], "gsm": row["gsm"], "donor_id": row.get("donor_id", prefix), "library_id": prefix,
        "compartment": row["compartment"], "disease_state": row["disease_state"], "included_cells": cells_by_prefix[prefix],
        "total_umi_included_cells": total_umi[prefix], "matrix_status": row.get("matrix_status", ""),
        "prior_processing_note": "GEO-retained sample-prefixed matrix; prior cell filtering disclosed in ledger",
    } for prefix, row in sorted(ledger.items())]
    prefix = args.matrix.stem.replace(".tsv", "")
    write_csv(args.output_dir / f"{prefix}_module_scores.csv", score_rows, list(score_rows[0]))
    write_csv(args.output_dir / f"{prefix}_module_gene_pseudobulk.csv", gene_rows, list(gene_rows[0]))
    write_csv(args.output_dir / f"{prefix}_module_mapping_audit.csv", mapping_rows, list(mapping_rows[0]))
    write_csv(args.output_dir / f"{prefix}_library_pseudobulk_ledger.csv", library_rows, list(library_rows[0]))
    write_csv(args.output_dir / f"{prefix}_module_score_parameters.csv", [{
        "matrix": str(args.matrix), "ledger": str(args.ledger), "module_config": str(args.module_config),
        "min_mapped_fraction": args.min_mapped_fraction, "genes_seen": genes_seen, "libraries": len(ledger),
        "formula": "mean over mapped genes of log1p(1e6 * sample-level gene sum / sample-level total supplied UMI)",
        "score_direction": score_direction, "inference_unit": "presumed donor/library; cells nested",
        "matrix_boundary": "GEO-retained dense matrix; no cell-level inference",
    }], ["matrix", "ledger", "module_config", "min_mapped_fraction", "genes_seen", "libraries", "formula", "score_direction", "inference_unit", "matrix_boundary"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
