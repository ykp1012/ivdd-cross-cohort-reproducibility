"""Apply the locked conservative AF/NP annotation rules to QC-passing cells.

This script streams a 10x TAR archive a second time and only stores per-cell
panel aggregates. It uses GEO anatomical source as the primary label and never
relabels an AF library as NP (or the reverse) from expression alone.
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

from qc_10x_archives import gsm_members, load_ledger, member_map, read_barcodes, read_features


def read_qc_passes(path: Path) -> dict[tuple[str, str], bool]:
    result: dict[tuple[str, str], bool] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[(row["gsm"], row["barcode"])] = row["qc_pass"].lower() == "true"
    return result


def classify(
    source: str,
    scores: dict[str, float],
    detected: dict[str, int],
    config: dict[str, object],
) -> tuple[str, str, bool, str]:
    rules = config["rules"]
    assert isinstance(rules, dict)
    min_pos = int(rules["minimum_positive_detected"])
    min_resident = float(rules["minimum_resident_score"])
    min_nonresident = float(rules["minimum_nonresident_score"])
    min_nonresident_detected = int(rules["minimum_nonresident_detected"])
    margin = float(rules["resident_margin"])
    panel_pass = {
        panel: detected[panel] >= min_nonresident_detected and scores[panel] >= min_nonresident
        for panel in scores
    }
    nonresident = [panel for panel in ("immune_exclusion", "endothelial_exclusion", "mural_exclusion", "erythroid_exclusion") if panel_pass[panel]]
    source_panel = "AF_support" if source.upper() == "AF" else "NP_support" if source.upper() == "NP" else ""
    alternative_panel = "NP_support" if source_panel == "AF_support" else "AF_support"
    source_support = bool(
        source_panel
        and detected[source_panel] >= min_pos
        and scores[source_panel] >= min_resident
    )
    source_disagreement = bool(source_panel and scores[alternative_panel] > scores[source_panel] + margin)
    if source_support and nonresident:
        return "mixed_or_nonresident", "source and nonresident panels pass", source_disagreement, "supporting"
    if nonresident:
        label = nonresident[0].replace("_exclusion", "")
        return f"nonresident_{label}", "fixed exclusion panel passes", source_disagreement, "not_applicable"
    if source_panel:
        if source_support:
            return f"source_{source.upper()}_nonexcluded", "source-labelled and source-supporting", source_disagreement, "supporting"
        return (
            f"source_{source.upper()}_nonexcluded",
            "source-labelled; support_insufficient",
            source_disagreement,
            "support_insufficient",
        )
    return "ambiguous", "missing anatomical source label", source_disagreement, "not_applicable"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("ledger", type=Path)
    parser.add_argument("qc_cells", type=Path)
    parser.add_argument("panel_config", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.panel_config.read_text(encoding="utf-8"))
    panels = config["panels"]
    assert isinstance(panels, dict)
    panel_names = list(panels)
    gene_to_panels: defaultdict[str, list[str]] = defaultdict(list)
    for panel, genes in panels.items():
        assert isinstance(genes, list)
        for gene in genes:
            gene_to_panels[str(gene).upper()].append(panel)
    ledger = load_ledger(args.ledger)
    qc_passes = read_qc_passes(args.qc_cells)
    per_cell_rows: list[dict[str, object]] = []
    donor_counts: defaultdict[tuple[str, str, str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    mapping_rows: list[dict[str, object]] = []

    with tarfile.open(args.archive, "r") as tar:
        members = member_map(tar)
        gsms = sorted(gsm for gsm in ledger if any(name.startswith(f"{gsm}_") for name in members))
        for gsm in gsms:
            metadata = ledger[gsm]
            parts = gsm_members(members, gsm)
            barcodes = read_barcodes(parts["barcodes"], tar)
            genes, _ = read_features(parts["features"], tar)
            gene_panels = [gene_to_panels.get(gene.upper(), []) for gene in genes]
            present = {panel: 0 for panel in panel_names}
            panel_rows: dict[str, dict[str, int]] = {panel: {} for panel in panel_names}
            for feature, gene in enumerate(genes):
                for panel in gene_to_panels.get(gene.upper(), []):
                    if gene.upper() not in panel_rows[panel]:
                        panel_rows[panel][gene.upper()] = len(panel_rows[panel])
            present = {panel: len(panel_rows[panel]) for panel in panel_names}
            mapping_rows.extend(
                {
                    "dataset": metadata.get("dataset", ""), "gsm": gsm, "panel": panel,
                    "configured_gene_count": len(panels[panel]), "measured_gene_count": present[panel],
                    "measured_fraction": f"{present[panel] / len(panels[panel]):.6f}",
                }
                for panel in panel_names
            )
            total = np.zeros(len(barcodes), dtype=np.int64)
            # The target gene matrices are small (fixed panels by one library),
            # avoiding a dense whole-transcriptome allocation.
            panel_counts = {
                panel: np.zeros((present[panel], len(barcodes)), dtype=np.uint32)
                for panel in panel_names
            }
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
                rows, cols, _ = (int(value) for value in line.split())
                if rows != len(genes) or cols != len(barcodes):
                    raise ValueError(f"Dimension mismatch in {parts['matrix'].name}")
                for raw in handle:
                    if not raw.strip() or raw.startswith("%"):
                        continue
                    feature, cell, value = (int(value) for value in raw.split()[:3])
                    feature -= 1
                    cell -= 1
                    if value <= 0:
                        continue
                    total[cell] += value
                    gene_symbol = genes[feature].upper()
                    for panel in gene_panels[feature]:
                        panel_counts[panel][panel_rows[panel][gene_symbol], cell] += value
            safe_total = np.maximum(total, 1)
            panel_scores = {
                panel: np.log1p(panel_counts[panel] * (1_000_000.0 / safe_total)[None, :]).mean(axis=0)
                if present[panel] else np.zeros(len(barcodes), dtype=float)
                for panel in panel_names
            }
            panel_detected = {
                panel: (panel_counts[panel] > 0).sum(axis=0)
                if present[panel] else np.zeros(len(barcodes), dtype=int)
                for panel in panel_names
            }
            donor_id = metadata.get("donor_id", "") or metadata.get("patient_id", "") or gsm
            source = metadata.get("compartment", "")
            key = (metadata.get("dataset", ""), donor_id, source, metadata.get("disease_state", ""))
            for cell, barcode in enumerate(barcodes):
                if not qc_passes.get((gsm, barcode), False):
                    continue
                scores = {panel: float(panel_scores[panel][cell]) for panel in panel_names}
                detected = {panel: int(panel_detected[panel][cell]) for panel in panel_names}
                label, reason, source_disagreement, support_status = classify(source, scores, detected, config)
                include = label == f"source_{source.upper()}_nonexcluded"
                donor_counts[key]["qc_cells"] += 1
                donor_counts[key][label] += 1
                record: dict[str, object] = {
                    "dataset": metadata.get("dataset", ""), "gsm": gsm, "donor_id": donor_id,
                    "compartment_source": source, "disease_state": metadata.get("disease_state", ""),
                    "barcode": barcode, "annotation_label": label, "annotation_reason": reason,
                    "source_marker_discordant": source_disagreement, "source_support_status": support_status,
                    "compartment_pseudobulk_include": include,
                }
                for panel in panel_names:
                    record[f"{panel}_score"] = f"{scores[panel]:.6f}"
                    record[f"{panel}_detected_genes"] = detected[panel]
                per_cell_rows.append(record)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cell_path = args.output_dir / f"{args.archive.stem}_cell_annotation.csv.gz"
    fields = list(per_cell_rows[0]) if per_cell_rows else ["gsm"]
    with gzip.open(cell_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(per_cell_rows)
    donor_rows: list[dict[str, object]] = []
    for key, counts in sorted(donor_counts.items()):
        dataset, donor_id, source, disease = key
        source_nonexcluded = counts.get(f"source_{source.upper()}_nonexcluded", 0)
        donor_rows.append({
            "dataset": dataset, "donor_id": donor_id, "compartment_source": source,
            "disease_state": disease, "qc_cells": counts.get("qc_cells", 0),
            "source_nonexcluded_cells": source_nonexcluded,
            "source_nonexcluded_threshold_20_pass": source_nonexcluded >= 20,
            "source_nonexcluded_threshold_30_pass": source_nonexcluded >= 30,
            "source_nonexcluded_threshold_50_pass": source_nonexcluded >= 50,
            "ambiguous_cells": counts.get("ambiguous", 0),
            "mixed_or_nonresident_cells": counts.get("mixed_or_nonresident", 0),
            "nonresident_immune_cells": counts.get("nonresident_immune", 0),
            "nonresident_endothelial_cells": counts.get("nonresident_endothelial", 0),
            "nonresident_mural_cells": counts.get("nonresident_mural", 0),
            "nonresident_erythroid_cells": counts.get("nonresident_erythroid", 0),
        })
    with (args.output_dir / f"{args.archive.stem}_donor_annotation_sensitivity.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(donor_rows[0]) if donor_rows else ["donor_id"])
        writer.writeheader()
        writer.writerows(donor_rows)
    with (args.output_dir / f"{args.archive.stem}_marker_mapping.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(mapping_rows[0]) if mapping_rows else ["gsm"])
        writer.writeheader()
        writer.writerows(mapping_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
