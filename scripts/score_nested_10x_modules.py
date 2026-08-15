"""QC, source-restricted selection, and module scoring for nested 10x-like TARs."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import shutil
import tarfile
import tempfile
from collections import defaultdict
from pathlib import Path


def parse_config(module_path: Path, panel_path: Path) -> tuple[dict[str, list[str]], dict[str, list[str]], dict[str, float | int]]:
    module_cfg = json.loads(module_path.read_text(encoding="utf-8"))
    panel_cfg = json.loads(panel_path.read_text(encoding="utf-8"))
    modules = {m["module_id"]: [g.upper() for g in m["genes"]] for m in module_cfg["modules"]}
    panels = {k: [g.upper() for g in v] for k, v in panel_cfg["panels"].items()}
    return modules, panels, panel_cfg["rules"]


def extract_nested_to_disk(
    outer: tarfile.TarFile, member: tarfile.TarInfo, workdir: Path
) -> tuple[list[str], list[str], Path]:
    """Extract one nested matrix to disk, keeping peak RAM bounded."""
    nested_path = workdir / "nested.tar.gz"
    stream = outer.extractfile(member)
    if stream is None:
        raise OSError(member.name)
    with stream, nested_path.open("wb") as handle:
        shutil.copyfileobj(stream, handle, length=1024 * 1024)
    matrix_path = workdir / "matrix.mtx"
    with tarfile.open(nested_path, mode="r:gz") as inner:
        names = {Path(m.name).name: m for m in inner.getmembers() if m.isfile()}
        if not {"genes.tsv", "barcodes.tsv", "matrix.mtx"}.issubset(names):
            raise OSError(f"Incomplete nested matrix: {member.name}")
        genes_stream = inner.extractfile(names["genes.tsv"])
        barcodes_stream = inner.extractfile(names["barcodes.tsv"])
        matrix_stream = inner.extractfile(names["matrix.mtx"])
        if genes_stream is None or barcodes_stream is None or matrix_stream is None:
            raise OSError(f"Incomplete nested matrix: {member.name}")
        with genes_stream:
            genes = [line.decode("utf-8", errors="replace").rstrip("\r\n").split("\t")[-1].upper() for line in genes_stream if line.strip()]
        with barcodes_stream:
            barcodes = [line.decode("utf-8", errors="replace").rstrip("\r\n").split("\t", 1)[0] for line in barcodes_stream if line.strip()]
        with matrix_stream, matrix_path.open("wb") as handle:
            shutil.copyfileobj(matrix_stream, handle, length=1024 * 1024)
    return genes, barcodes, matrix_path


def iter_entries(matrix_path: Path, n_genes: int, n_cells: int):
    """Yield Matrix Market entries from a disk-backed binary stream."""
    handle = matrix_path.open("rb")
    first = handle.readline().strip()
    try:
        if not first.startswith(b"%%MatrixMarket matrix coordinate"):
            raise ValueError("Unexpected Matrix Market header")
        line = handle.readline().strip()
        while line.startswith(b"%"):
            line = handle.readline().strip()
        if not line:
            raise ValueError("Missing Matrix Market dimensions")
        rows, cols, _ = (int(value) for value in line.split())
        if (rows, cols) != (n_genes, n_cells):
            raise ValueError(f"Matrix dimensions {rows}x{cols} do not match {n_genes}x{n_cells}")
        for raw in handle:
            if raw.strip() and not raw.startswith(b"%"):
                fields = raw.split()
                if len(fields) < 3:
                    raise ValueError(f"Malformed Matrix Market entry in {matrix_path}")
                feature, cell, value = (int(value) for value in fields[:3])
                yield feature - 1, cell - 1, value
    finally:
        handle.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("ledger", type=Path)
    parser.add_argument("module_config", type=Path)
    parser.add_argument("--panel-config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-umi", type=int, default=500)
    parser.add_argument("--max-mt-pct", type=float, default=20.0)
    parser.add_argument("--min-mapped-fraction", type=float, default=0.80)
    args = parser.parse_args()
    modules, panels, rules = parse_config(args.module_config, args.panel_config)
    with args.ledger.open(newline="", encoding="utf-8") as handle:
        ledger = {row["gsm"]: row for row in csv.DictReader(handle)}
    panel_genes = {panel: set(genes) for panel, genes in panels.items()}
    module_genes = {module: set(genes) for module, genes in modules.items()}
    annotation_rows: list[dict[str, object]] = []
    score_rows: list[dict[str, object]] = []
    library_rows: list[dict[str, object]] = []
    donor_rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="ivdd_nested_") as temp_root, tarfile.open(args.archive, "r") as outer:
        members = sorted((m for m in outer.getmembers() if m.isfile()), key=lambda m: m.name)
        for member_index, member in enumerate(members):
            gsm = member.name.split("_", 1)[0]
            if gsm not in ledger:
                raise ValueError(f"{gsm} absent from ledger")
            metadata = ledger[gsm]
            workdir = Path(temp_root) / str(member_index)
            workdir.mkdir(parents=True, exist_ok=True)
            genes, barcodes, matrix_path = extract_nested_to_disk(outer, member, workdir)
            mt = [gene.startswith("MT-") for gene in genes]
            total = [0] * len(barcodes); detected = [0] * len(barcodes); mt_total = [0] * len(barcodes)
            for feature, cell, value in iter_entries(matrix_path, len(genes), len(barcodes)):
                if value <= 0: continue
                total[cell] += value; detected[cell] += 1
                if mt[feature]: mt_total[cell] += value
            source = metadata.get("compartment", "NP")
            donor = metadata.get("presumed_donor_or_library_key", gsm)
            selected: set[int] = set(); counts = defaultdict(int)
            panel_counts = {panel: [[0] * len(barcodes) for _ in panel_genes[panel]] for panel in panel_genes}
            panel_index = {panel: {gene: i for i, gene in enumerate(sorted(panel_genes[panel]))} for panel in panel_genes}
            gene_to_panels = defaultdict(list)
            for panel, geneset in panel_genes.items():
                for gene in geneset: gene_to_panels[gene].append(panel)
            for feature, cell, value in iter_entries(matrix_path, len(genes), len(barcodes)):
                if value <= 0: continue
                for panel in gene_to_panels.get(genes[feature], []):
                    panel_counts[panel][panel_index[panel][genes[feature]]][cell] += value
            for cell, barcode in enumerate(barcodes):
                pct = 0.0 if total[cell] == 0 else 100 * mt_total[cell] / total[cell]
                reasons = []
                if detected[cell] < args.min_genes: reasons.append("detected_genes_below_min")
                if total[cell] < args.min_umi: reasons.append("umi_below_min")
                if pct >= args.max_mt_pct: reasons.append("mt_pct_at_or_above_max")
                if reasons:
                    counts["qc_excluded"] += 1
                    continue
                scores = {panel: sum(math.log1p(1_000_000 * row[cell] / max(1, total[cell])) for row in panel_counts[panel]) / max(1, len(panel_counts[panel])) for panel in panel_counts}
                det = {panel: sum(row[cell] > 0 for row in panel_counts[panel]) for panel in panel_counts}
                nonresident = [panel for panel in ("immune_exclusion", "endothelial_exclusion", "mural_exclusion", "erythroid_exclusion") if det[panel] >= int(rules["minimum_nonresident_detected"]) and scores[panel] >= float(rules["minimum_nonresident_score"])]
                if nonresident:
                    label = "mixed_or_nonresident" if scores.get("NP_support", 0) >= float(rules["minimum_resident_score"]) and det.get("NP_support", 0) >= int(rules["minimum_positive_detected"]) else "nonresident"
                    counts[label] += 1
                    include = False
                else:
                    label = "source_NP_nonexcluded"; counts[label] += 1; selected.add(cell); include = True
                annotation_rows.append({"dataset": metadata.get("dataset", ""), "gsm": gsm, "donor_id": donor, "compartment_source": source, "severity_group": metadata.get("severity_group", ""), "barcode": barcode, "total_umi": total[cell], "detected_genes": detected[cell], "pct_mt": f"{pct:.6f}", "annotation_label": label, "compartment_pseudobulk_include": include, "immune_score": f"{scores['immune_exclusion']:.6f}", "mural_score": f"{scores['mural_exclusion']:.6f}"})
            total_selected_umi = 0
            module_counts = {module: {gene: 0 for gene in geneset} for module, geneset in module_genes.items()}
            feature_modules = defaultdict(list)
            for feature, gene in enumerate(genes):
                for module, geneset in module_genes.items():
                    if gene in geneset: feature_modules[feature].append((module, gene))
            for feature, cell, value in iter_entries(matrix_path, len(genes), len(barcodes)):
                if cell not in selected or value <= 0: continue
                total_selected_umi += value
                for module, gene in feature_modules.get(feature, []): module_counts[module][gene] += value
            library_rows.append({"dataset": metadata.get("dataset", ""), "gsm": gsm, "donor_id": donor, "compartment": source, "severity_group": metadata.get("severity_group", ""), "included_cells": len(selected), "total_umi_included_cells": total_selected_umi})
            for module, geneset in module_genes.items():
                mapped = [gene for gene in modules[module] if gene in set(genes)]
                fraction = len(mapped) / len(geneset)
                score = sum(math.log1p(1_000_000 * module_counts[module][gene] / max(1, total_selected_umi)) for gene in mapped) / len(mapped) if mapped and fraction >= args.min_mapped_fraction and total_selected_umi else float("nan")
                score_rows.append({"dataset": metadata.get("dataset", ""), "gsm": gsm, "donor_id": donor, "compartment": source, "severity_group": metadata.get("severity_group", ""), "module_id": module, "module_score_log1p_cpm": "" if math.isnan(score) else f"{score:.8f}", "mapped_fraction": f"{fraction:.6f}", "included_cells": len(selected), "total_umi_included_cells": total_selected_umi, "score_status": "score_available" if not math.isnan(score) else "mapping_below_minimum_or_zero_umi"})
            donor_rows.append({"dataset": metadata.get("dataset", ""), "donor_id": donor, "compartment": source, "severity_group": metadata.get("severity_group", ""), "total_cells": len(barcodes), "qc_cells": len(barcodes) - counts["qc_excluded"], "source_nonexcluded_cells": len(selected), "source_nonexcluded_threshold_20_pass": len(selected) >= 20, "source_nonexcluded_threshold_30_pass": len(selected) >= 30, "source_nonexcluded_threshold_50_pass": len(selected) >= 50, "mixed_or_nonresident_cells": counts["mixed_or_nonresident"] + counts["nonresident"]})
    args.output_dir.mkdir(parents=True, exist_ok=True); prefix = args.archive.stem
    def write(path, rows):
        with path.open("w", newline="", encoding="utf-8") as h:
            w = csv.DictWriter(h, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with gzip.open(args.output_dir / f"{prefix}_cell_annotation.csv.gz", "wt", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(annotation_rows[0])); w.writeheader(); w.writerows(annotation_rows)
    write(args.output_dir / f"{prefix}_module_scores.csv", score_rows); write(args.output_dir / f"{prefix}_library_pseudobulk_ledger.csv", library_rows); write(args.output_dir / f"{prefix}_donor_annotation_sensitivity.csv", donor_rows)
    return 0


if __name__ == "__main__": raise SystemExit(main())
