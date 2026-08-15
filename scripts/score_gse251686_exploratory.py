"""Audit-gated, exploratory module scoring for GSE251686.

The scorer consumes only stream-integrity-passing libraries from the audited
ledger.  It never opens the malformed GSM7986002 payload, validates the
remaining nested Matrix Market streams against the machine-readable audits,
and records Ensembl-to-symbol mapping and duplicate-symbol handling.  The
output is deliberately a separate exploratory score set and is not a
confirmatory cohort or an input to the default cross-cohort summary.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import shutil
import tarfile
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


BAD_GSM = "GSM7986002"
REQUIRED_LEDGER = {
    "gsm",
    "outer_matrix_member",
    "compartment",
    "severity_group",
    "analysis_inclusion",
    "stream_audit_text_integrity_pass",
    "stream_audit_matrix_payload_legal",
}
REQUIRED_STREAM = {
    "gsm",
    "matrix_rows",
    "matrix_columns",
    "matrix_nnz_header",
    "coordinate_lines_observed",
    "valid_coordinate_lines",
    "malformed_coordinate_lines",
    "out_of_range_coordinates",
    "negative_values",
    "zero_values",
    "nul_bytes",
    "dimension_match",
    "line_count_matches_header",
    "valid_count_matches_header",
    "text_integrity_pass",
}
REQUIRED_NESTED = {
    "gsm",
    "features",
    "barcodes",
    "matrix_rows",
    "matrix_columns",
    "matrix_nnz",
    "dimension_check_pass",
}
REQUIRED_IDENTIFIER = {
    "gsm",
    "unique_ensembl_feature_ids",
    "duplicate_ensembl_feature_ids",
    "unique_gene_symbols",
    "duplicate_gene_symbols",
    "unique_barcodes",
    "duplicate_barcodes",
    "recommended_cross_library_feature_key",
    "identifier_check_pass",
}


def parse_bool(value: str, field: str, key: str) -> bool:
    normalized = str(value).strip().lower()
    if normalized not in {"true", "false"}:
        raise ValueError(f"{field} for {key} must be true or false, got {value!r}")
    return normalized == "true"


def read_unique_csv(path: Path, key: str, required: set[str]) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        missing = required - fields
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        rows: dict[str, dict[str, str]] = {}
        for row in reader:
            value = str(row.get(key, "")).strip()
            if not value or value in rows:
                raise ValueError(f"{path}: duplicate or empty {key}: {value!r}")
            rows[value] = {name: str(item or "").strip() for name, item in row.items()}
    if not rows:
        raise ValueError(f"{path}: no rows")
    return rows


def read_modules(path: Path) -> dict[str, list[str]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    modules: dict[str, list[str]] = {}
    for item in payload.get("modules", []):
        module_id = str(item["module_id"]).strip()
        genes = [str(gene).strip().upper() for gene in item["genes"]]
        if not module_id or not genes or len(genes) != len(set(genes)):
            raise ValueError(f"Invalid module definition: {item!r}")
        if module_id in modules:
            raise ValueError(f"Duplicate module ID: {module_id}")
        modules[module_id] = genes
    if not modules:
        raise ValueError(f"No modules found in {path}")
    return modules


def read_panels(path: Path) -> tuple[dict[str, list[str]], dict[str, object]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    panels = {
        str(name): [str(gene).strip().upper() for gene in genes]
        for name, genes in payload["panels"].items()
    }
    rules = payload["rules"]
    return panels, rules


def validate_audit_alignment(
    ledger: dict[str, dict[str, str]],
    stream: dict[str, dict[str, str]],
    nested: dict[str, dict[str, str]],
    identifiers: dict[str, dict[str, str]],
) -> list[str]:
    all_sets = {name: set(rows) for name, rows in {
        "ledger": ledger,
        "stream": stream,
        "nested": nested,
        "identifier": identifiers,
    }.items()}
    expected = all_sets["ledger"]
    for name, values in all_sets.items():
        if values != expected:
            raise ValueError(
                f"GSE251686 audit GSM mismatch for {name}: "
                f"missing={sorted(expected - values)}, extra={sorted(values - expected)}"
            )
    selected: list[str] = []
    for gsm, row in sorted(ledger.items()):
        text_pass = parse_bool(row["stream_audit_text_integrity_pass"], "ledger integrity", gsm)
        payload_pass = parse_bool(row["stream_audit_matrix_payload_legal"], "ledger payload", gsm)
        stream_row = stream[gsm]
        nested_row = nested[gsm]
        identifier_row = identifiers[gsm]
        stream_pass = parse_bool(stream_row["text_integrity_pass"], "stream integrity", gsm)
        if gsm == BAD_GSM:
            if text_pass or payload_pass or stream_pass:
                raise ValueError(f"{BAD_GSM} is permanently excluded but an audit marks it as passing")
            continue
        if not (text_pass and payload_pass and stream_pass):
            raise ValueError(f"Non-excluded {gsm} lacks both required integrity passes")
        if row["analysis_inclusion"] != "included after stream-integrity audit":
            raise ValueError(f"{gsm} has passing integrity flags but is not marked included")
        if nested_row["features"] != stream_row["features_from_genes_tsv"]:
            raise ValueError(f"{gsm}: nested feature count disagrees with stream audit")
        if nested_row["barcodes"] != stream_row["barcodes_from_barcodes_tsv"]:
            raise ValueError(f"{gsm}: nested barcode count disagrees with stream audit")
        if nested_row["matrix_rows"] != stream_row["matrix_rows"] or nested_row["matrix_columns"] != stream_row["matrix_columns"]:
            raise ValueError(f"{gsm}: nested dimensions disagree with stream audit")
        if nested_row["matrix_nnz"] != stream_row["matrix_nnz_header"]:
            raise ValueError(f"{gsm}: nested nnz disagrees with stream audit")
        if not parse_bool(nested_row["dimension_check_pass"], "nested dimension", gsm):
            raise ValueError(f"{gsm}: nested dimension audit failed")
        if not parse_bool(identifier_row["identifier_check_pass"], "identifier audit", gsm):
            raise ValueError(f"{gsm}: identifier audit failed")
        if identifier_row["recommended_cross_library_feature_key"].lower() != "ensembl feature id":
            raise ValueError(f"{gsm}: cross-library feature key is not Ensembl feature ID")
        if int(identifier_row["duplicate_ensembl_feature_ids"] or 0) != 0:
            raise ValueError(f"{gsm}: duplicate Ensembl feature IDs are not allowed")
        if int(identifier_row["duplicate_barcodes"] or 0) != 0:
            raise ValueError(f"{gsm}: duplicate barcodes are not allowed")
        selected.append(gsm)
    if not selected:
        raise ValueError("No GSE251686 library passed the hard integrity gate")
    return selected


def extract_nested(
    outer: tarfile.TarFile, member: tarfile.TarInfo, workdir: Path
) -> tuple[Path, Path, Path]:
    nested_path = workdir / "nested.tar.gz"
    source = outer.extractfile(member)
    if source is None:
        raise OSError(member.name)
    with source, nested_path.open("wb") as target:
        shutil.copyfileobj(source, target, length=1024 * 1024)
    paths: dict[str, Path] = {}
    with tarfile.open(nested_path, mode="r:gz") as inner:
        members: dict[str, tarfile.TarInfo] = {}
        for item in inner.getmembers():
            if item.isfile():
                base = Path(item.name).name
                if base in members:
                    raise ValueError(f"{member.name}: duplicate nested member {base}")
                members[base] = item
        required = {"genes.tsv", "barcodes.tsv", "matrix.mtx"}
        if not required <= set(members):
            raise ValueError(f"{member.name}: missing nested matrix members")
        for base in required:
            source_member = inner.extractfile(members[base])
            if source_member is None:
                raise OSError(f"{member.name}: {base}")
            target = workdir / base
            with source_member, target.open("wb") as handle:
                shutil.copyfileobj(source_member, handle, length=1024 * 1024)
            paths[base] = target
    return paths["genes.tsv"], paths["barcodes.tsv"], paths["matrix.mtx"]


def read_features(path: Path) -> tuple[list[str], list[str], dict[str, int]]:
    feature_ids: list[str] = []
    symbols: list[str] = []
    with path.open("rb") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip(b"\r\n").split(b"\t")
            if len(fields) < 2:
                raise ValueError(f"{path}: malformed genes.tsv row {line_number}")
            feature_id = fields[0].decode("utf-8", "strict").strip()
            symbol = fields[-1].decode("utf-8", "strict").strip().upper()
            if not feature_id:
                raise ValueError(f"{path}: empty Ensembl feature ID at row {line_number}")
            feature_ids.append(feature_id)
            symbols.append(symbol)
    if len(feature_ids) != len(set(feature_ids)):
        raise ValueError(f"{path}: duplicate Ensembl feature IDs")
    symbol_counts = Counter(symbol for symbol in symbols if symbol)
    return feature_ids, symbols, symbol_counts


def read_barcodes(path: Path) -> list[str]:
    with path.open("rb") as handle:
        barcodes = [line.rstrip(b"\r\n").split(b"\t", 1)[0].decode("utf-8", "strict") for line in handle if line.strip()]
    if len(barcodes) != len(set(barcodes)):
        raise ValueError(f"{path}: duplicate barcodes")
    return barcodes


def iter_matrix(path: Path, n_features: int, n_cells: int, expected_nnz: int):
    with path.open("rb") as handle:
        first = handle.readline()
        if not first.startswith(b"%%MatrixMarket matrix coordinate"):
            raise ValueError(f"{path}: invalid Matrix Market header")
        dimensions = handle.readline()
        while dimensions.startswith(b"%"):
            dimensions = handle.readline()
        try:
            rows, columns, declared_nnz = (int(value) for value in dimensions.split())
        except ValueError as exc:
            raise ValueError(f"{path}: invalid Matrix Market dimensions") from exc
        if (rows, columns, declared_nnz) != (n_features, n_cells, expected_nnz):
            raise ValueError(f"{path}: dimensions/nnz disagree with audit")
        entries = 0
        for raw in handle:
            if not raw.strip() or raw.startswith(b"%"):
                continue
            if b"\x00" in raw:
                raise ValueError(f"{path}: NUL byte encountered in a library that passed audit")
            fields = raw.split()
            if len(fields) < 3:
                raise ValueError(f"{path}: malformed coordinate line")
            try:
                feature, cell, value = (int(item) for item in fields[:3])
            except ValueError as exc:
                raise ValueError(f"{path}: non-integer coordinate line") from exc
            if not (1 <= feature <= n_features and 1 <= cell <= n_cells) or value < 0:
                raise ValueError(f"{path}: illegal coordinate/value")
            entries += 1
            yield feature - 1, cell - 1, value
        if entries != expected_nnz:
            raise ValueError(f"{path}: read {entries} coordinates, expected {expected_nnz}")


def classify_cell(
    source: str,
    scores: dict[str, float],
    detected: dict[str, int],
    rules: dict[str, object],
) -> str:
    min_pos = int(rules["minimum_positive_detected"])
    min_resident = float(rules["minimum_resident_score"])
    min_nonresident = float(rules["minimum_nonresident_score"])
    min_nonresident_detected = int(rules["minimum_nonresident_detected"])
    source_panel = "NP_support" if source.upper() == "NP" else "AF_support"
    nonresident_panels = (
        "immune_exclusion",
        "endothelial_exclusion",
        "mural_exclusion",
        "erythroid_exclusion",
    )
    nonresident = [
        panel for panel in nonresident_panels
        if detected[panel] >= min_nonresident_detected and scores[panel] >= min_nonresident
    ]
    resident = detected[source_panel] >= min_pos and scores[source_panel] >= min_resident
    if nonresident and resident:
        return "mixed_or_nonresident"
    if nonresident:
        return f"nonresident_{nonresident[0].replace('_exclusion', '')}"
    return f"source_{source.upper()}_nonexcluded"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"Cannot write empty output: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("ledger", type=Path)
    parser.add_argument("module_config", type=Path)
    parser.add_argument("--stream-audit", type=Path, required=True)
    parser.add_argument("--nested-audit", type=Path, required=True)
    parser.add_argument("--identifier-audit", type=Path, required=True)
    parser.add_argument("--panel-config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-umi", type=int, default=500)
    parser.add_argument("--max-mt-pct", type=float, default=20.0)
    parser.add_argument("--min-retained-cells", type=int, choices=(20, 30, 50), default=30)
    parser.add_argument("--min-mapped-fraction", type=float, default=0.80)
    args = parser.parse_args()
    if not 0 < args.min_mapped_fraction <= 1:
        raise ValueError("--min-mapped-fraction must be in (0, 1]")

    ledger = read_unique_csv(args.ledger, "gsm", REQUIRED_LEDGER)
    stream = read_unique_csv(args.stream_audit, "gsm", REQUIRED_STREAM)
    nested = read_unique_csv(args.nested_audit, "gsm", REQUIRED_NESTED)
    identifiers = read_unique_csv(args.identifier_audit, "gsm", REQUIRED_IDENTIFIER)
    selected_gsms = validate_audit_alignment(ledger, stream, nested, identifiers)
    modules = read_modules(args.module_config)
    panels, rules = read_panels(args.panel_config)
    panel_names = list(panels)
    module_names = list(modules)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    score_rows: list[dict[str, object]] = []
    mapping_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    library_rows: list[dict[str, object]] = []
    annotation_rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="ivdd_gse251686_score_") as temp_root:
        with tarfile.open(args.archive, mode="r") as outer:
            members = {member.name: member for member in outer.getmembers() if member.isfile()}
            for gsm in selected_gsms:
                row = ledger[gsm]
                member_name = row["outer_matrix_member"]
                if member_name not in members:
                    raise ValueError(f"{gsm}: ledger outer member is absent from archive")
                # BAD_GSM is intentionally never passed to extractfile().
                workdir = Path(temp_root) / gsm
                workdir.mkdir(parents=True, exist_ok=True)
                genes_path, barcodes_path, matrix_path = extract_nested(outer, members[member_name], workdir)
                feature_ids, symbols, symbol_counts = read_features(genes_path)
                barcodes = read_barcodes(barcodes_path)
                gsm_stream = stream[gsm]
                n_features = len(feature_ids)
                n_cells = len(barcodes)
                expected_nnz = int(gsm_stream["matrix_nnz_header"])
                if (n_features, n_cells) != (int(gsm_stream["matrix_rows"]), int(gsm_stream["matrix_columns"])):
                    raise ValueError(f"{gsm}: extracted dimensions disagree with stream audit")

                panel_symbols = {
                    panel: sorted(set(symbol for symbol in symbols if symbol in set(genes)))
                    for panel, genes in panels.items()
                }
                panel_counts = {
                    panel: {symbol: [0] * n_cells for symbol in measured}
                    for panel, measured in panel_symbols.items()
                }
                feature_panels: dict[int, list[tuple[str, str]]] = defaultdict(list)
                for feature_index, symbol in enumerate(symbols):
                    if not symbol:
                        continue
                    for panel, measured in panel_symbols.items():
                        if symbol in measured:
                            feature_panels[feature_index].append((panel, symbol))
                totals = [0] * n_cells
                detected_genes = [0] * n_cells
                mt_totals = [0] * n_cells
                for feature, cell, value in iter_matrix(matrix_path, n_features, n_cells, expected_nnz):
                    if value <= 0:
                        continue
                    totals[cell] += value
                    detected_genes[cell] += 1
                    if symbols[feature].startswith("MT-"):
                        mt_totals[cell] += value
                    for panel, symbol in feature_panels.get(feature, []):
                        panel_counts[panel][symbol][cell] += value

                include: list[bool] = [False] * n_cells
                labels: list[str] = ["qc_excluded"] * n_cells
                for cell, barcode in enumerate(barcodes):
                    mt_pct = 100.0 * mt_totals[cell] / totals[cell] if totals[cell] else 0.0
                    qc_reasons: list[str] = []
                    if detected_genes[cell] < args.min_genes:
                        qc_reasons.append("detected_genes_below_min")
                    if totals[cell] < args.min_umi:
                        qc_reasons.append("umi_below_min")
                    if mt_pct >= args.max_mt_pct:
                        qc_reasons.append("mt_pct_at_or_above_max")
                    if not qc_reasons:
                        panel_scores = {
                            panel: sum(
                                math.log1p(1_000_000.0 * panel_counts[panel][symbol][cell] / max(1, totals[cell]))
                                for symbol in panel_symbols[panel]
                            ) / len(panel_symbols[panel]) if panel_symbols[panel] else 0.0
                            for panel in panel_names
                        }
                        panel_detected = {
                            panel: sum(panel_counts[panel][symbol][cell] > 0 for symbol in panel_symbols[panel])
                            for panel in panel_names
                        }
                        label = classify_cell("NP", panel_scores, panel_detected, rules)
                        labels[cell] = label
                        include[cell] = label == "source_NP_nonexcluded"
                    annotation_rows.append({
                        "dataset": "GSE251686",
                        "gsm": gsm,
                        "presumed_sample_library_key": gsm,
                        "compartment_source": "NP",
                        "severity_group": row["severity_group"],
                        "barcode": barcode,
                        "total_umi": totals[cell],
                        "detected_genes": detected_genes[cell],
                        "pct_mt": f"{mt_pct:.6f}",
                        "qc_pass": not qc_reasons,
                        "qc_reason": ";".join(qc_reasons) if qc_reasons else "pass",
                        "annotation_label": labels[cell],
                        "compartment_pseudobulk_include": include[cell],
                    })

                retained_cells = sum(include)
                selected_total_umi = 0
                module_counts = {module: {gene: 0 for gene in genes} for module, genes in modules.items()}
                module_feature_map: dict[int, list[tuple[str, str]]] = defaultdict(list)
                for feature_index, symbol in enumerate(symbols):
                    for module, genes in modules.items():
                        if symbol in genes:
                            module_feature_map[feature_index].append((module, symbol))
                for feature, cell, value in iter_matrix(matrix_path, n_features, n_cells, expected_nnz):
                    if not include[cell] or value <= 0:
                        continue
                    selected_total_umi += value
                    for module, symbol in module_feature_map.get(feature, []):
                        module_counts[module][symbol] += value

                if retained_cells < args.min_retained_cells:
                    raise ValueError(
                        f"{gsm}: {retained_cells} source-restricted cells below the "
                        f"{args.min_retained_cells}-cell eligibility threshold"
                    )
                for module, genes in modules.items():
                    mapped = [gene for gene in genes if gene in symbol_counts]
                    feature_ids_by_symbol = {
                        gene: [feature_ids[index] for index, symbol in enumerate(symbols) if symbol == gene]
                        for gene in mapped
                    }
                    fraction = len(mapped) / len(genes)
                    duplicate_rows = sum(max(0, symbol_counts[gene] - 1) for gene in mapped)
                    mapping_rows.append({
                        "dataset": "GSE251686",
                        "gsm": gsm,
                        "module_id": module,
                        "configured_genes": len(genes),
                        "mapped_genes": len(mapped),
                        "mapped_fraction": f"{fraction:.6f}",
                        "mapping_pass": fraction >= args.min_mapped_fraction,
                        "mapping_key": "Ensembl feature ID with explicit uppercase symbol audit",
                        "mapped_gene_symbols": ";".join(mapped),
                        "mapped_ensembl_feature_ids": ";".join(
                            f"{gene}:{'|'.join(feature_ids_by_symbol[gene])}" for gene in mapped
                        ),
                        "duplicate_symbol_rows_in_module": duplicate_rows,
                        "duplicate_symbol_policy": "sum all Ensembl feature rows sharing the symbol; retain audit",
                    })
                    for gene in genes:
                        gene_rows.append({
                            "dataset": "GSE251686",
                            "gsm": gsm,
                            "module_id": module,
                            "gene_symbol": gene,
                            "feature_rows_mapped": symbol_counts.get(gene, 0),
                            "pseudobulk_count": module_counts[module][gene],
                            "total_umi_included_cells": selected_total_umi,
                        })
                    if fraction < args.min_mapped_fraction or selected_total_umi <= 0:
                        score = float("nan")
                        status = "mapping_below_minimum_or_zero_selected_umi"
                    else:
                        score = sum(
                            math.log1p(1_000_000.0 * module_counts[module][gene] / selected_total_umi)
                            for gene in mapped
                        ) / len(mapped)
                        status = "score_available"
                    score_rows.append({
                        "dataset": "GSE251686",
                        "gsm": gsm,
                        "presumed_sample_library_key": gsm,
                        "compartment": "NP",
                        "severity_group": row["severity_group"],
                        "module_id": module,
                        "module_score_log1p_cpm": "" if math.isnan(score) else f"{score:.8f}",
                        "score_status": status,
                        "mapped_fraction": f"{fraction:.6f}",
                        "included_cells": retained_cells,
                        "total_umi_included_cells": selected_total_umi,
                        "analysis_role": "incomplete non-balanced exploratory NP severity direction check only",
                        "confirmatory_eligible": "false",
                    })
                library_rows.append({
                    "dataset": "GSE251686",
                    "gsm": gsm,
                    "presumed_sample_library_key": gsm,
                    "compartment": "NP",
                    "severity_group": row["severity_group"],
                    "pre_qc_cells": n_cells,
                    "qc_passing_cells": sum(label != "qc_excluded" for label in labels),
                    "source_restricted_cells": retained_cells,
                    "source_restricted_threshold_20_pass": retained_cells >= 20,
                    "source_restricted_threshold_30_pass": retained_cells >= 30,
                    "source_restricted_threshold_50_pass": retained_cells >= 50,
                    "selected_umi": selected_total_umi,
                    "stream_integrity_pass": "true",
                    "identifier_audit_pass": "true",
                    "duplicate_gene_symbol_rows": sum(max(0, count - 1) for count in symbol_counts.values()),
                    "analysis_role": "incomplete non-balanced exploratory NP severity direction check only",
                    "confirmatory_eligible": "false",
                })

    write_csv(output_dir / "GSE251686_exploratory_module_scores.csv", score_rows)
    write_csv(output_dir / "GSE251686_exploratory_module_mapping_audit.csv", mapping_rows)
    write_csv(output_dir / "GSE251686_exploratory_module_gene_scores.csv", gene_rows)
    write_csv(output_dir / "GSE251686_exploratory_library_ledger.csv", library_rows)
    with gzip.open(output_dir / "GSE251686_exploratory_cell_annotation.csv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(annotation_rows[0]))
        writer.writeheader()
        writer.writerows(annotation_rows)
    parameter_rows = [
        {"parameter": "analysis_script", "value": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
        {"parameter": "archive", "value": str(args.archive.resolve()), "sha256": sha256(args.archive)},
        {"parameter": "ledger", "value": str(args.ledger.resolve()), "sha256": sha256(args.ledger)},
        {"parameter": "stream_audit", "value": str(args.stream_audit.resolve()), "sha256": sha256(args.stream_audit)},
        {"parameter": "nested_audit", "value": str(args.nested_audit.resolve()), "sha256": sha256(args.nested_audit)},
        {"parameter": "identifier_audit", "value": str(args.identifier_audit.resolve()), "sha256": sha256(args.identifier_audit)},
        {"parameter": "module_config", "value": str(args.module_config.resolve()), "sha256": sha256(args.module_config)},
        {"parameter": "panel_config", "value": str(args.panel_config.resolve()), "sha256": sha256(args.panel_config)},
        {"parameter": "selected_gsms", "value": ";".join(selected_gsms), "sha256": ""},
        {"parameter": "excluded_gsms", "value": BAD_GSM, "sha256": ""},
        {"parameter": "min_genes", "value": str(args.min_genes), "sha256": ""},
        {"parameter": "min_umi", "value": str(args.min_umi), "sha256": ""},
        {"parameter": "max_mt_pct", "value": str(args.max_mt_pct), "sha256": ""},
        {"parameter": "min_retained_cells", "value": str(args.min_retained_cells), "sha256": ""},
        {"parameter": "min_mapped_fraction", "value": str(args.min_mapped_fraction), "sha256": ""},
        {"parameter": "formula", "value": "mean over mapped symbols of log1p(1e6 * selected-cell pseudobulk symbol count / selected-cell total UMI)", "sha256": ""},
        {"parameter": "duplicate_feature_rule", "value": "Ensembl IDs are primary; all feature rows sharing an uppercase symbol are summed and counted in the mapping audit", "sha256": ""},
        {"parameter": "inference_boundary", "value": "presumed sample/library key; cells nested; exploratory direction display only; never confirmatory", "sha256": ""},
    ]
    write_csv(output_dir / "GSE251686_exploratory_score_parameters.csv", parameter_rows)
    (output_dir / "README.md").write_text(
        "# GSE251686 exploratory module scores\n\n"
        "This output was generated only from the five libraries whose independent "
        "stream-integrity and identifier audits passed. `GSM7986002` was never "
        "opened by the scorer and remains permanently excluded. Ensembl feature "
        "IDs are the primary cross-library key; module symbols are mapped with an "
        "explicit duplicate-symbol audit, and duplicate feature rows are summed.\n\n"
        "The result is an incomplete, non-balanced mild n=2 versus severe n=3 "
        "presumed sample/library direction check. It is not a validation cohort, "
        "does not enter the default 20-effect summary, and cannot support cell-level, "
        "causal, therapeutic, or universal-program claims.\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
