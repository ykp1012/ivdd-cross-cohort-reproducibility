"""Audit the native clinical-comparison subset of GSE245147.

GSE245147 contains a clean three-versus-three Degenerated/NO_Degenerated NP
comparison plus separate P2/P8 passage and DMSO/H-151 treatment arms.  This
script records the explicit six-sample selection and proves that the other
arms are excluded before score-level analysis.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
from pathlib import Path


BIOPROJECT = re.compile(r"BioProject:\s*https?://[^/]+/bioproject/([^\s\"]+)", re.IGNORECASE)
BIOSAMPLE = re.compile(r"BioSample:\s*https?://[^/]+/biosample/([^\s\"]+)", re.IGNORECASE)
NATIVE_TITLE = re.compile(r"^(Degenerated|No-degenerated) nucleus pulposus cells #[123] for RNA-seq$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def read_modules(path: Path) -> dict[str, list[str]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return {
        str(module["module_id"]): [str(gene).upper() for gene in module["genes"]]
        for module in payload["modules"]
    }


def main() -> int:
    parser = __import__("argparse").ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    raw_dir = root / "data" / "raw" / "geo_candidates" / "GSE245147"
    ledger_path = root / "data" / "derived" / "geo_candidate_audit" / "GSE245147_sample_ledger.csv"
    matrix_path = raw_dir / "GSE245147_Degenerated_NO_Degenerated_RPKM.txt.gz"
    soft_path = raw_dir / "GSE245147_family.soft.gz"
    module_path = root / "config" / "program_modules.json"
    required = [ledger_path, matrix_path, soft_path, module_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required GSE245147 audit input(s):\n" + "\n".join(missing))

    with ledger_path.open(newline="", encoding="utf-8-sig") as handle:
        ledger = list(csv.DictReader(handle))
    if len(ledger) != 18:
        raise ValueError(f"Expected 18 GSE245147 ledger rows, found {len(ledger)}")

    selection_rows: list[dict[str, object]] = []
    selected_gsm: list[str] = []
    for row in ledger:
        title = row["sample_title"].strip()
        native = bool(NATIVE_TITLE.fullmatch(title))
        reason = "native_Degenerated_vs_NO_Degenerated_clinical_comparison" if native else (
            "excluded_passage_arm_P2_or_P8" if re.match(r"^P[28] nucleus", title) and "treated" not in title else
            "excluded_treatment_arm_DMSO_or_H-151"
        )
        if native:
            selected_gsm.append(row["gsm"].strip())
        selection_rows.append({
            "dataset": "GSE245147",
            "gsm": row["gsm"].strip(),
            "sample_title": title,
            "selected_for_meta_analysis": str(native).lower(),
            "selection_reason": reason,
            "source_disease_state": row.get("disease_state", ""),
            "tissue": row.get("tissue", ""),
        })
    if len(selected_gsm) != 6 or sum("Degenerated" in row["sample_title"] for row in selection_rows if row["selected_for_meta_analysis"] == "true") != 3:
        raise ValueError(f"Native subset is not exactly three versus three: {selected_gsm}")
    selected_columns = [
        "Degenerated_1", "Degenerated_2", "Degenerated_3",
        "NO_Degenerated_1", "NO_Degenerated_2", "NO_Degenerated_3",
    ]

    modules = read_modules(module_path)
    observed_symbols: set[str] = set()
    feature_rows = 0
    malformed_rows = 0
    nonfinite_values = 0
    negative_values = 0
    matrix_header: list[str] = []
    selected_column_indexes: list[int] = []
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        first = handle.readline().rstrip("\r\n")
        matrix_header = first.split("\t")
        expected_prefix = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
        if matrix_header[: len(expected_prefix)] != expected_prefix:
            raise ValueError(f"Unexpected GSE245147 matrix prefix: {matrix_header[:6]}")
        missing_columns = [column for column in selected_columns if column not in matrix_header]
        if missing_columns:
            raise ValueError(f"Missing native matrix column(s): {missing_columns}")
        selected_column_indexes = [matrix_header.index(column) for column in selected_columns]
        for raw_line in handle:
            fields = raw_line.rstrip("\r\n").split("\t")
            if not fields:
                continue
            if len(fields) != len(matrix_header):
                malformed_rows += 1
                continue
            feature = fields[0].strip().upper()
            if feature:
                observed_symbols.add(feature)
            feature_rows += 1
            for index in selected_column_indexes:
                try:
                    value = float(fields[index])
                except ValueError:
                    nonfinite_values += 1
                    continue
                if not math.isfinite(value):
                    nonfinite_values += 1
                elif value < 0:
                    negative_values += 1

    mapping_rows: list[dict[str, object]] = []
    for module_id, genes in modules.items():
        mapped = [gene for gene in genes if gene in observed_symbols]
        missing = [gene for gene in genes if gene not in observed_symbols]
        fraction = len(mapped) / len(genes)
        mapping_rows.append({
            "dataset": "GSE245147",
            "subset": "Degenerated_1-3_vs_NO_Degenerated_1-3",
            "module_id": module_id,
            "configured_genes": len(genes),
            "mapped_genes": len(mapped),
            "mapped_fraction": f"{fraction:.6f}",
            "mapping_pass_at_0_80": str(fraction >= 0.80).lower(),
            "mapped_gene_symbols": ";".join(mapped),
            "missing_gene_symbols": ";".join(missing),
            "annotation_method": "Exact uppercase match to Geneid field; TNF is absent from this RPKM matrix",
        })

    bioprojects: set[str] = set()
    biosamples: set[str] = set()
    with gzip.open(soft_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            project_match = BIOPROJECT.search(line)
            if project_match:
                bioprojects.add(project_match.group(1))
            sample_match = BIOSAMPLE.search(line)
            if sample_match:
                biosamples.add(sample_match.group(1))
    if bioprojects != {"PRJNA1027236"}:
        raise ValueError(f"Unexpected GSE245147 BioProject set: {sorted(bioprojects)}")
    independence = {
        "dataset": "GSE245147",
        "bioproject": "PRJNA1027236",
        "pubmed": "38488012",
        "unique_biosamples_exposed": len(biosamples),
        "native_comparison_biosamples": len(selected_gsm),
        "default_cohort_bioproject_overlap": "not observed at local SOFT/BioProject level",
        "patient_level_independence": "not verifiable: public GEO metadata expose no patient IDs, age, sex, or disc level",
        "related_candidate_boundary": "GSE266883 is excluded because the human design and author chain are highly similar; do not pool both",
        "analysis_boundary": "score-level RPKM only; no raw-count model, treatment effect, passage effect, causal or confirmatory inference",
    }

    matrix_audit = {
        "dataset": "GSE245147",
        "subset": "Degenerated_1-3_vs_NO_Degenerated_1-3",
        "matrix_file": str(matrix_path.relative_to(root)).replace("\\", "/"),
        "matrix_sha256": sha256(matrix_path),
        "matrix_bytes": matrix_path.stat().st_size,
        "matrix_header_columns": len(matrix_header),
        "selected_sample_columns": ";".join(selected_columns),
        "excluded_sample_arm_count": 12,
        "feature_rows": feature_rows,
        "unique_geneid_values": len(observed_symbols),
        "malformed_rows": malformed_rows,
        "selected_values_finite": str(nonfinite_values == 0).lower(),
        "selected_values_nonnegative": str(negative_values == 0).lower(),
        "module_mapping_all_pass": str(all(row["mapping_pass_at_0_80"] == "true" for row in mapping_rows)).lower(),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "sample_selection.csv", selection_rows, list(selection_rows[0]))
    write_csv(output_dir / "matrix_integrity.csv", [matrix_audit], list(matrix_audit))
    write_csv(output_dir / "module_mapping.csv", mapping_rows, list(mapping_rows[0]))
    (output_dir / "independence_ledger.json").write_text(json.dumps(independence, indent=2) + "\n", encoding="utf-8")

    generated = [
        output_dir / "sample_selection.csv",
        output_dir / "matrix_integrity.csv",
        output_dir / "module_mapping.csv",
        output_dir / "independence_ledger.json",
    ]
    manifest = {
        "schema_version": 1,
        "purpose": "Audit the selected native clinical-comparison subset before score-level inclusion.",
        "selected_gsm": selected_gsm,
        "excluded_design_arms": "P2/P8 passage and DMSO/H-151 treatment samples",
        "input_sha256": {
            str(path.relative_to(root)).replace("\\", "/"): sha256(path)
            for path in [ledger_path, matrix_path, soft_path, module_path]
        },
        "generated_artifact_sha256": {
            path.name: sha256(path) for path in generated
        },
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"GSE245147 native-comparison audit completed: {output_dir / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
