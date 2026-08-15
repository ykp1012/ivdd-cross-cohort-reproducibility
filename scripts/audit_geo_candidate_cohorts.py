"""Audit candidate human NP GEO cohorts without modifying raw inputs.

The script produces a compact, reproducible evidence package for three
candidate cohorts.  It checks the exact sample-column mapping, scans every
processed matrix for basic structural integrity, audits locked-module feature
coverage where gene symbols are available, and records the accession-level
independence boundary.  It intentionally does not calculate differential
expression or alter the frozen default IVDD analysis.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Iterable


ENSEMBL_SYMBOL = re.compile(r"^(ENSG\d+(?:\.\d+)?)\((.+)\)$")
BIOSAMPLE = re.compile(r"BioSample:\s*https?://[^/]+/biosample/([^\s\"]+)", re.IGNORECASE)
BIOPROJECT = re.compile(r"BioProject:\s*https?://[^/]+/bioproject/([^\s\"]+)", re.IGNORECASE)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_modules(path: Path) -> dict[str, list[str]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return {
        str(module["module_id"]): [str(gene).upper() for gene in module["genes"]]
        for module in payload["modules"]
    }


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_gzip_lines(path: Path) -> Iterable[str]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        yield from handle


def write_csv_new(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing audit artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def unquote(value: str) -> str:
    return value.strip().strip('"')


def parse_soft_relations(path: Path) -> tuple[set[str], set[str], list[str]]:
    biosamples: set[str] = set()
    bioprojects: set[str] = set()
    identical_rna: list[str] = []
    for raw in read_gzip_lines(path):
        line = raw.rstrip("\r\n")
        match = BIOSAMPLE.search(line)
        if match:
            biosamples.add(match.group(1))
        match = BIOPROJECT.search(line)
        if match:
            bioprojects.add(match.group(1))
        if "identical RNA" in line:
            identical_rna.append(line)
    return biosamples, bioprojects, identical_rna


def audit_tabular_matrix(
    path: Path,
    *,
    expected_columns: list[str],
    feature_field: str,
    gene_symbol: callable | None = None,
    stop_marker: str | None = None,
    allow_column_reorder: bool = False,
    numeric_start_index: int = 1,
) -> dict[str, object]:
    """Stream a gzip TSV, retaining only integrity and identifier summaries."""
    lines = iter(read_gzip_lines(path))
    if stop_marker is not None:
        for raw in lines:
            if raw.rstrip("\r\n") == stop_marker:
                break
        else:
            raise ValueError(f"Could not find table marker {stop_marker!r} in {path}")

    try:
        header = next(lines).rstrip("\r\n").split("\t")
    except StopIteration as exc:
        raise ValueError(f"No table header in {path}") from exc
    header = [unquote(value) for value in header]
    if not header or header[0] != feature_field:
        raise ValueError(f"{path}: expected first field {feature_field!r}, got {header[:1]!r}")
    header_columns = header[1:]
    columns_match = header_columns == expected_columns
    if allow_column_reorder:
        columns_match = len(header_columns) == len(expected_columns) and set(header_columns) == set(expected_columns)
    if not columns_match:
        raise ValueError(
            f"{path}: sample columns differ from expected. got={header_columns!r}, expected={expected_columns!r}"
        )

    feature_ids: set[str] = set()
    duplicate_feature_ids = 0
    malformed_rows = 0
    non_finite_values = 0
    negative_values = 0
    non_integer_values = 0
    feature_rows = 0
    symbols: set[str] = set()
    for raw in lines:
        line = raw.rstrip("\r\n")
        if not line:
            continue
        if stop_marker is not None and line == "!series_matrix_table_end":
            break
        values = line.split("\t")
        # GEO's GSE167931 text tables retain a terminal tab on data rows but
        # not on the header.  It is formatting, not an extra sample column.
        if len(values) == len(header) + 1 and values[-1] == "":
            values.pop()
        if len(values) != len(header):
            malformed_rows += 1
            continue
        feature_id = unquote(values[0])
        if feature_id in feature_ids:
            duplicate_feature_ids += 1
        feature_ids.add(feature_id)
        if gene_symbol is not None:
            symbol = gene_symbol(values)
            if symbol:
                symbols.add(symbol.upper())
        feature_rows += 1
        for value in values[numeric_start_index:]:
            try:
                numeric = float(value)
            except ValueError:
                non_finite_values += 1
                continue
            if not math.isfinite(numeric):
                non_finite_values += 1
            elif numeric < 0:
                negative_values += 1
            elif not numeric.is_integer():
                non_integer_values += 1
    return {
        "header": header,
        "feature_rows": feature_rows,
        "unique_feature_ids": len(feature_ids),
        "duplicate_feature_ids": duplicate_feature_ids,
        "malformed_rows": malformed_rows,
        "non_finite_values": non_finite_values,
        "negative_values": negative_values,
        "non_integer_values": non_integer_values,
        "gene_symbols": symbols,
    }


def module_mapping_rows(
    dataset: str,
    modules: dict[str, list[str]],
    symbols: set[str] | None,
    *,
    annotation_method: str,
    annotation_status: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for module_id, genes in modules.items():
        if symbols is None:
            rows.append({
                "dataset": dataset,
                "module_id": module_id,
                "configured_genes": len(genes),
                "mapped_genes": "",
                "mapped_fraction": "",
                "mapping_status": annotation_status,
                "mapping_pass_at_0_80": "",
                "mapped_gene_symbols": "",
                "missing_gene_symbols": ";".join(genes),
                "annotation_method": annotation_method,
            })
            continue
        mapped = [gene for gene in genes if gene in symbols]
        missing = [gene for gene in genes if gene not in symbols]
        fraction = len(mapped) / len(genes)
        rows.append({
            "dataset": dataset,
            "module_id": module_id,
            "configured_genes": len(genes),
            "mapped_genes": len(mapped),
            "mapped_fraction": f"{fraction:.6f}",
            "mapping_status": annotation_status,
            "mapping_pass_at_0_80": str(fraction >= 0.80).lower(),
            "mapped_gene_symbols": ";".join(mapped),
            "missing_gene_symbols": ";".join(missing),
            "annotation_method": annotation_method,
        })
    return rows


def expected_sample_rows(
    dataset: str,
    matrix_file: str,
    ledger: list[dict[str, str]],
    column_builder: callable,
    group_builder: callable,
    contrast: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for sample in ledger:
        title = sample["sample_title"]
        rows.append({
            "dataset": dataset,
            "matrix_file": matrix_file,
            "gsm": sample["gsm"],
            "sample_title": title,
            "matrix_column": column_builder(sample),
            "comparison_group": group_builder(sample),
            "tissue": sample["tissue"],
            "recorded_disease_state": sample["disease_state"],
            "recorded_grade": sample["degeneration_grade"],
            "contrast_interpretation": contrast,
            "mapping_rule": "Exact header-to-ledger mapping after documented title normalization",
            "mapping_status": "pass",
            "identity_unit": "GSM/sample key; no patient identifier exposed in GEO metadata",
        })
    return rows


def find_default_bioprojects(project_root: Path) -> set[str]:
    paths = [
        project_root / "data" / "raw" / "GSE229711_family.soft.gz",
        project_root / "data" / "raw" / "GSE230808_family.soft.gz",
        project_root / "data" / "raw" / "GSE153066_family.soft.gz",
        project_root / "data" / "raw" / "GSE244889_family.soft.gz",
        project_root / "data" / "raw" / "GSE165722_family.soft.gz",
    ]
    result: set[str] = set()
    for path in paths:
        if path.exists():
            _, projects, _ = parse_soft_relations(path)
            result.update(projects)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    raw = root / "data" / "raw" / "geo_candidates"
    ledgers = output_dir
    modules = read_modules(root / "config" / "program_modules.json")

    files = {
        "GSE186542": {
            "matrix": raw / "GSE186542" / "GSE186542_gene_expression.txt.gz",
            "soft": raw / "GSE186542" / "GSE186542_family.soft.gz",
            "series": raw / "GSE186542" / "GSE186542_series_matrix.txt.gz",
            "ledger": ledgers / "GSE186542_sample_ledger.csv",
        },
        "GSE167931": {
            "tpm": raw / "GSE167931" / "GSE167931_AllSamplesTPMValue.txt.gz",
            "fpkm": raw / "GSE167931" / "GSE167931_AllSamplesFPKMValue.txt.gz",
            "soft": raw / "GSE167931" / "GSE167931_family.soft.gz",
            "series": raw / "GSE167931" / "GSE167931_series_matrix.txt.gz",
            "ledger": ledgers / "GSE167931_sample_ledger.csv",
        },
        "GSE56081": {
            "matrix": raw / "GSE56081" / "GSE56081_series_matrix.txt.gz",
            "soft": raw / "GSE56081" / "GSE56081_family.soft.gz",
            "ledger": ledgers / "GSE56081_sample_ledger.csv",
        },
    }
    missing = [str(path) for item in files.values() for path in item.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required candidate inputs:\n" + "\n".join(missing))

    ledger_186542 = read_csv_rows(files["GSE186542"]["ledger"])
    ledger_167931 = read_csv_rows(files["GSE167931"]["ledger"])
    ledger_56081 = read_csv_rows(files["GSE56081"]["ledger"])

    columns_186542 = [f"{row['sample_title']}_count" for row in ledger_186542]
    matrix_186542 = audit_tabular_matrix(
        files["GSE186542"]["matrix"],
        expected_columns=["gene_name", "description", "locus", *columns_186542],
        feature_field="gene_id",
        gene_symbol=lambda values: values[1].strip(),
        allow_column_reorder=True,
        numeric_start_index=4,
    )
    # The source supplementary matrix orders columns differently from the SOFT
    # sample order; same set is the correct identity contract.
    if set(matrix_186542["header"][4:]) != set(columns_186542):
        raise ValueError("GSE186542 matrix does not contain exactly the six expected count columns")

    columns_167931 = [row["sample_title"].replace("-", "_") for row in ledger_167931]

    def gse167931_symbol(values: list[str]) -> str:
        match = ENSEMBL_SYMBOL.match(values[0].strip())
        return match.group(2).strip() if match else ""

    matrix_167931_tpm = audit_tabular_matrix(
        files["GSE167931"]["tpm"],
        expected_columns=columns_167931,
        feature_field="GeneID",
        gene_symbol=gse167931_symbol,
    )
    matrix_167931_fpkm = audit_tabular_matrix(
        files["GSE167931"]["fpkm"],
        expected_columns=columns_167931,
        feature_field="GeneID",
        gene_symbol=gse167931_symbol,
    )

    columns_56081 = [row["gsm"] for row in ledger_56081]
    matrix_56081 = audit_tabular_matrix(
        files["GSE56081"]["matrix"],
        expected_columns=columns_56081,
        feature_field="ID_REF",
        stop_marker="!series_matrix_table_begin",
    )

    matrix_rows = [
        {
            "dataset": "GSE186542",
            "matrix_file": files["GSE186542"]["matrix"].relative_to(root).as_posix(),
            "sha256": sha256(files["GSE186542"]["matrix"]),
            "feature_field": "gene_id; gene_name",
            "feature_rows": matrix_186542["feature_rows"],
            "sample_columns": len(columns_186542),
            "sample_column_names": ";".join(matrix_186542["header"][4:]),
            "header_mapping": "pass: set of <sample_title>_count columns exactly matches all six ledger titles; file order differs from SOFT order",
            "numeric_integrity": "pass" if not any(matrix_186542[key] for key in ["malformed_rows", "non_finite_values", "negative_values", "non_integer_values"]) else "review",
            "duplicate_feature_ids": matrix_186542["duplicate_feature_ids"],
            "malformed_rows": matrix_186542["malformed_rows"],
            "non_finite_values": matrix_186542["non_finite_values"],
            "negative_values": matrix_186542["negative_values"],
            "non_integer_values": matrix_186542["non_integer_values"],
            "value_scale": "integer count-like values; GEO also states FPKM/TMM processing",
            "raw_count_model_eligibility": "not yet: count provenance is internally inconsistent and needs SRA/recount confirmation",
        },
        {
            "dataset": "GSE167931",
            "matrix_file": files["GSE167931"]["tpm"].relative_to(root).as_posix(),
            "sha256": sha256(files["GSE167931"]["tpm"]),
            "feature_field": "GeneID formatted as ENSG...(SYMBOL)",
            "feature_rows": matrix_167931_tpm["feature_rows"],
            "sample_columns": len(columns_167931),
            "sample_column_names": ";".join(matrix_167931_tpm["header"][1:]),
            "header_mapping": "pass: hyphen-to-underscore normalization of all nine SOFT sample titles",
            "numeric_integrity": "pass" if not any(matrix_167931_tpm[key] for key in ["malformed_rows", "non_finite_values", "negative_values"]) else "review",
            "duplicate_feature_ids": matrix_167931_tpm["duplicate_feature_ids"],
            "malformed_rows": matrix_167931_tpm["malformed_rows"],
            "non_finite_values": matrix_167931_tpm["non_finite_values"],
            "negative_values": matrix_167931_tpm["negative_values"],
            "non_integer_values": matrix_167931_tpm["non_integer_values"],
            "value_scale": "TPM; normalized expression",
            "raw_count_model_eligibility": "no: use log2(TPM + 1) score-level or limma-style analysis only",
        },
        {
            "dataset": "GSE167931",
            "matrix_file": files["GSE167931"]["fpkm"].relative_to(root).as_posix(),
            "sha256": sha256(files["GSE167931"]["fpkm"]),
            "feature_field": "GeneID formatted as ENSG...(SYMBOL)",
            "feature_rows": matrix_167931_fpkm["feature_rows"],
            "sample_columns": len(columns_167931),
            "sample_column_names": ";".join(matrix_167931_fpkm["header"][1:]),
            "header_mapping": "pass: hyphen-to-underscore normalization of all nine SOFT sample titles",
            "numeric_integrity": "pass" if not any(matrix_167931_fpkm[key] for key in ["malformed_rows", "non_finite_values", "negative_values"]) else "review",
            "duplicate_feature_ids": matrix_167931_fpkm["duplicate_feature_ids"],
            "malformed_rows": matrix_167931_fpkm["malformed_rows"],
            "non_finite_values": matrix_167931_fpkm["non_finite_values"],
            "negative_values": matrix_167931_fpkm["negative_values"],
            "non_integer_values": matrix_167931_fpkm["non_integer_values"],
            "value_scale": "FPKM/RPKM-family normalized expression",
            "raw_count_model_eligibility": "no: use log2(FPKM + 1) score-level or limma-style analysis only",
        },
        {
            "dataset": "GSE56081",
            "matrix_file": files["GSE56081"]["matrix"].relative_to(root).as_posix(),
            "sha256": sha256(files["GSE56081"]["matrix"]),
            "feature_field": "ID_REF Arraystar probe ID",
            "feature_rows": matrix_56081["feature_rows"],
            "sample_columns": len(columns_56081),
            "sample_column_names": ";".join(matrix_56081["header"][1:]),
            "header_mapping": "pass: exact GSM-to-matrix-column match for all ten samples",
            "numeric_integrity": "pass" if not any(matrix_56081[key] for key in ["malformed_rows", "non_finite_values", "negative_values"]) else "review",
            "duplicate_feature_ids": matrix_56081["duplicate_feature_ids"],
            "malformed_rows": matrix_56081["malformed_rows"],
            "non_finite_values": matrix_56081["non_finite_values"],
            "negative_values": matrix_56081["negative_values"],
            "non_integer_values": matrix_56081["non_integer_values"],
            "value_scale": "log2 quantile-normalized microarray intensity",
            "raw_count_model_eligibility": "no: microarray data; platform-specific probe-to-gene annotation is required before module analysis",
        },
    ]

    sample_rows = []
    sample_rows.extend(expected_sample_rows(
        "GSE186542", "GSE186542_gene_expression.txt.gz", ledger_186542,
        lambda row: f"{row['sample_title']}_count",
        lambda row: "early_stage_Pfirrmann_I_to_III" if row["disease_state"] == "control" else "advanced_stage_Pfirrmann_IV_to_V",
        "advanced-stage minus early-stage degeneration; series design is the authority for Pfirrmann bands",
    ))
    sample_rows.extend(expected_sample_rows(
        "GSE167931", "GSE167931_AllSamplesTPMValue.txt.gz", ledger_167931,
        lambda row: row["sample_title"].replace("-", "_"),
        lambda row: "normal" if row["sample_title"].startswith("NP-C") else "degenerated",
        "degenerated NP cells minus normal NP cells",
    ))
    sample_rows.extend(expected_sample_rows(
        "GSE56081", "GSE56081_series_matrix.txt.gz", ledger_56081,
        lambda row: row["gsm"],
        lambda row: "control_grade_I" if "control" in row["sample_title"].lower() else "degenerated_grade_IV_to_V",
        "degenerated grade IV/V NP tissue minus control grade I NP tissue",
    ))

    mapping_rows = []
    mapping_rows.extend(module_mapping_rows(
        "GSE186542", modules, matrix_186542["gene_symbols"],
        annotation_method="Exact uppercase match to supplementary gene_name column",
        annotation_status="evaluated",
    ))
    mapping_rows.extend(module_mapping_rows(
        "GSE167931", modules, matrix_167931_tpm["gene_symbols"],
        annotation_method="Exact uppercase match to symbol parsed from ENSG...(SYMBOL) GeneID",
        annotation_status="evaluated",
    ))
    mapping_rows.extend(module_mapping_rows(
        "GSE56081", modules, None,
        annotation_method="Not attempted: local GPL15314 table exposes custom probe IDs and blank ENSEMBL_ID/ACCESSION_STRING fields for the inspected rows",
        annotation_status="unresolved_platform_annotation",
    ))

    default_projects = find_default_bioprojects(root)
    independence_rows = []
    for dataset, label, parent, role, note in [
        (
            "GSE186542", "PRJNA774342", "none stated in GEO series metadata", "conditional external NP cohort",
            "Six unique BioSample accessions are exposed; no patient ID, age, sex, or disc level is exposed. Independent at BioProject/accession level only.",
        ),
        (
            "GSE167931", "PRJNA705603", "none stated in GEO series metadata", "conditional external NP cohort",
            "Nine unique BioSample accessions are exposed; no patient ID, age, sex, or disc level is exposed. Independent at BioProject/accession level only; literature-level donor overlap remains untestable from GEO metadata.",
        ),
        (
            "GSE56081", "PRJNA242356", "SubSeries of GSE67567", "secondary microarray cohort pending annotation",
            "SOFT identifies identical-RNA alternatives for every GSM. Count this biological sample set once and do not treat its alternative subseries as independent cohorts.",
        ),
    ]:
        biosamples, projects, identical_rna = parse_soft_relations(files[dataset]["soft"])
        if projects != {label}:
            raise ValueError(f"{dataset}: expected BioProject {label}, found {sorted(projects)}")
        independence_rows.append({
            "dataset": dataset,
            "bioproject": label,
            "bioproject_distinct_from_default_project_set": str(label not in default_projects).lower(),
            "default_project_set_seen_in_local_SOFT": ";".join(sorted(default_projects)),
            "parent_or_superseries_relation": parent,
            "biosample_count": len(biosamples),
            "biosample_ids": ";".join(sorted(biosamples)),
            "identical_rna_alternatives_reported": len(identical_rna),
            "independence_assessment": "provisionally independent from current default cohorts" if label not in default_projects else "not independent from current default cohorts",
            "role": role,
            "limitation": note,
        })

    assessment_rows = [
        {
            "dataset": "GSE186542",
            "species_tissue": "Homo sapiens nucleus pulposus tissue",
            "assay": "bulk RNA-seq",
            "comparison": "early-stage Pfirrmann I-III (n=3) versus advanced-stage IV-V (n=3)",
            "processed_matrix": "GSE186542_gene_expression.txt.gz; six count-labelled columns",
            "module_mapping": "all four locked modules pass 80% coverage (100% observed)",
            "independence": "distinct BioProject PRJNA774342; no stated reuse in its GEO series metadata",
            "recommendation": "conditional include as a small external score-level comparison; do not call the supplied values raw-count-model eligible until the FPKM/TMM versus raw-count provenance conflict is resolved",
            "analysis_boundary": "advanced-stage minus early-stage; not a healthy-versus-diseased contrast; donor metadata are incomplete",
        },
        {
            "dataset": "GSE167931",
            "species_tissue": "Homo sapiens nucleus pulposus cells",
            "assay": "bulk RNA-seq of isolated NP cells",
            "comparison": "normal (n=4) versus degenerated (n=5)",
            "processed_matrix": "paired TPM and FPKM matrices; nine matching sample columns",
            "module_mapping": "all four locked modules pass 80% coverage (100% observed)",
            "independence": "distinct BioProject PRJNA705603; no stated reuse in its GEO series metadata",
            "recommendation": "conditional include as an external normalized-expression score cohort using a predeclared log2(TPM + 1) or log2(FPKM + 1) scale; do not use DESeq2/edgeR raw-count inference",
            "analysis_boundary": "degenerated minus normal; donor metadata and publication-level overlap checks remain incomplete",
        },
        {
            "dataset": "GSE56081",
            "species_tissue": "Homo sapiens nucleus pulposus tissue",
            "assay": "Arraystar Human LncRNA microarray V2.0",
            "comparison": "control grade I (n=5) versus degenerated grade IV/V (n=5)",
            "processed_matrix": "log2 quantile-normalized 10-column microarray matrix",
            "module_mapping": "unresolved: custom GPL15314 probe IDs require a defensible probe-to-gene annotation before locked-module scoring",
            "independence": "distinct BioProject PRJNA242356 but a GSE67567 subseries with identical-RNA alternatives",
            "recommendation": "retain as a secondary conventional microarray candidate; exclude from current locked-module extension until probe annotation is independently audited",
            "analysis_boundary": "count GSE56081 biological samples once; do not count identical-RNA alternatives as independent replication",
        },
    ]

    input_rows: list[dict[str, object]] = []
    for dataset, item in files.items():
        for kind, path in item.items():
            input_rows.append({
                "dataset": dataset,
                "input_kind": kind,
                "path": path.relative_to(root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            })

    write_csv_new(
        output_dir / "candidate_matrix_integrity_audit.csv", matrix_rows, list(matrix_rows[0]),
    )
    write_csv_new(
        output_dir / "candidate_sample_matrix_mapping.csv", sample_rows, list(sample_rows[0]),
    )
    write_csv_new(
        output_dir / "candidate_module_mapping_audit.csv", mapping_rows, list(mapping_rows[0]),
    )
    write_csv_new(
        output_dir / "candidate_independence_ledger.csv", independence_rows, list(independence_rows[0]),
    )
    write_csv_new(
        output_dir / "candidate_cohort_assessment.csv", assessment_rows, list(assessment_rows[0]),
    )
    write_csv_new(
        output_dir / "candidate_input_artifact_hashes.csv", input_rows, list(input_rows[0]),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
