"""Audit and score two external human NP bulk-RNA cohorts with locked modules.

The script only supports the explicitly audited GSE186542 count matrix and
GSE167931 FPKM/TPM matrices.  Each GEO sample is retained as a presumed
donor-level observation.  It deliberately does not infer donor identity beyond
the publicly supplied sample identifiers, perform gene-level differential
expression testing, or treat technical values as cells.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd


MIN_MAPPED_FRACTION = 0.80
MODULE_IDS = (
    "ecm_collagen_remodeling",
    "inflammatory_nfkb",
    "hypoxia_oxidative_stress",
    "disc_matrix_homeostasis",
)


@dataclass(frozen=True)
class CohortMatrix:
    cohort_id: str
    dataset: str
    relative_path: str
    gene_column: str
    gene_parser: Callable[[str], str]
    sample_to_gsm: dict[str, str]
    sample_to_group: dict[str, str]
    source_value_scale: str
    analysis_role: str
    matrix_status: str
    score_formula: str
    non_sample_columns: tuple[str, ...] = ()


def identity(value: str) -> str:
    return value.strip().upper()


def ensembl_parenthetical_symbol(value: str) -> str:
    match = re.search(r"\(([^()]+)\)\s*$", value)
    return match.group(1).strip().upper() if match else ""


GSE186542_COLUMNS = {
    "NP_control_1_count": "GSM5655437",
    "NP_control_2_count": "GSM5655438",
    "NP_control_3_count": "GSM5655439",
    "NP_degenerated_1_count": "GSM5655440",
    "NP_degenerated_2_count": "GSM5655441",
    "NP_degenerated_3_count": "GSM5655442",
}

GSE167931_COLUMNS = {
    "NP_C_3": "GSM5115822",
    "NP_C_4": "GSM5115823",
    "NP_C_7": "GSM5115824",
    "NP_C_8": "GSM5115825",
    "NP_E_1": "GSM5115826",
    "NP_E_2": "GSM5115827",
    "NP_E_3": "GSM5115828",
    "NP_E_4": "GSM5115829",
    "NP_E_5": "GSM5115830",
}

GSE245147_COLUMNS = {
    "Degenerated_1": "GSM7837568",
    "Degenerated_2": "GSM7837569",
    "Degenerated_3": "GSM7837570",
    "NO_Degenerated_1": "GSM7837571",
    "NO_Degenerated_2": "GSM7837572",
    "NO_Degenerated_3": "GSM7837573",
}

# GSE167931 reports two different Ensembl records with the symbol SOD2.  The
# canonical protein-coding record is retained; the alternate record is logged
# rather than summed into the locked module score.
CANONICAL_ENSEMBL_FOR_DUPLICATE_SYMBOL = {"SOD2": "ENSG00000112096"}


def group_from_control_or_degenerated(column: str) -> str:
    return "early_control" if "control" in column.lower() else "advanced_degenerated"


def group_from_normal_or_degenerated(column: str) -> str:
    return "normal" if "_C_" in column else "degenerated"


def group_from_native_degenerated_comparison(column: str) -> str:
    return "degenerated" if column.startswith("Degenerated_") else "no_degenerated"


def matrix_specs() -> tuple[CohortMatrix, ...]:
    return (
        CohortMatrix(
            cohort_id="GSE186542_external_count_support",
            dataset="GSE186542",
            relative_path=(
                "data/raw/geo_candidates/GSE186542/GSE186542_gene_expression.txt.gz"
            ),
            gene_column="gene_name",
            gene_parser=identity,
            sample_to_gsm=GSE186542_COLUMNS,
            sample_to_group={
                column: group_from_control_or_degenerated(column)
                for column in GSE186542_COLUMNS
            },
            source_value_scale="GEO-deposited raw gene counts",
            analysis_role="post_hoc_external_count_level_score_support_only",
            matrix_status=(
                "GEO supplementary gene-by-sample count matrix; six exact matrix-column-to-GSM "
                "matches; no publicly supplied patient ID, age, or sex"
            ),
            score_formula=(
                "mean over mapped genes of log1p(1e6 * raw gene count / supplied sample total count)"
            ),
        ),
        CohortMatrix(
            cohort_id="GSE167931_external_fpkm_support",
            dataset="GSE167931",
            relative_path=(
                "data/raw/geo_candidates/GSE167931/GSE167931_AllSamplesFPKMValue.txt.gz"
            ),
            gene_column="GeneID",
            gene_parser=ensembl_parenthetical_symbol,
            sample_to_gsm=GSE167931_COLUMNS,
            sample_to_group={
                column: group_from_normal_or_degenerated(column)
                for column in GSE167931_COLUMNS
            },
            source_value_scale="GEO-deposited FPKM",
            analysis_role="post_hoc_external_processed_score_support_only",
            matrix_status=(
                "GEO supplementary FPKM matrix; nine exact matrix-column-to-GSM matches; "
                "public series metadata has no patient ID, age, or sex"
            ),
            score_formula="mean over mapped genes of log1p(FPKM)",
        ),
        CohortMatrix(
            cohort_id="GSE245147_external_native_comparison_support",
            dataset="GSE245147",
            relative_path=(
                "data/raw/geo_candidates/GSE245147/"
                "GSE245147_Degenerated_NO_Degenerated_RPKM.txt.gz"
            ),
            gene_column="Geneid",
            gene_parser=identity,
            sample_to_gsm=GSE245147_COLUMNS,
            sample_to_group={
                column: group_from_native_degenerated_comparison(column)
                for column in GSE245147_COLUMNS
            },
            source_value_scale="GEO-deposited RPKM",
            analysis_role="post_hoc_external_processed_score_support_only_native_clinical_subset",
            matrix_status=(
                "GEO supplementary RPKM matrix; the six native Degenerated/NO_Degenerated columns map exactly "
                "to six GSM accessions; P2/P8 passage and DMSO/H-151 treatment columns are explicitly excluded"
            ),
            score_formula="mean over mapped genes of log1p(RPKM) for the native clinical-comparison subset",
            non_sample_columns=("Chr", "Start", "End", "Strand", "Length"),
        ),
        CohortMatrix(
            cohort_id="GSE167931_external_tpm_sensitivity",
            dataset="GSE167931",
            relative_path=(
                "data/raw/geo_candidates/GSE167931/GSE167931_AllSamplesTPMValue.txt.gz"
            ),
            gene_column="GeneID",
            gene_parser=ensembl_parenthetical_symbol,
            sample_to_gsm=GSE167931_COLUMNS,
            sample_to_group={
                column: group_from_normal_or_degenerated(column)
                for column in GSE167931_COLUMNS
            },
            source_value_scale="GEO-deposited TPM",
            analysis_role="data_processing_sensitivity_not_an_additional_cohort",
            matrix_status=(
                "GEO supplementary TPM matrix paired to the same nine GSE167931 samples; "
                "not eligible to enter a meta-analysis together with the FPKM representation"
            ),
            score_formula="mean over mapped genes of log1p(TPM)",
        ),
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_modules(path: Path) -> dict[str, list[str]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    modules = {
        str(item["module_id"]): [str(gene).strip().upper() for gene in item["genes"]]
        for item in payload["modules"]
    }
    if tuple(modules) != MODULE_IDS:
        raise ValueError("The locked module configuration has unexpected module IDs or ordering")
    for module_id, genes in modules.items():
        if not genes or len(genes) != len(set(genes)):
            raise ValueError(f"Invalid locked gene set: {module_id}")
    return modules


def read_ledger(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    expected = {"dataset", "gsm", "sample_title", "raw_characteristics"}
    if not rows or not expected.issubset(rows[0]):
        raise ValueError(f"Candidate ledger is malformed: {path}")
    ledger = {row["gsm"].strip(): row for row in rows}
    if len(ledger) != len(rows) or any(not gsm for gsm in ledger):
        raise ValueError(f"Candidate ledger has duplicate or empty GSM values: {path}")
    return ledger


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def read_expression_matrix(path: Path, gene_column: str) -> tuple[pd.DataFrame, int]:
    """Read a tabular GEO matrix, allowing only audited trailing-empty fields."""
    trailing_empty_rows = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if not header or gene_column not in header:
            raise ValueError(
                f"{path}: the declared gene identifier column {gene_column!r} is absent; "
                f"header beginning={header[:5]!r}"
            )
        expected_width = len(header)
        for line_number, raw in enumerate(handle, start=2):
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) == expected_width:
                continue
            if len(fields) == expected_width + 1 and fields[-1] == "":
                trailing_empty_rows += 1
                continue
            raise ValueError(
                f"{path}: row {line_number} has {len(fields)} fields; expected {expected_width}, "
                "or one explicitly empty trailing field"
            )
    frame = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        usecols=list(range(expected_width)),
        dtype={gene_column: str},
    )
    return frame, trailing_empty_rows


def log1p_cpm_scores(values: pd.DataFrame, mapped_symbols: list[str]) -> pd.DataFrame:
    totals = values.sum(axis=0)
    if (totals <= 0).any() or not np.isfinite(totals.to_numpy(dtype=float)).all():
        raise ValueError("A count-matrix sample has a non-positive or non-finite total")
    return np.log1p(values.loc[mapped_symbols].multiply(1_000_000.0 / totals, axis=1)).mean(axis=0)


def log1p_processed_scores(values: pd.DataFrame, mapped_symbols: list[str]) -> pd.DataFrame:
    return np.log1p(values.loc[mapped_symbols]).mean(axis=0)


def score_matrix(
    spec: CohortMatrix,
    project_root: Path,
    modules: dict[str, list[str]],
    ledgers: dict[str, dict[str, dict[str, str]]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    matrix_path = project_root / spec.relative_path
    if not matrix_path.is_file():
        raise FileNotFoundError(f"Required GEO matrix is missing: {matrix_path}")
    frame, trailing_empty_rows = read_expression_matrix(matrix_path, spec.gene_column)
    if spec.gene_column not in frame.columns:
        raise ValueError(f"{matrix_path}: missing gene identifier column {spec.gene_column!r}")
    expected_columns = list(spec.sample_to_gsm)
    sample_columns = [column for column in expected_columns if column in frame.columns]
    if set(sample_columns) != set(expected_columns) or len(sample_columns) != len(expected_columns):
        raise ValueError(
            f"{matrix_path}: matrix sample columns do not exactly match its audited mapping; "
            f"observed={sample_columns}, expected={expected_columns}"
        )
    if frame.empty:
        raise ValueError(f"{matrix_path}: no gene rows")
    values = frame[sample_columns].apply(pd.to_numeric, errors="coerce")
    if not np.isfinite(values.to_numpy(dtype=float)).all() or (values.to_numpy(dtype=float) < 0).any():
        raise ValueError(f"{matrix_path}: non-finite or negative expression value")

    symbols = frame[spec.gene_column].map(spec.gene_parser)
    if symbols.eq("").all():
        raise ValueError(f"{matrix_path}: no parsable gene symbols")
    feature_ids = frame[spec.gene_column].astype(str).str.extract(r"^(ENSG\d+)", expand=False).fillna("").str.upper()
    valid_symbols = symbols.ne("")
    values = values.loc[valid_symbols].copy()
    symbols = symbols.loc[valid_symbols].reset_index(drop=True)
    feature_ids = feature_ids.loc[valid_symbols].reset_index(drop=True)
    values.index = symbols.to_numpy()
    duplicate_symbol_counts = values.index.value_counts()
    duplicate_locked = sorted(
        symbol
        for symbol, count in duplicate_symbol_counts.items()
        if count > 1 and any(symbol in genes for genes in modules.values())
    )
    duplicate_resolution: list[str] = []
    keep_rows = np.ones(len(values), dtype=bool)
    for symbol in duplicate_locked:
        canonical_feature_id = CANONICAL_ENSEMBL_FOR_DUPLICATE_SYMBOL.get(symbol)
        positions = np.flatnonzero(values.index.to_numpy() == symbol)
        canonical_positions = [
            position for position in positions if feature_ids.iloc[position] == canonical_feature_id
        ]
        if canonical_feature_id is None or len(canonical_positions) != 1:
            raise ValueError(
                f"{matrix_path}: duplicated locked-gene symbol {symbol!r} has no uniquely declared "
                "canonical Ensembl resolution"
            )
        keep_rows[positions] = False
        keep_rows[canonical_positions[0]] = True
        alternate_ids = [feature_ids.iloc[position] for position in positions if position != canonical_positions[0]]
        duplicate_resolution.append(
            f"{symbol}:{canonical_feature_id}_retained;{','.join(alternate_ids)}_excluded"
        )
    values = values.iloc[keep_rows]
    values = values.loc[~values.index.duplicated(keep="first")]

    ledger = ledgers[spec.dataset]
    identity_rows: list[dict[str, object]] = []
    for column in expected_columns:
        gsm = spec.sample_to_gsm[column]
        ledger_row = ledger.get(gsm)
        if ledger_row is None:
            raise ValueError(f"{spec.dataset}: matrix column {column} maps to an absent ledger GSM {gsm}")
        if ledger_row["dataset"].strip() != spec.dataset:
            raise ValueError(f"{spec.dataset}: ledger dataset mismatch for {gsm}")
        identity_rows.append(
            {
                "cohort_id": spec.cohort_id,
                "dataset": spec.dataset,
                "gsm": gsm,
                "matrix_column": column,
                "sample_title": ledger_row["sample_title"],
                "disease_state": spec.sample_to_group[column],
                "contrast_arm": (
                    "comparison"
                    if spec.sample_to_group[column] in {"early_control", "normal"}
                    else "target"
                ),
                "identity_status": "pass_exact_predeclared_matrix_column_to_GSM",
                "public_person_level_metadata": (
                    "not reported in GEO sample characteristics; GSM remains a presumed donor/sample key"
                ),
            }
        )

    mapping_rows: list[dict[str, object]] = []
    score_rows: list[dict[str, object]] = []
    for module_id, genes in modules.items():
        mapped = [gene for gene in genes if gene in values.index]
        mapped_fraction = len(mapped) / len(genes)
        mapping_pass = mapped_fraction >= MIN_MAPPED_FRACTION
        mapping_rows.append(
            {
                "cohort_id": spec.cohort_id,
                "dataset": spec.dataset,
                "matrix_file": spec.relative_path,
                "source_value_scale": spec.source_value_scale,
                "module_id": module_id,
                "configured_genes": len(genes),
                "mapped_genes": len(mapped),
                "mapped_fraction": f"{mapped_fraction:.6f}",
                "mapping_pass": str(mapping_pass).lower(),
                "mapped_gene_symbols": ";".join(mapped),
                "duplicate_locked_gene_symbols": ";".join(duplicate_resolution),
                "matrix_unique_gene_symbols": len(values.index),
                "matrix_sample_columns": len(sample_columns),
            }
        )
        if not mapping_pass:
            raise ValueError(
                f"{spec.cohort_id}/{module_id}: locked-gene mapping is below {MIN_MAPPED_FRACTION:.2f}"
            )
        if spec.source_value_scale.endswith("counts"):
            scores = log1p_cpm_scores(values, mapped)
            totals = values.sum(axis=0)
        else:
            scores = log1p_processed_scores(values, mapped)
            totals = pd.Series(np.nan, index=sample_columns, dtype=float)
        for column in expected_columns:
            score_rows.append(
                {
                    "dataset": spec.dataset,
                    "gsm": spec.sample_to_gsm[column],
                    "donor_id": spec.sample_to_gsm[column],
                    "compartment": "NP",
                    "disease_state": spec.sample_to_group[column],
                    "module_id": module_id,
                    "module_score_log1p_cpm": f"{float(scores[column]):.8f}",
                    "score_status": "score_available",
                    "mapped_fraction": f"{mapped_fraction:.6f}",
                    "included_cells": "",
                    "total_umi_included_cells": "" if math.isnan(float(totals[column])) else f"{float(totals[column]):.8f}",
                    "matrix_status": spec.matrix_status,
                    "analysis_role": spec.analysis_role,
                    "source_value_scale": spec.source_value_scale,
                    "score_formula": spec.score_formula,
                    "matrix_column": column,
                }
            )

    matrix_audit = {
        "cohort_id": spec.cohort_id,
        "dataset": spec.dataset,
        "matrix_file": spec.relative_path,
        "matrix_sha256": sha256(matrix_path),
        "matrix_bytes": matrix_path.stat().st_size,
        "source_value_scale": spec.source_value_scale,
        "gene_identifier_column": spec.gene_column,
        "gene_rows": len(frame),
        "unique_nonempty_gene_symbols": len(values.index),
        "matrix_sample_columns": len(sample_columns),
        "expected_sample_columns": len(expected_columns),
        "matrix_column_mapping_status": "pass_exact_set_match",
        "numeric_values_finite_and_nonnegative": "true",
        "rows_with_explicit_empty_trailing_field_dropped": trailing_empty_rows,
        "duplicate_locked_gene_symbols": ";".join(duplicate_resolution),
        "all_locked_modules_mapping_pass": "true",
    }
    return score_rows, mapping_rows, identity_rows, matrix_audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--module-config", type=Path, default=Path("config/program_modules.json"))
    parser.add_argument("--candidate-ledger-dir", type=Path, default=Path("data/derived/geo_candidate_audit"))
    parser.add_argument("--output-root", type=Path, default=Path("data/derived/module_scores_external"))
    parser.add_argument("--audit-output-dir", type=Path, default=Path("data/derived/geo_candidate_audit"))
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    module_config = args.module_config if args.module_config.is_absolute() else project_root / args.module_config
    candidate_ledger_dir = (
        args.candidate_ledger_dir if args.candidate_ledger_dir.is_absolute() else project_root / args.candidate_ledger_dir
    )
    output_root = args.output_root if args.output_root.is_absolute() else project_root / args.output_root
    audit_output_dir = (
        args.audit_output_dir if args.audit_output_dir.is_absolute() else project_root / args.audit_output_dir
    )
    modules = load_modules(module_config)
    ledgers = {
        dataset: read_ledger(candidate_ledger_dir / f"{dataset}_sample_ledger.csv")
        for dataset in ("GSE186542", "GSE167931", "GSE245147")
    }

    all_matrix_audits: list[dict[str, object]] = []
    all_identity_rows: list[dict[str, object]] = []
    all_mapping_rows: list[dict[str, object]] = []
    generated_paths: list[Path] = []
    for spec in matrix_specs():
        scores, mappings, identities, matrix_audit = score_matrix(spec, project_root, modules, ledgers)
        matrix_stem = Path(spec.relative_path).name.removesuffix(".txt.gz")
        cohort_dir = output_root / spec.dataset
        score_path = cohort_dir / f"{matrix_stem}_module_scores.csv"
        mapping_path = cohort_dir / f"{matrix_stem}_module_mapping_audit.csv"
        write_csv(score_path, scores)
        write_csv(mapping_path, mappings)
        generated_paths.extend((score_path, mapping_path))
        all_matrix_audits.append(matrix_audit)
        all_identity_rows.extend(identities)
        all_mapping_rows.extend(mappings)

    audit_output_dir.mkdir(parents=True, exist_ok=True)
    matrix_audit_path = audit_output_dir / "external_human_np_candidate_matrix_audit.csv"
    identity_audit_path = audit_output_dir / "external_human_np_candidate_identity_audit.csv"
    mapping_audit_path = audit_output_dir / "external_human_np_candidate_module_mapping_audit.csv"
    write_csv(matrix_audit_path, all_matrix_audits)
    write_csv(identity_audit_path, all_identity_rows)
    write_csv(mapping_audit_path, all_mapping_rows)
    generated_paths.extend((matrix_audit_path, identity_audit_path, mapping_audit_path))

    manifest = {
        "schema_version": 1,
        "purpose": "Audited scoring of GSE186542 and GSE167931 external human NP cohorts with project-locked modules.",
        "inference_boundary": (
            "Each GSM is a presumed donor/sample key. All outputs are exploratory score-level support only; "
            "no causal, age-adjusted, confirmatory, or cell-level inference."
        ),
        "included_meta_analysis_representations": [
        "GSE186542_external_count_support",
        "GSE167931_external_fpkm_support",
        "GSE245147_external_native_comparison_support",
        ],
        "excluded_as_duplicate_representation": ["GSE167931_external_tpm_sensitivity"],
        "module_config_sha256": sha256(module_config),
        "input_sha256": {
            spec.relative_path: sha256(project_root / spec.relative_path) for spec in matrix_specs()
        },
        "generated_artifact_sha256": {str(path.relative_to(project_root)).replace("\\", "/"): sha256(path) for path in generated_paths},
    }
    manifest_path = audit_output_dir / "external_human_np_candidate_scoring_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Audited external NP cohort scoring completed: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
