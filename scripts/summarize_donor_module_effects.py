"""Summarize donor/library-level IVDD module-score contrasts.

This script deliberately operates only on already aggregated module-score
tables.  It does not read cell-level expression values and it does not emit
cell-level tests or p-values.  For every cohort x compartment x module
contrast, the estimand is the unweighted difference in mean donor/library
module score (target severity group minus comparison group).

It writes donor/library group descriptives, Welch analytic confidence
intervals where their prerequisites are met, donor-bootstrap percentile
confidence intervals, leave-one-key-out effects, and a descriptive
cross-cohort direction-alignment table.  The latter is not a meta-analysis or
a replication adjudication.

The default current-project specification includes only score files already
available in this repository.  GSE230809 is retained as one exploratory
parent project, not two independent cohorts.  GSE244889 is directional only,
and GSE165722 is a normalized-count score-level direction/stability cohort.
Future cohorts can be supplied through --contrast-spec after their own score
tables and provenance gates have passed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import t as student_t


SPEC_COLUMNS = [
    "cohort_id",
    "cohort_role",
    "inference_key_type",
    "score_files",
    "dataset_values",
    "compartment",
    "group_column",
    "comparison_value",
    "target_value",
    "comparison_label",
    "target_label",
    "contrast_label",
    "severity_orientation",
    "analysis_role",
    "confirmatory_eligible",
    "notes",
]

REQUIRED_SCORE_COLUMNS = {
    "dataset",
    "gsm",
    "donor_id",
    "compartment",
    "module_id",
    "module_score_log1p_cpm",
}

SUMMARY_CONTRACTS = {
    "default_current_project": {
        "manifest_value": "default_current_project_cross_cohort_summary_only",
        "title": "Donor/Library Module-Effect Summary",
        "authoritative_default": True,
        "boundary": (
            "This is the sole authoritative frozen default descriptive summary. It contains 20 effects "
            "from the four locked cohort roles and remains separate from S7, S8, and S9."
        ),
    },
    "post_hoc_external_expansion": {
        "manifest_value": "post_hoc_external_expansion_intermediate_summary_only",
        "title": "Post Hoc External-Expansion Donor/Library Module-Effect Summary",
        "authoritative_default": False,
        "boundary": (
            "This non-authoritative 24-effect staging summary supplies study-level inputs to S8. "
            "It does not replace or revise data/derived/donor_module_effect_summary/."
        ),
    },
    "source_family_replacement": {
        "manifest_value": "source_family_replacement_intermediate_summary_only",
        "title": "Source-Family Replacement Donor/Library Module-Effect Summary",
        "authoritative_default": False,
        "boundary": (
            "This non-authoritative 24-effect staging summary supplies the GSE245147-for-GSE167931 "
            "replacement inputs to S9. It does not replace or revise the frozen default summary."
        ),
    },
    "custom_sensitivity": {
        "manifest_value": "custom_non_authoritative_descriptive_summary",
        "title": "Non-Authoritative Donor/Library Module-Effect Sensitivity Summary",
        "authoritative_default": False,
        "boundary": (
            "This is a custom descriptive sensitivity/audit summary. It must not be substituted for "
            "data/derived/donor_module_effect_summary/."
        ),
    },
}


def current_project_spec() -> list[dict[str, str]]:
    """Return the locked, auditable contrast map for currently scored cohorts."""
    discovery_files = ";".join(
        [
            "data/derived/module_scores_recomputed/GSE229711/GSE229711_RAW_module_scores.csv",
            "data/derived/module_scores_recomputed/GSE230808/GSE230808_RAW_module_scores.csv",
        ]
    )
    discovery_note = (
        "One parent project; age and disease status are fully confounded, and the "
        "healthy group has n=3 per compartment. Exploratory effect-size and "
        "stability display only."
    )
    rows = []
    for compartment in ("AF", "NP"):
        rows.append(
            {
                "cohort_id": "GSE230809_discovery",
                "cohort_role": "exploratory_discovery_parent_project",
                "inference_key_type": "donor_id",
                "score_files": discovery_files,
                "dataset_values": "GSE229711;GSE230808",
                "compartment": compartment,
                "group_column": "disease_state",
                "comparison_value": "healthy",
                "target_value": "diseased",
                "comparison_label": "healthy_low_grade",
                "target_label": "advanced_degeneration_associated",
                "contrast_label": "advanced_degeneration_associated_minus_healthy_low_grade",
                "severity_orientation": "higher_severity_minus_lower_severity",
                "analysis_role": "exploratory_effect_size_and_stability_only",
                "confirmatory_eligible": "false",
                "notes": discovery_note,
            }
        )
    rows.append(
        {
            "cohort_id": "GSE244889_directional",
            "cohort_role": "external_directional_support_only",
            "inference_key_type": "presumed_donor_or_library_key",
            "score_files": (
                "data/derived/module_scores_recomputed/GSE244889/"
                "GSE244889_RAW_module_scores.csv"
            ),
            "dataset_values": "GSE244889",
            "compartment": "NP",
            "group_column": "disease_state",
            "comparison_value": "MDD",
            "target_value": "SDD",
            "comparison_label": "mild_MDD",
            "target_label": "severe_SDD",
            "contrast_label": "severe_SDD_minus_mild_MDD",
            "severity_orientation": "higher_severity_minus_lower_severity",
            "analysis_role": "directional_effect_size_and_stability_only",
            "confirmatory_eligible": "false",
            "notes": (
                "Four MDD versus three SDD title-derived presumed donor/library keys; "
                "directional support only, not confirmatory validation."
            ),
        }
    )
    rows.append(
        {
            "cohort_id": "GSE153066_support",
            "cohort_role": "external_count_level_support_only",
            "inference_key_type": "presumed_donor_or_library_key",
            "score_files": (
                "data/derived/module_scores_external/GSE153066/"
                "GSE153066_AllSample.counts_module_scores.csv"
            ),
            "dataset_values": "GSE153066",
            "compartment": "NP",
            "group_column": "disease_state",
            "comparison_value": "relatively normal",
            "target_value": "degenerated",
            "comparison_label": "relatively_normal",
            "target_label": "degenerated",
            "contrast_label": "degenerated_minus_relatively_normal",
            "severity_orientation": "higher_severity_minus_lower_severity",
            "analysis_role": "external_support_effect_size_and_stability_only",
            "confirmatory_eligible": "false",
            "notes": (
                "Eight relatively-normal versus eight degenerated presumed sample/library keys. "
                "The GEO-retained dense count matrix follows contributor cell filtering; clinical "
                "source and age are confounded with disease status. External support only, not an "
                "age-independent or causal contrast."
            ),
        }
    )
    rows.append(
        {
            "cohort_id": "GSE165722_score_level",
            "cohort_role": "external_normalized_count_score_level_direction_only",
            "inference_key_type": "presumed_donor_level_sample_key",
            "score_files": (
                "data/derived/module_scores_external/GSE165722/"
                "GSE165722_RAW_module_scores.csv"
            ),
            "dataset_values": "GSE165722",
            "compartment": "NP",
            "group_column": "disease_state",
            "comparison_value": "mild",
            "target_value": "severe",
            "comparison_label": "mild",
            "target_label": "severe",
            "contrast_label": "severe_minus_mild",
            "severity_orientation": "higher_severity_minus_lower_severity",
            "analysis_role": "score_level_direction_and_stability_only",
            "confirmatory_eligible": "false",
            "notes": (
                "Four mild versus four severe presumed sample keys, with severity assigned from the "
                "source-publication grouping because GEO and source-publication grades differ. GEO "
                "describes supplied values as normalized counts; this cohort supports only score-level "
                "direction and leave-one-key-out stability, not raw-count, cell-level, confirmatory, "
                "causal, or age-adjusted inference."
            ),
        }
    )
    return rows


def current_project_score_availability() -> list[dict[str, str]]:
    """Describe all audited human cohorts without treating absent scores as null effects."""
    return [
        {
            "cohort_id": "GSE230809_discovery",
            "dataset_or_parent_project": "GSE229711 + GSE230808 (GSE230809 parent project)",
            "compartment_scope": "AF;NP",
            "sample_key_type": "donor_id",
            "audited_group_structure": "healthy n=3 per compartment; diseased AF n=10, NP n=8",
            "module_score_availability": "available",
            "score_table_paths": (
                "data/derived/module_scores_recomputed/GSE229711/"
                "GSE229711_RAW_module_scores.csv;"
                "data/derived/module_scores_recomputed/GSE230808/"
                "GSE230808_RAW_module_scores.csv"
            ),
            "effect_summary_included": "true",
            "analysis_boundary": "one exploratory parent project; age/disease confounded; no confirmatory inference",
        },
        {
            "cohort_id": "GSE244889_directional",
            "dataset_or_parent_project": "GSE244889 scRNA-seq",
            "compartment_scope": "NP",
            "sample_key_type": "presumed donor/library key",
            "audited_group_structure": "MDD n=4; SDD n=3",
            "module_score_availability": "available",
            "score_table_paths": (
                "data/derived/module_scores_recomputed/GSE244889/"
                "GSE244889_RAW_module_scores.csv"
            ),
            "effect_summary_included": "true",
            "analysis_boundary": "directional effect/stability display only; no confirmatory validation",
        },
        {
            "cohort_id": "GSE153066_support",
            "dataset_or_parent_project": "GSE153066",
            "compartment_scope": "NP",
            "sample_key_type": "presumed donor/library key",
            "audited_group_structure": "relatively normal n=8; degenerated n=8",
            "module_score_availability": "available",
            "score_table_paths": (
                "data/derived/module_scores_external/GSE153066/"
                "GSE153066_AllSample.counts_module_scores.csv"
            ),
            "effect_summary_included": "true",
            "analysis_boundary": (
                "external count-level support with contributor-retained cells; clinical source and age are confounded; "
                "no age-independent or causal interpretation"
            ),
        },
        {
            "cohort_id": "GSE165722_score_level",
            "dataset_or_parent_project": "GSE165722",
            "compartment_scope": "NP",
            "sample_key_type": "presumed donor-level sample key",
            "audited_group_structure": "mild n=4; severe n=4",
            "module_score_availability": "available",
            "score_table_paths": (
                "data/derived/module_scores_external/GSE165722/"
                "GSE165722_RAW_module_scores.csv"
            ),
            "effect_summary_included": "true",
            "analysis_boundary": (
                "normalized-count score-level direction/stability only; source-publication severity grouping used "
                "because GEO and source-publication grades differ; no raw-count, cell-level, confirmatory, "
                "causal, or age-adjusted inference"
            ),
        },
        {
            "cohort_id": "GSE251686_exploratory",
            "dataset_or_parent_project": "GSE251686",
            "compartment_scope": "NP",
            "sample_key_type": "presumed sample/library key",
            "audited_group_structure": "mild n=2; severe n=3 after exclusion of GSM7986002",
            "module_score_availability": "available_separate_exploratory_output",
            "score_table_paths": (
                "data/derived/GSE251686_exploratory_scores/"
                "GSE251686_exploratory_module_scores.csv"
            ),
            "effect_summary_included": "false",
            "analysis_boundary": (
                "separate isolated exploratory score/effect display only; deliberately excluded from the "
                "default 20-effect summary; GSM7986002 remains permanently excluded"
            ),
        },
    ]


def parse_bool(value: Any, field: str) -> bool:
    normalized = str(value).strip().lower()
    if normalized not in {"true", "false"}:
        raise ValueError(f"{field} must be true or false, got {value!r}")
    return normalized == "true"


def parse_delimited(value: Any) -> list[str]:
    return [item.strip() for item in str(value).split(";") if item.strip()]


def resolve_path(value: str, project_root: Path, spec_directory: Path | None) -> Path:
    """Resolve an input path without silently accepting an ambiguous location."""
    given = Path(value)
    candidates = [given] if given.is_absolute() else [project_root / given]
    if spec_directory is not None and not given.is_absolute():
        candidates.append(spec_directory / given)
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    searched = ", ".join(str(candidate) for candidate in candidates)
    raise FileNotFoundError(f"Score table not found. Searched: {searched}")


def load_contrast_spec(path: Path) -> list[dict[str, str]]:
    # utf-8-sig also accepts ordinary UTF-8 while handling CSVs written by
    # Excel or PowerShell with a leading UTF-8 byte-order mark.
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = set(reader.fieldnames or [])
        missing = set(SPEC_COLUMNS) - fieldnames
        if missing:
            raise ValueError(f"Contrast specification missing columns: {sorted(missing)}")
        rows = [{column: str(row.get(column, "")).strip() for column in SPEC_COLUMNS} for row in reader]
    if not rows:
        raise ValueError(f"Contrast specification is empty: {path}")
    return rows


def validate_spec(rows: list[dict[str, str]]) -> None:
    seen: set[tuple[str, str]] = set()
    for row in rows:
        required = [
            "cohort_id",
            "cohort_role",
            "inference_key_type",
            "score_files",
            "dataset_values",
            "compartment",
            "group_column",
            "comparison_value",
            "target_value",
            "comparison_label",
            "target_label",
            "contrast_label",
            "severity_orientation",
            "analysis_role",
        ]
        missing = [field for field in required if not row[field]]
        if missing:
            raise ValueError(f"Incomplete contrast specification {row!r}: missing {missing}")
        if row["comparison_value"] == row["target_value"]:
            raise ValueError(f"{row['cohort_id']}/{row['compartment']}: comparison and target values match")
        parse_bool(row["confirmatory_eligible"], "confirmatory_eligible")
        key = (row["cohort_id"], row["compartment"])
        if key in seen:
            raise ValueError(f"Duplicate cohort/compartment contrast specification: {key}")
        seen.add(key)


def ci_percentiles(confidence_level: float) -> tuple[float, float]:
    tail = (1.0 - confidence_level) / 2.0
    return tail, 1.0 - tail


def numeric_or_blank(value: float | int | None, digits: int = 8) -> str:
    if value is None or not math.isfinite(float(value)):
        return ""
    return f"{float(value):.{digits}f}"


def direction(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "not_estimable"
    if value > 0:
        return "positive_higher_target_severity"
    if value < 0:
        return "negative_lower_target_severity"
    return "zero"


def safe_sample_sd(values: np.ndarray) -> float | None:
    if len(values) < 2:
        return None
    value = float(np.std(values, ddof=1))
    return value if math.isfinite(value) else None


def analytic_difference_ci(
    target: np.ndarray, comparison: np.ndarray, confidence_level: float
) -> tuple[float | None, float | None, float | None, float | None, str]:
    """Welch t confidence interval; return an explicit status when unavailable."""
    if len(target) < 2 or len(comparison) < 2:
        return None, None, None, None, "not_estimable_group_n_lt_2"
    target_var = float(np.var(target, ddof=1))
    comparison_var = float(np.var(comparison, ddof=1))
    if not (math.isfinite(target_var) and math.isfinite(comparison_var)):
        return None, None, None, None, "not_estimable_nonfinite_variance"
    se_sq = target_var / len(target) + comparison_var / len(comparison)
    if not math.isfinite(se_sq) or se_sq <= 0:
        return None, None, None, None, "not_estimable_zero_or_nonpositive_standard_error"
    denominator = (target_var**2) / (len(target) ** 2 * (len(target) - 1))
    denominator += (comparison_var**2) / (len(comparison) ** 2 * (len(comparison) - 1))
    if not math.isfinite(denominator) or denominator <= 0:
        return None, None, None, None, "not_estimable_welch_df"
    df = se_sq**2 / denominator
    if not math.isfinite(df) or df <= 0:
        return None, None, None, None, "not_estimable_welch_df"
    critical = float(student_t.ppf((1.0 + confidence_level) / 2.0, df))
    if not math.isfinite(critical):
        return None, None, None, None, "not_estimable_welch_critical_value"
    effect = float(np.mean(target) - np.mean(comparison))
    half_width = critical * math.sqrt(se_sq)
    return effect - half_width, effect + half_width, math.sqrt(se_sq), df, "available_welch_t"


def stable_seed(seed: int, cohort_id: str, compartment: str, module_id: str) -> int:
    payload = f"{seed}|{cohort_id}|{compartment}|{module_id}".encode("utf-8")
    digest = hashlib.sha256(payload).digest()
    return int.from_bytes(digest[:8], "little", signed=False)


def bootstrap_difference_ci(
    target: np.ndarray,
    comparison: np.ndarray,
    confidence_level: float,
    n_bootstrap: int,
    seed: int,
) -> tuple[float | None, float | None, str]:
    """Independent donor/library bootstrap, never a cell-level resample."""
    if len(target) < 2 or len(comparison) < 2:
        return None, None, "not_estimable_group_n_lt_2"
    if n_bootstrap < 1:
        return None, None, "not_run_n_bootstrap_lt_1"
    rng = np.random.default_rng(seed)
    target_samples = rng.choice(target, size=(n_bootstrap, len(target)), replace=True)
    comparison_samples = rng.choice(comparison, size=(n_bootstrap, len(comparison)), replace=True)
    differences = target_samples.mean(axis=1) - comparison_samples.mean(axis=1)
    if not np.isfinite(differences).all():
        return None, None, "not_estimable_nonfinite_bootstrap_draw"
    lower_q, upper_q = ci_percentiles(confidence_level)
    lower, upper = np.quantile(differences, [lower_q, upper_q], method="linear")
    return float(lower), float(upper), "available_percentile_donor_bootstrap"


def mean_ci(values: np.ndarray, confidence_level: float) -> tuple[float | None, float | None, str]:
    if len(values) < 2:
        return None, None, "not_estimable_group_n_lt_2"
    sd = safe_sample_sd(values)
    if sd is None:
        return None, None, "not_estimable_nonfinite_sd"
    se = sd / math.sqrt(len(values))
    if not math.isfinite(se) or se <= 0:
        return None, None, "not_estimable_zero_or_nonpositive_standard_error"
    critical = float(student_t.ppf((1.0 + confidence_level) / 2.0, len(values) - 1))
    estimate = float(np.mean(values))
    return estimate - critical * se, estimate + critical * se, "available_t"


def rows_to_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def optional_numeric_column(frame: pd.DataFrame, column: str) -> pd.Series:
    """Return an all-missing numeric series when an optional score-table field is absent."""
    if column not in frame.columns:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_run_manifest(
    output_dir: Path,
    artifact_names: list[str],
    specs: list[dict[str, str]],
    parameter_rows: list[dict[str, Any]],
    effect_rows: list[dict[str, Any]],
    summary_contract: str,
) -> None:
    """Make a result directory self-identifying and resistant to stale snapshots."""
    manifest_path = output_dir / "run_manifest.json"
    artifacts = {
        name: file_sha256(output_dir / name)
        for name in artifact_names
        if name != manifest_path.name and name != "run_artifacts.csv" and (output_dir / name).is_file()
    }
    contract = SUMMARY_CONTRACTS[summary_contract]
    cohort_ids = sorted({str(row["cohort_id"]) for row in effect_rows})
    manifest = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "result_contract": contract["manifest_value"],
        "authoritative_default_summary": contract["authoritative_default"],
        "result_directory": str(output_dir),
        "summary_effect_count": len(effect_rows),
        "summary_cohort_ids": cohort_ids,
        "default_authoritative_result_directory": "data/derived/donor_module_effect_summary/",
        "all_effects_confirmatory_eligible": all(
            str(row["confirmatory_eligible"]).strip().lower() == "true" for row in effect_rows
        ),
        "excluded_separate_exploratory_package": {
            "cohort_id": "GSE251686_exploratory",
            "score_package": "data/derived/GSE251686_exploratory_scores/",
            "reason": (
                "A separate mild n=2 versus severe n=3 exploratory score/effect package exists, but it is "
                "intentionally excluded from this summary and is not confirmatory evidence."
            ),
        },
        "no_p_values_or_multiple_testing_adjustment": True,
        "contrast_specification": specs,
        "input_parameter_rows": parameter_rows,
        "generated_artifact_sha256": artifacts,
    }
    if contract["authoritative_default"]:
        manifest["authoritative_result_directory"] = str(output_dir)
        manifest["default_summary_effect_count"] = len(effect_rows)
        manifest["default_summary_cohort_ids"] = cohort_ids
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run_parameter_rows(
    specs: list[dict[str, str]],
    project_root: Path,
    spec_directory: Path | None,
    confidence_level: float,
    n_bootstrap: int,
    random_seed: int,
    min_mapped_fraction: float,
    summary_contract: str,
) -> list[dict[str, Any]]:
    """Record the exact score-table bytes and non-defaultable summary settings."""
    source_files: list[Path] = []
    for spec in specs:
        for value in parse_delimited(spec["score_files"]):
            path = resolve_path(value, project_root, spec_directory)
            if path not in source_files:
                source_files.append(path)
    rows: list[dict[str, Any]] = []
    for path in source_files:
        rows.append(
            {
                "parameter_class": "score_table_input",
                "parameter_name": "score_table",
                "parameter_value": str(path),
                "sha256": file_sha256(path),
                "notes": "Exact input table used for donor/library-level summarization.",
            }
        )
    rows.extend(
        [
            {
                "parameter_class": "analysis_script",
                "parameter_name": "summarize_donor_module_effects.py",
                "parameter_value": str(project_root / "scripts/summarize_donor_module_effects.py"),
                "sha256": file_sha256(project_root / "scripts/summarize_donor_module_effects.py"),
                "notes": "Exact implementation used to create this summary directory.",
            },
            {
                "parameter_class": "analysis_setting",
                "parameter_name": "summary_contract",
                "parameter_value": summary_contract,
                "sha256": "",
                "notes": SUMMARY_CONTRACTS[summary_contract]["boundary"],
            },
            {
                "parameter_class": "analysis_setting",
                "parameter_name": "confidence_level",
                "parameter_value": f"{confidence_level:.8f}",
                "sha256": "",
                "notes": "Used for Welch and percentile bootstrap intervals.",
            },
            {
                "parameter_class": "analysis_setting",
                "parameter_name": "bootstrap_replicates",
                "parameter_value": str(n_bootstrap),
                "sha256": "",
                "notes": "Independent resampling within donor/library contrast arms.",
            },
            {
                "parameter_class": "analysis_setting",
                "parameter_name": "random_seed",
                "parameter_value": str(random_seed),
                "sha256": "",
                "notes": "Cohort/compartment/module-specific deterministic seeds are derived from this root seed.",
            },
            {
                "parameter_class": "analysis_setting",
                "parameter_name": "min_mapped_fraction",
                "parameter_value": f"{min_mapped_fraction:.8f}",
                "sha256": "",
                "notes": "Score rows below this pre-specified module mapping threshold are excluded with an audit row.",
            },
        ]
    )
    return rows


def strict_identity_keys(
    path: Path,
    required_columns: list[str],
    ledger_group_column: str = "disease_state",
) -> pd.DataFrame:
    """Load and normalize a sample ledger identity crosswalk without silent remapping."""
    frame = pd.read_csv(path, dtype=str)
    ledger_required_columns = [
        ledger_group_column if column == "disease_state" else column for column in required_columns
    ]
    missing = set(ledger_required_columns) - set(frame.columns)
    if missing:
        raise ValueError(f"{path}: ledger missing identity columns {sorted(missing)}")
    keys = frame[ledger_required_columns].fillna("").copy()
    if ledger_group_column != "disease_state":
        keys = keys.rename(columns={ledger_group_column: "disease_state"})
    keys = keys[required_columns]
    empty_columns = [column for column in required_columns if keys[column].str.strip().eq("").any()]
    if empty_columns:
        raise ValueError(f"{path}: ledger has empty identity values in {empty_columns}")
    duplicates = keys.duplicated(keep=False)
    if duplicates.any():
        raise ValueError(
            f"{path}: ledger has duplicate identity keys:\n" + keys.loc[duplicates].to_csv(index=False)
        )
    return keys


def canonical_score_table_identity_keys(
    specs: list[dict[str, str]], project_root: Path, spec_directory: Path | None
) -> pd.DataFrame:
    """Obtain input-table sample identities before module-specific score filtering."""
    identity_columns = ["dataset", "gsm", "donor_id", "compartment", "disease_state"]
    selected: list[pd.DataFrame] = []
    for spec in specs:
        dataset_values = set(parse_delimited(spec["dataset_values"]))
        for file_value in parse_delimited(spec["score_files"]):
            path = resolve_path(file_value, project_root, spec_directory)
            frame = pd.read_csv(path, dtype=str)
            missing = set(identity_columns) - set(frame.columns)
            if missing:
                raise ValueError(f"{path}: score table missing identity columns {sorted(missing)}")
            if spec["group_column"] not in frame.columns:
                raise ValueError(
                    f"{path}: score table lacks contrast group column {spec['group_column']!r}"
                )
            selected_mask = (
                frame["dataset"].isin(dataset_values)
                & frame["compartment"].eq(spec["compartment"])
                & frame[spec["group_column"]].isin([spec["comparison_value"], spec["target_value"]])
            )
            candidate = frame.loc[selected_mask, identity_columns].fillna("").copy()
            candidate["cohort_id"] = spec["cohort_id"]
            selected.append(candidate)
    if not selected:
        raise ValueError("No score-table identities available for current-project ledger crosswalk")
    result = pd.concat(selected, ignore_index=True)
    for column in identity_columns:
        if result[column].str.strip().eq("").any():
            raise ValueError(f"Score table contains empty identity value in {column}")
    # Rows repeat per module.  Only conflicting metadata for a given GSM is an error.
    result = result.drop_duplicates()
    inconsistent = result.duplicated(["cohort_id", "dataset", "gsm"], keep=False)
    if inconsistent.any():
        raise ValueError(
            "Score table has conflicting sample identity metadata for the same cohort/dataset/GSM:\n"
            + result.loc[inconsistent].sort_values(["cohort_id", "dataset", "gsm"]).to_csv(index=False)
        )
    return result


def current_project_ledger_crosswalk(
    specs: list[dict[str, str]],
    retained_scores: pd.DataFrame,
    project_root: Path,
    spec_directory: Path | None,
) -> list[dict[str, Any]]:
    """Require exact input-score-table-to-ledger identity agreement for fixed current inputs."""
    identity_columns = ["dataset", "gsm", "donor_id", "compartment", "disease_state"]
    input_identities = canonical_score_table_identity_keys(specs, project_root, spec_directory)
    checks = [
        (
            "GSE230809_discovery",
            project_root / "data/derived/GSE230809_discovery_raw_data_ledger.csv",
            "GSE230809 discovery raw-data ledger",
            "disease_state",
        ),
        (
            "GSE244889_directional",
            project_root / "data/derived/GSE244889_scrna_sample_ledger.csv",
            "GSE244889 scRNA sample ledger",
            "disease_state",
        ),
        (
            "GSE153066_support",
            project_root / "data/derived/GSE153066_donor_ledger.csv",
            "GSE153066 presumed donor/library ledger",
            "disease_state",
        ),
        (
            "GSE165722_score_level",
            project_root / "data/derived/GSE165722_donor_ledger.csv",
            "GSE165722 donor ledger using source-publication severity grouping",
            "source_publication_severity_group",
        ),
    ]
    rows: list[dict[str, Any]] = []
    for cohort_id, ledger_path, ledger_label, ledger_group_column in checks:
        if not ledger_path.is_file():
            raise FileNotFoundError(f"Required current-project ledger not found: {ledger_path}")
        ledger_keys = strict_identity_keys(ledger_path, identity_columns, ledger_group_column)
        score_keys = input_identities.loc[
            input_identities["cohort_id"].eq(cohort_id), identity_columns
        ].copy()
        if score_keys.empty:
            raise ValueError(f"{cohort_id}: no score-table identities for required score-to-ledger crosswalk")
        retained_keys = retained_scores.loc[
            retained_scores["cohort_id"].eq(cohort_id), identity_columns
        ].fillna("").drop_duplicates()
        merged = score_keys.merge(ledger_keys, on=identity_columns, how="outer", indicator=True)
        score_only = int(merged["_merge"].eq("left_only").sum())
        ledger_only = int(merged["_merge"].eq("right_only").sum())
        matched = int(merged["_merge"].eq("both").sum())
        status = "pass_exact_identity_crosswalk" if not score_only and not ledger_only else "fail_identity_crosswalk"
        rows.append(
            {
                "cohort_id": cohort_id,
                "ledger_label": ledger_label,
                "ledger_path": str(ledger_path),
                "identity_columns": ";".join(identity_columns),
                "score_group_column": "disease_state",
                "ledger_group_column": ledger_group_column,
                "unique_score_table_sample_keys": len(score_keys),
                "unique_retained_score_sample_keys": len(retained_keys),
                "unique_ledger_sample_keys": len(ledger_keys),
                "matched_keys": matched,
                "score_table_only_keys": score_only,
                "ledger_only_keys": ledger_only,
                "crosswalk_status": status,
            }
        )
        if status != "pass_exact_identity_crosswalk":
            raise ValueError(
                f"{cohort_id}: score-table/ledger identity mismatch: "
                f"score_table_only={score_only}, ledger_only={ledger_only}"
            )
    return rows


def collect_scores(
    specs: list[dict[str, str]], project_root: Path, spec_directory: Path | None, min_mapped_fraction: float
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    """Read score tables, retain valid donor/library score records, and audit all exclusions."""
    selected: list[pd.DataFrame] = []
    audit_rows: list[dict[str, Any]] = []
    for spec_index, spec in enumerate(specs):
        file_values = parse_delimited(spec["score_files"])
        dataset_values = set(parse_delimited(spec["dataset_values"]))
        frames: list[pd.DataFrame] = []
        for file_value in file_values:
            path = resolve_path(file_value, project_root, spec_directory)
            frame = pd.read_csv(path, dtype=str)
            missing = REQUIRED_SCORE_COLUMNS - set(frame.columns)
            if missing:
                raise ValueError(f"{path}: score table missing columns {sorted(missing)}")
            frame["_score_file"] = str(path)
            frames.append(frame)
        raw = pd.concat(frames, ignore_index=True)
        matched_dataset = raw["dataset"].isin(dataset_values)
        matched_compartment = raw["compartment"].eq(spec["compartment"])
        selected_group = raw[spec["group_column"]].isin(
            [spec["comparison_value"], spec["target_value"]]
        ) if spec["group_column"] in raw.columns else None
        if selected_group is None:
            raise ValueError(
                f"{spec['cohort_id']}/{spec['compartment']}: group column "
                f"{spec['group_column']!r} is absent from supplied score tables"
            )
        candidate_mask = matched_dataset & matched_compartment & selected_group
        candidate = raw.loc[candidate_mask].copy()
        if candidate.empty:
            raise ValueError(
                f"{spec['cohort_id']}/{spec['compartment']}: no score rows match the selected contrast"
            )
        status = candidate.get("score_status", pd.Series("score_available", index=candidate.index)).fillna("")
        mapped_fraction = pd.to_numeric(
            candidate.get("mapped_fraction", pd.Series("1", index=candidate.index)), errors="coerce"
        )
        score = pd.to_numeric(candidate["module_score_log1p_cpm"], errors="coerce")
        valid_status = status.eq("score_available")
        valid_mapping = mapped_fraction.ge(min_mapped_fraction)
        valid_score = np.isfinite(score)
        invalid_reasons = np.where(
            ~valid_status,
            "score_status_not_available",
            np.where(~valid_mapping, "mapped_fraction_below_threshold", np.where(~valid_score, "nonfinite_score", "")),
        )
        for reason, count in pd.Series(invalid_reasons).value_counts(dropna=False).items():
            audit_rows.append(
                {
                    "cohort_id": spec["cohort_id"],
                    "compartment": spec["compartment"],
                    "spec_index": spec_index,
                    "source_score_files": ";".join(str(path) for path in sorted({Path(value).name for value in file_values})),
                    "raw_rows_loaded": len(raw),
                    "rows_matching_dataset": int(matched_dataset.sum()),
                    "rows_matching_compartment": int((matched_dataset & matched_compartment).sum()),
                    "rows_matching_contrast_groups": len(candidate),
                    "row_disposition": "retained" if reason == "" else "excluded",
                    "reason": "score_available_and_mapping_pass" if reason == "" else str(reason),
                    "row_count": int(count),
                }
            )
        candidate = candidate.loc[valid_status & valid_mapping & valid_score].copy()
        candidate["module_score"] = score.loc[candidate.index].astype(float)
        candidate["mapped_fraction_numeric"] = mapped_fraction.loc[candidate.index].astype(float)
        candidate["severity_group"] = np.where(
            candidate[spec["group_column"]].eq(spec["comparison_value"]),
            spec["comparison_label"],
            spec["target_label"],
        )
        candidate["contrast_arm"] = np.where(
            candidate[spec["group_column"]].eq(spec["comparison_value"]), "comparison", "target"
        )
        for column, value in spec.items():
            candidate[column] = value
        candidate["spec_index"] = spec_index
        selected.append(candidate)
    combined = pd.concat(selected, ignore_index=True)
    observation_columns = ["cohort_id", "compartment", "donor_id"]
    if combined["donor_id"].fillna("").str.strip().eq("").any():
        raise ValueError("A retained score row has an empty donor/library inference key")
    arm_assignments = combined[observation_columns + ["contrast_arm"]].drop_duplicates()
    conflicting_arms = arm_assignments.duplicated(observation_columns, keep=False)
    if conflicting_arms.any():
        raise ValueError(
            "A donor/library inference key is assigned to both contrast arms within the same "
            "cohort and compartment. Resolve sample metadata or aggregate technical libraries before "
            "donor-level summarization.\n"
            + arm_assignments.loc[conflicting_arms].sort_values(observation_columns).to_csv(index=False)
        )
    duplicate_columns = ["cohort_id", "donor_id", "compartment", "module_id"]
    duplicates = combined.duplicated(duplicate_columns, keep=False)
    if duplicates.any():
        problematic = combined.loc[duplicates, duplicate_columns + ["dataset", "gsm", "_score_file"]]
        raise ValueError(
            "More than one score record for the same donor/library x cohort x compartment x module. "
            "Aggregate technical libraries before this donor-level summary.\n"
            + problematic.to_csv(index=False)
        )
    return combined, audit_rows


def build_group_descriptives(
    scores: pd.DataFrame, confidence_level: float
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    group_columns = [
        "cohort_id", "cohort_role", "analysis_role", "inference_key_type", "compartment", "module_id",
        "severity_orientation", "contrast_label", "severity_group", "contrast_arm",
    ]
    for values, frame in scores.groupby(group_columns, sort=True, dropna=False):
        (
            cohort_id, cohort_role, analysis_role, inference_key_type, compartment, module_id,
            severity_orientation, contrast_label, severity_group, contrast_arm,
        ) = values
        score_values = frame["module_score"].to_numpy(dtype=float)
        lower, upper, ci_status = mean_ci(score_values, confidence_level)
        cells = optional_numeric_column(frame, "included_cells")
        umi = optional_numeric_column(frame, "total_umi_included_cells")
        rows.append(
            {
                "cohort_id": cohort_id,
                "cohort_role": cohort_role,
                "analysis_role": analysis_role,
                "inference_key_type": inference_key_type,
                "compartment": compartment,
                "module_id": module_id,
                "severity_orientation": severity_orientation,
                "contrast_label": contrast_label,
                "severity_group": severity_group,
                "contrast_arm": contrast_arm,
                "n_donor_or_library_keys": int(len(frame)),
                "mean_module_score_log1p_cpm": numeric_or_blank(float(np.mean(score_values))),
                "sd_module_score_log1p_cpm": numeric_or_blank(safe_sample_sd(score_values)),
                "median_module_score_log1p_cpm": numeric_or_blank(float(np.median(score_values))),
                "min_module_score_log1p_cpm": numeric_or_blank(float(np.min(score_values))),
                "max_module_score_log1p_cpm": numeric_or_blank(float(np.max(score_values))),
                "mean_ci_lower": numeric_or_blank(lower),
                "mean_ci_upper": numeric_or_blank(upper),
                "mean_ci_status": ci_status,
                "median_included_cells": numeric_or_blank(float(cells.median()) if cells.notna().any() else None, 3),
                "min_included_cells": numeric_or_blank(float(cells.min()) if cells.notna().any() else None, 3),
                "max_included_cells": numeric_or_blank(float(cells.max()) if cells.notna().any() else None, 3),
                "median_total_umi_included_cells": numeric_or_blank(float(umi.median()) if umi.notna().any() else None, 3),
                "analysis_note": "Scores are unweighted donor/library observations; nested cells are not treated as replicates.",
            }
        )
    return rows


def build_effects_and_lodo(
    scores: pd.DataFrame,
    confidence_level: float,
    n_bootstrap: int,
    random_seed: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    effect_rows: list[dict[str, Any]] = []
    lodo_rows: list[dict[str, Any]] = []
    grouping = [
        "cohort_id", "cohort_role", "analysis_role", "inference_key_type", "compartment", "module_id",
        "comparison_label", "target_label", "contrast_label", "severity_orientation", "confirmatory_eligible", "notes",
    ]
    for values, frame in scores.groupby(grouping, sort=True, dropna=False):
        (
            cohort_id, cohort_role, analysis_role, inference_key_type, compartment, module_id,
            comparison_label, target_label, contrast_label, severity_orientation, confirmatory_eligible, notes,
        ) = values
        comparison_frame = frame.loc[frame["contrast_arm"].eq("comparison")].copy()
        target_frame = frame.loc[frame["contrast_arm"].eq("target")].copy()
        comparison = comparison_frame["module_score"].to_numpy(dtype=float)
        target = target_frame["module_score"].to_numpy(dtype=float)
        if len(comparison) == 0 or len(target) == 0:
            effect = None
            effect_status = "not_estimable_one_or_both_contrast_groups_empty"
        else:
            effect = float(np.mean(target) - np.mean(comparison))
            effect_status = "available_unweighted_target_minus_comparison_mean_difference"
        analytic_lower, analytic_upper, analytic_se, analytic_df, analytic_status = analytic_difference_ci(
            target, comparison, confidence_level
        )
        effect_seed = stable_seed(random_seed, cohort_id, compartment, module_id)
        bootstrap_lower, bootstrap_upper, bootstrap_status = bootstrap_difference_ci(
            target, comparison, confidence_level, n_bootstrap, effect_seed
        )

        lodo_full_direction = direction(effect)
        lodo_matching = 0
        lodo_direction_defined = 0
        lodo_available = 0
        lodo_effects: list[float] = []
        full_frame = pd.concat([comparison_frame, target_frame], ignore_index=True)
        for _, excluded in full_frame.sort_values(["contrast_arm", "donor_id"]).iterrows():
            retained_comparison = comparison_frame.loc[
                comparison_frame["donor_id"].ne(excluded["donor_id"]), "module_score"
            ].to_numpy(dtype=float)
            retained_target = target_frame.loc[
                target_frame["donor_id"].ne(excluded["donor_id"]), "module_score"
            ].to_numpy(dtype=float)
            if len(retained_comparison) == 0 or len(retained_target) == 0:
                lodo_effect = None
                lodo_status = "not_estimable_one_or_both_remaining_groups_empty"
            else:
                lodo_effect = float(np.mean(retained_target) - np.mean(retained_comparison))
                lodo_status = "available_unweighted_target_minus_comparison_mean_difference"
                lodo_available += 1
                lodo_effects.append(lodo_effect)
            lodo_direction = direction(lodo_effect)
            is_matching: bool | None
            if lodo_full_direction in {"not_estimable", "zero"} or lodo_direction in {"not_estimable", "zero"}:
                is_matching = None
            else:
                lodo_direction_defined += 1
                is_matching = lodo_direction == lodo_full_direction
                lodo_matching += int(is_matching)
            lodo_rows.append(
                {
                    "cohort_id": cohort_id,
                    "cohort_role": cohort_role,
                    "analysis_role": analysis_role,
                    "inference_key_type": inference_key_type,
                    "compartment": compartment,
                    "module_id": module_id,
                    "contrast_label": contrast_label,
                    "excluded_donor_or_library_key": excluded["donor_id"],
                    "excluded_arm": excluded["contrast_arm"],
                    "remaining_comparison_n": len(retained_comparison),
                    "remaining_target_n": len(retained_target),
                    "leave_one_out_effect_target_minus_comparison": numeric_or_blank(lodo_effect),
                    "leave_one_out_effect_direction": lodo_direction,
                    "full_effect_direction": lodo_full_direction,
                    "direction_agrees_with_full_effect": "" if is_matching is None else str(is_matching).lower(),
                    "leave_one_out_status": lodo_status,
                }
            )
        if lodo_full_direction in {"not_estimable", "zero"}:
            retention_fraction = None
            lodo_stability_status = "not_assessable_full_effect_zero_or_nonfinite"
        elif lodo_direction_defined == 0:
            retention_fraction = None
            lodo_stability_status = "not_assessable_no_nonzero_leave_one_out_effect"
        else:
            retention_fraction = lodo_matching / lodo_direction_defined
            lodo_stability_status = (
                "direction_retained_at_or_above_0.80"
                if retention_fraction >= 0.80
                else "direction_retention_below_0.80"
            )
        effect_rows.append(
            {
                "cohort_id": cohort_id,
                "cohort_role": cohort_role,
                "analysis_role": analysis_role,
                "inference_key_type": inference_key_type,
                "compartment": compartment,
                "module_id": module_id,
                "comparison_label": comparison_label,
                "target_label": target_label,
                "contrast_label": contrast_label,
                "severity_orientation": severity_orientation,
                "comparison_n": len(comparison),
                "target_n": len(target),
                "comparison_mean_module_score_log1p_cpm": numeric_or_blank(float(np.mean(comparison)) if len(comparison) else None),
                "target_mean_module_score_log1p_cpm": numeric_or_blank(float(np.mean(target)) if len(target) else None),
                "mean_difference_target_minus_comparison": numeric_or_blank(effect),
                "effect_direction": direction(effect),
                "effect_status": effect_status,
                "analytic_ci_method": "Welch_t_interval_for_unweighted_donor_or_library_mean_difference",
                "analytic_ci_lower": numeric_or_blank(analytic_lower),
                "analytic_ci_upper": numeric_or_blank(analytic_upper),
                "analytic_standard_error": numeric_or_blank(analytic_se),
                "analytic_welch_df": numeric_or_blank(analytic_df),
                "analytic_ci_status": analytic_status,
                "bootstrap_ci_method": "independent_percentile_donor_or_library_bootstrap",
                "bootstrap_ci_lower": numeric_or_blank(bootstrap_lower),
                "bootstrap_ci_upper": numeric_or_blank(bootstrap_upper),
                "bootstrap_ci_status": bootstrap_status,
                "bootstrap_replicates": n_bootstrap,
                "bootstrap_seed": effect_seed,
                "leave_one_out_total_keys": len(full_frame),
                "leave_one_out_available_effects": lodo_available,
                "leave_one_out_nonzero_directional_effects": lodo_direction_defined,
                "leave_one_out_matching_full_direction": lodo_matching,
                "leave_one_out_direction_retention_fraction": numeric_or_blank(retention_fraction),
                "leave_one_out_stability_status": lodo_stability_status,
                "confirmatory_eligible": confirmatory_eligible,
                "cohort_limitations": notes,
                "inference_boundary": (
                    "Descriptive donor/library-level effect only; no cell-level test, p-value, causal, "
                    "age-independent, or therapeutic interpretation."
                ),
            }
        )
    return effect_rows, lodo_rows


def build_direction_consistency(effect_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Describe effect-sign alignment across independently labelled cohort IDs."""
    rows: list[dict[str, Any]] = []
    grouped: defaultdict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in effect_rows:
        grouped[(row["module_id"], row["compartment"], row["severity_orientation"])].append(row)
    for (module_id, compartment, severity_orientation), effects in sorted(grouped.items()):
        available = [row for row in effects if row["effect_status"].startswith("available")]
        positive = sum(row["effect_direction"] == "positive_higher_target_severity" for row in available)
        negative = sum(row["effect_direction"] == "negative_lower_target_severity" for row in available)
        zero = sum(row["effect_direction"] == "zero" for row in available)
        nonzero = positive + negative
        if nonzero < 2:
            consensus_direction = "not_assessable_fewer_than_two_nonzero_cohort_effects"
            consistency_status = "insufficient_cohort_effects_for_direction_alignment"
        elif positive == nonzero:
            consensus_direction = "positive_higher_target_severity"
            consistency_status = "available_effects_directionally_aligned"
        elif negative == nonzero:
            consensus_direction = "negative_lower_target_severity"
            consistency_status = "available_effects_directionally_aligned"
        else:
            consensus_direction = "discordant"
            consistency_status = "available_effects_directionally_discordant"
        eligible = [row for row in available if parse_bool(row["confirmatory_eligible"], "confirmatory_eligible")]
        if len(eligible) < 2:
            replication_status = (
                "not_assessed_fewer_than_two_confirmatory_eligible_independent_validation_contrasts; "
                "no_meta_analysis_run"
            )
        else:
            replication_status = "not_adjudicated_by_this_script_no_meta_analysis_run"
        effect_pairs = ";".join(
            f"{row['cohort_id']}={row['mean_difference_target_minus_comparison']}" for row in available
        )
        rows.append(
            {
                "module_id": module_id,
                "compartment": compartment,
                "severity_orientation": severity_orientation,
                "n_independent_cohort_effects_available": len(available),
                "n_positive": positive,
                "n_negative": negative,
                "n_zero": zero,
                "direction_consistency_fraction_among_nonzero": numeric_or_blank(
                    max(positive, negative) / nonzero if nonzero else None
                ),
                "consensus_direction": consensus_direction,
                "direction_consistency_status": consistency_status,
                "cohort_ids": ";".join(row["cohort_id"] for row in available),
                "cohort_roles": ";".join(row["cohort_role"] for row in available),
                "cohort_effects_target_minus_comparison": effect_pairs,
                "confirmatory_eligible_validation_effects": len(eligible),
                "replication_assessment": replication_status,
                "interpretation_boundary": (
                    "Sign alignment is descriptive. It is not a random-effects meta-analysis, does not "
                    "establish replication, and does not upgrade exploratory or directional cohorts."
                ),
            }
        )
    return rows


def write_readme(
    output_dir: Path,
    confidence_level: float,
    n_bootstrap: int,
    min_mapped_fraction: float,
    effect_rows: list[dict[str, Any]],
    availability_written: bool,
    ledger_crosswalk_written: bool,
    summary_contract: str,
) -> None:
    available = sum(row["effect_status"].startswith("available") for row in effect_rows)
    contract = SUMMARY_CONTRACTS[summary_contract]
    text = f"""# {contract["title"]}

Generated by `scripts/summarize_donor_module_effects.py`.

## Result contract

{contract["boundary"]}

## Estimand

Each effect is the unweighted mean module-score difference of target minus
comparison severity group, with one donor or presumed donor/library key as one
observation. Cells are nested observations and were not used as independent
replicates. This directory itself includes no p-values, cell-level tests,
meta-analysis, or replication adjudication. Downstream S7/S8/S9 exploratory
meta-analyses are separately packaged and do not change this directory's contract.

## Current run

- Available cohort/compartment/module effects: {available}
- Confidence level: {confidence_level:.3f}
- Bootstrap replicates: {n_bootstrap}
- Minimum mapped-gene fraction accepted from score table: {min_mapped_fraction:.3f}
"""
    if summary_contract == "default_current_project":
        text += """
- GSE230809 is one exploratory parent project and is not split into discovery
  and validation cohorts.
- GSE244889 remains a directional 4-versus-3 presumed donor/library-key
  comparison; it is not confirmatory validation.
- GSE153066 is represented as an external NP count-level support contrast,
  with the prior contributor cell filtering and clinical-source/age
  confounding explicitly retained.
- GSE165722 is represented only as a 4-versus-4 NP normalized-count
  score-level direction/stability contrast. Its source-publication severity
  grouping is used because it differs from GEO grades; it supplies no
  raw-count or confirmatory inference.
- GSE251686 has a separate audited exploratory mild n=2 versus severe n=3
  score/effect package, but is deliberately excluded from this default
  20-effect cross-cohort summary. Its permanently malformed `GSM7986002`
  record remains excluded and the separate package is not validation,
  replication, or confirmatory evidence.
"""
    elif summary_contract == "post_hoc_external_expansion":
        text += """
- GSE186542 and the GSE167931 FPKM representation are post hoc additions.
- GSE167931 TPM is the same-sample processing sensitivity and is not counted as
  another cohort.
"""
    elif summary_contract == "source_family_replacement":
        text += """
- GSE245147 replaces GSE167931 for source-family sensitivity only; the two are
  never pooled as independent cohorts.
- Only the native GSE245147 Degenerated versus No-degenerated columns enter;
  passage and treatment arms remain excluded.
"""
    text += """

## Current-run artifacts

Only the files listed in `run_artifacts.csv` belong to this invocation.  If
the output directory was reused, any unlisted file may be from an older run
and must not be interpreted as part of this result set.

- `resolved_contrast_spec.csv`: exact cohort/condition map used.
- `run_parameters.csv`: input score-table SHA-256 values and summary settings.
- `score_ingestion_audit.csv`: retained and excluded score-table rows.
- `donor_module_group_descriptives.csv`: group-level donor/library score summaries.
- `donor_module_effects.csv`: target-minus-comparison effects and two descriptive CI forms.
- `donor_module_leave_one_out.csv`: leave-one-donor/library-key effects.
- `cross_cohort_direction_consistency.csv`: descriptive sign alignment only.
- `run_manifest.json`: result contract, exact contrast map, input
  hashes, and generated-artifact hashes (excluding the self-referential run
  ledger).
"""
    if availability_written:
        text += """
- `cohort_score_availability.csv`: audited current-project cohorts and whether
  a compatible score table was available; absent tables are not interpreted as
  null effects.
"""
    if ledger_crosswalk_written:
        text += """
- `score_ledger_identity_crosswalk.csv`: exact score-to-sample-ledger identity
  check for the current fixed-project inputs.
"""
    text += """
- `run_artifacts.csv`: machine-readable list of artifacts generated in this invocation.

Confidence intervals can be unavailable for groups with fewer than two keys or
zero/nonfinite variance; the corresponding status fields state why. The
analytical interval is Welch t based and the bootstrap interval resamples
donor/library keys independently within the two groups.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--use-current-project-inputs",
        action="store_true",
        help=(
            "use the fixed map for currently available GSE230809, GSE244889, "
            "GSE153066, and GSE165722 score tables"
        ),
    )
    source.add_argument(
        "--contrast-spec",
        type=Path,
        help="CSV with the exact columns written to resolved_contrast_spec.csv",
    )
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--confidence-level", type=float, default=0.95)
    parser.add_argument("--bootstrap-replicates", type=int, default=10_000)
    parser.add_argument("--random-seed", type=int, default=20_260_814)
    parser.add_argument("--min-mapped-fraction", type=float, default=0.80)
    parser.add_argument(
        "--summary-contract",
        choices=sorted(SUMMARY_CONTRACTS),
        help=(
            "result-directory contract; defaults to default_current_project for the fixed current "
            "inputs and custom_sensitivity for an explicit contrast specification"
        ),
    )
    args = parser.parse_args()
    if not 0 < args.confidence_level < 1:
        raise ValueError("--confidence-level must be between 0 and 1")
    if args.bootstrap_replicates < 0:
        raise ValueError("--bootstrap-replicates cannot be negative")
    if not 0 < args.min_mapped_fraction <= 1:
        raise ValueError("--min-mapped-fraction must be in (0, 1]")

    summary_contract = args.summary_contract or (
        "default_current_project" if args.use_current_project_inputs else "custom_sensitivity"
    )

    project_root = args.project_root.resolve()
    if args.use_current_project_inputs:
        specs = current_project_spec()
        spec_directory = None
        availability_rows = current_project_score_availability()
    else:
        spec_path = args.contrast_spec.resolve()
        specs = load_contrast_spec(spec_path)
        spec_directory = spec_path.parent
        availability_rows = []
    validate_spec(specs)
    scores, ingestion_audit = collect_scores(specs, project_root, spec_directory, args.min_mapped_fraction)
    ledger_crosswalk_rows = (
        current_project_ledger_crosswalk(specs, scores, project_root, spec_directory)
        if args.use_current_project_inputs
        else []
    )
    group_rows = build_group_descriptives(scores, args.confidence_level)
    effect_rows, lodo_rows = build_effects_and_lodo(
        scores, args.confidence_level, args.bootstrap_replicates, args.random_seed
    )
    direction_rows = build_direction_consistency(effect_rows)

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    generated_artifacts = [
        "README.md",
        "resolved_contrast_spec.csv",
        "run_parameters.csv",
        "score_ingestion_audit.csv",
        "donor_module_group_descriptives.csv",
        "donor_module_effects.csv",
        "donor_module_leave_one_out.csv",
        "cross_cohort_direction_consistency.csv",
        "run_artifacts.csv",
        "run_manifest.json",
    ]
    if availability_rows:
        generated_artifacts.append("cohort_score_availability.csv")
    if ledger_crosswalk_rows:
        generated_artifacts.append("score_ledger_identity_crosswalk.csv")
    rows_to_csv(output_dir / "resolved_contrast_spec.csv", specs, SPEC_COLUMNS)
    if availability_rows:
        rows_to_csv(
            output_dir / "cohort_score_availability.csv",
            availability_rows,
            [
                "cohort_id", "dataset_or_parent_project", "compartment_scope", "sample_key_type",
                "audited_group_structure", "module_score_availability", "score_table_paths",
                "effect_summary_included", "analysis_boundary",
            ],
        )
    if ledger_crosswalk_rows:
        rows_to_csv(
            output_dir / "score_ledger_identity_crosswalk.csv",
            ledger_crosswalk_rows,
            [
                "cohort_id", "ledger_label", "ledger_path", "identity_columns", "score_group_column",
                "ledger_group_column", "unique_score_table_sample_keys",
                "unique_retained_score_sample_keys", "unique_ledger_sample_keys", "matched_keys",
                "score_table_only_keys", "ledger_only_keys", "crosswalk_status",
            ],
        )
    parameter_rows = run_parameter_rows(
        specs,
        project_root,
        spec_directory,
        args.confidence_level,
        args.bootstrap_replicates,
        args.random_seed,
        args.min_mapped_fraction,
        summary_contract,
    )
    rows_to_csv(
        output_dir / "run_parameters.csv",
        parameter_rows,
        ["parameter_class", "parameter_name", "parameter_value", "sha256", "notes"],
    )
    rows_to_csv(
        output_dir / "score_ingestion_audit.csv",
        ingestion_audit,
        [
            "cohort_id", "compartment", "spec_index", "source_score_files", "raw_rows_loaded",
            "rows_matching_dataset", "rows_matching_compartment", "rows_matching_contrast_groups",
            "row_disposition", "reason", "row_count",
        ],
    )
    rows_to_csv(
        output_dir / "donor_module_group_descriptives.csv",
        group_rows,
        [
            "cohort_id", "cohort_role", "analysis_role", "inference_key_type", "compartment", "module_id",
            "severity_orientation", "contrast_label", "severity_group", "contrast_arm", "n_donor_or_library_keys",
            "mean_module_score_log1p_cpm", "sd_module_score_log1p_cpm", "median_module_score_log1p_cpm",
            "min_module_score_log1p_cpm", "max_module_score_log1p_cpm", "mean_ci_lower", "mean_ci_upper",
            "mean_ci_status", "median_included_cells", "min_included_cells", "max_included_cells",
            "median_total_umi_included_cells", "analysis_note",
        ],
    )
    rows_to_csv(
        output_dir / "donor_module_effects.csv",
        effect_rows,
        [
            "cohort_id", "cohort_role", "analysis_role", "inference_key_type", "compartment", "module_id",
            "comparison_label", "target_label", "contrast_label", "severity_orientation", "comparison_n", "target_n",
            "comparison_mean_module_score_log1p_cpm", "target_mean_module_score_log1p_cpm",
            "mean_difference_target_minus_comparison", "effect_direction", "effect_status", "analytic_ci_method",
            "analytic_ci_lower", "analytic_ci_upper", "analytic_standard_error", "analytic_welch_df",
            "analytic_ci_status", "bootstrap_ci_method", "bootstrap_ci_lower", "bootstrap_ci_upper",
            "bootstrap_ci_status", "bootstrap_replicates", "bootstrap_seed", "leave_one_out_total_keys",
            "leave_one_out_available_effects", "leave_one_out_nonzero_directional_effects",
            "leave_one_out_matching_full_direction", "leave_one_out_direction_retention_fraction",
            "leave_one_out_stability_status", "confirmatory_eligible", "cohort_limitations", "inference_boundary",
        ],
    )
    rows_to_csv(
        output_dir / "donor_module_leave_one_out.csv",
        lodo_rows,
        [
            "cohort_id", "cohort_role", "analysis_role", "inference_key_type", "compartment", "module_id",
            "contrast_label", "excluded_donor_or_library_key", "excluded_arm", "remaining_comparison_n",
            "remaining_target_n", "leave_one_out_effect_target_minus_comparison", "leave_one_out_effect_direction",
            "full_effect_direction", "direction_agrees_with_full_effect", "leave_one_out_status",
        ],
    )
    rows_to_csv(
        output_dir / "cross_cohort_direction_consistency.csv",
        direction_rows,
        [
            "module_id", "compartment", "severity_orientation", "n_independent_cohort_effects_available",
            "n_positive", "n_negative", "n_zero", "direction_consistency_fraction_among_nonzero",
            "consensus_direction", "direction_consistency_status", "cohort_ids", "cohort_roles",
            "cohort_effects_target_minus_comparison", "confirmatory_eligible_validation_effects",
            "replication_assessment", "interpretation_boundary",
        ],
    )
    write_readme(
        output_dir,
        args.confidence_level,
        args.bootstrap_replicates,
        args.min_mapped_fraction,
        effect_rows,
        bool(availability_rows),
        bool(ledger_crosswalk_rows),
        summary_contract,
    )
    rows_to_csv(
        output_dir / "run_artifacts.csv",
        [
            {
                "artifact": name,
                "generated_by_this_invocation": "true",
                "purpose": "See README.md for the artifact contract.",
            }
            for name in generated_artifacts
        ],
        ["artifact", "generated_by_this_invocation", "purpose"],
    )
    write_run_manifest(
        output_dir,
        generated_artifacts,
        specs,
        parameter_rows,
        effect_rows,
        summary_contract,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
