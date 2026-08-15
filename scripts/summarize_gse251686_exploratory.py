"""Summarize the isolated GSE251686 exploratory module-score contrast.

This script is deliberately separate from the default IVDD cross-cohort
summary. It hard-gates the five stream-integrity-passing presumed
sample/library keys, keeps the malformed GSM7986002 out of every result, and
computes only descriptive severe-minus-mild score differences, uncertainty
intervals, and leave-one-key-out direction checks. It never emits p-values,
pooled estimates, or validation/replication decisions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import t as student_t


EXPECTED_GROUPS = {
    "GSM7986001": "mild",
    "GSM7986003": "mild",
    "GSM7986004": "severe",
    "GSM7986005": "severe",
    "GSM7986006": "severe",
}
EXCLUDED_GSM = "GSM7986002"
MODULE_ORDER = [
    "ecm_collagen_remodeling",
    "inflammatory_nfkb",
    "hypoxia_oxidative_stress",
    "disc_matrix_homeostasis",
]
REQUIRED_SCORE_COLUMNS = {
    "dataset",
    "gsm",
    "presumed_sample_library_key",
    "compartment",
    "severity_group",
    "module_id",
    "module_score_log1p_cpm",
    "score_status",
    "mapped_fraction",
    "included_cells",
    "total_umi_included_cells",
    "analysis_role",
    "confirmatory_eligible",
}
REQUIRED_LEDGER_COLUMNS = {
    "dataset",
    "gsm",
    "presumed_sample_library_key",
    "compartment",
    "severity_group",
    "source_restricted_cells",
    "source_restricted_threshold_20_pass",
    "source_restricted_threshold_30_pass",
    "source_restricted_threshold_50_pass",
    "stream_integrity_pass",
    "identifier_audit_pass",
    "analysis_role",
    "confirmatory_eligible",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    if not rows:
        raise ValueError(f"Refusing to create empty output: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


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


def stable_seed(root_seed: int, module_id: str) -> int:
    payload = f"{root_seed}|GSE251686_exploratory|NP|{module_id}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little", signed=False)


def ci_percentiles(confidence_level: float) -> tuple[float, float]:
    tail = (1.0 - confidence_level) / 2.0
    return tail, 1.0 - tail


def welch_interval(
    target: np.ndarray, comparison: np.ndarray, confidence_level: float
) -> tuple[float | None, float | None, float | None, float | None, str]:
    """Return a descriptive Welch interval without performing a hypothesis test."""
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


def bootstrap_interval(
    target: np.ndarray,
    comparison: np.ndarray,
    confidence_level: float,
    n_bootstrap: int,
    seed: int,
) -> tuple[float | None, float | None, str]:
    """Resample presumed keys independently within arms; never resample cells."""
    if len(target) < 2 or len(comparison) < 2:
        return None, None, "not_estimable_group_n_lt_2"
    if n_bootstrap < 1:
        return None, None, "not_run_n_bootstrap_lt_1"
    rng = np.random.default_rng(seed)
    target_draws = rng.choice(target, size=(n_bootstrap, len(target)), replace=True)
    comparison_draws = rng.choice(comparison, size=(n_bootstrap, len(comparison)), replace=True)
    differences = target_draws.mean(axis=1) - comparison_draws.mean(axis=1)
    if not np.isfinite(differences).all():
        return None, None, "not_estimable_nonfinite_bootstrap_draw"
    lower_q, upper_q = ci_percentiles(confidence_level)
    lower, upper = np.quantile(differences, [lower_q, upper_q], method="linear")
    return float(lower), float(upper), "available_percentile_presumed_key_bootstrap"


def mean_interval(values: np.ndarray, confidence_level: float) -> tuple[float | None, float | None, str]:
    if len(values) < 2:
        return None, None, "not_estimable_group_n_lt_2"
    sd = float(np.std(values, ddof=1))
    if not math.isfinite(sd):
        return None, None, "not_estimable_nonfinite_sd"
    se = sd / math.sqrt(len(values))
    if not math.isfinite(se) or se <= 0:
        return None, None, "not_estimable_zero_or_nonpositive_standard_error"
    critical = float(student_t.ppf((1.0 + confidence_level) / 2.0, len(values) - 1))
    mean = float(np.mean(values))
    return mean - critical * se, mean + critical * se, "available_t"


def normalized_bool(series: pd.Series, label: str) -> pd.Series:
    values = series.astype(str).str.strip().str.lower()
    invalid = ~values.isin({"true", "false"})
    if invalid.any():
        bad = sorted(values.loc[invalid].unique())
        raise ValueError(f"{label} contains non-boolean values: {bad}")
    return values.eq("true")


def load_score_parameters(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, dtype=str)
    required = {"parameter", "value", "sha256"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{path}: missing parameter columns {sorted(missing)}")
    if frame["parameter"].duplicated().any():
        duplicates = sorted(frame.loc[frame["parameter"].duplicated(), "parameter"].unique())
        raise ValueError(f"{path}: duplicate parameter rows {duplicates}")
    return frame.set_index("parameter")["value"].fillna("").astype(str).to_dict()


def load_and_validate_scores(path: Path, min_mapped_fraction: float) -> pd.DataFrame:
    scores = pd.read_csv(path, dtype=str)
    missing = REQUIRED_SCORE_COLUMNS - set(scores.columns)
    if missing:
        raise ValueError(f"{path}: missing score columns {sorted(missing)}")
    scores = scores.copy()
    if not scores["dataset"].eq("GSE251686").all():
        raise ValueError(f"{path}: score table contains a non-GSE251686 dataset")
    if not scores["compartment"].eq("NP").all():
        raise ValueError(f"{path}: score table contains a non-NP compartment")
    if not scores["score_status"].eq("score_available").all():
        raise ValueError(f"{path}: score table contains an unavailable module score")
    if not normalized_bool(scores["confirmatory_eligible"], "score confirmatory_eligible").eq(False).all():
        raise ValueError(f"{path}: exploratory score table cannot contain a confirmatory-eligible row")
    if not scores["analysis_role"].eq("incomplete non-balanced exploratory NP severity direction check only").all():
        raise ValueError(f"{path}: unexpected analysis role")
    if not scores["gsm"].eq(scores["presumed_sample_library_key"]).all():
        raise ValueError(f"{path}: GSM and presumed sample/library key disagree")
    if EXCLUDED_GSM in set(scores["gsm"]):
        raise ValueError(f"{path}: permanently excluded {EXCLUDED_GSM} appears in score rows")
    observed_groups = scores[["gsm", "severity_group"]].drop_duplicates()
    observed_group_map = dict(zip(observed_groups["gsm"], observed_groups["severity_group"], strict=True))
    if observed_group_map != EXPECTED_GROUPS:
        raise ValueError(
            f"{path}: expected exact selected GSE251686 keys/groups {EXPECTED_GROUPS}, got {observed_group_map}"
        )
    if scores.duplicated(["gsm", "module_id"]).any():
        raise ValueError(f"{path}: duplicate GSM/module score row")
    observed_modules = set(scores["module_id"])
    if observed_modules != set(MODULE_ORDER):
        raise ValueError(f"{path}: unexpected module set {sorted(observed_modules)}")
    expected_rows = len(EXPECTED_GROUPS) * len(MODULE_ORDER)
    if len(scores) != expected_rows:
        raise ValueError(f"{path}: expected {expected_rows} score rows, found {len(scores)}")
    scores["module_score"] = pd.to_numeric(scores["module_score_log1p_cpm"], errors="coerce")
    scores["mapped_fraction_numeric"] = pd.to_numeric(scores["mapped_fraction"], errors="coerce")
    scores["included_cells_numeric"] = pd.to_numeric(scores["included_cells"], errors="coerce")
    scores["included_umi_numeric"] = pd.to_numeric(scores["total_umi_included_cells"], errors="coerce")
    if not np.isfinite(scores["module_score"]).all():
        raise ValueError(f"{path}: non-finite module score")
    if not np.isfinite(scores["mapped_fraction_numeric"]).all() or (scores["mapped_fraction_numeric"] < min_mapped_fraction).any():
        raise ValueError(f"{path}: mapping fraction below configured minimum {min_mapped_fraction}")
    if (scores["included_cells_numeric"] <= 0).any() or (scores["included_umi_numeric"] <= 0).any():
        raise ValueError(f"{path}: non-positive included cell or UMI count")
    return scores


def load_and_validate_ledger(path: Path, scores: pd.DataFrame) -> pd.DataFrame:
    ledger = pd.read_csv(path, dtype=str)
    missing = REQUIRED_LEDGER_COLUMNS - set(ledger.columns)
    if missing:
        raise ValueError(f"{path}: missing library-ledger columns {sorted(missing)}")
    ledger = ledger.copy()
    if not ledger["dataset"].eq("GSE251686").all() or not ledger["compartment"].eq("NP").all():
        raise ValueError(f"{path}: ledger has an unexpected dataset or compartment")
    if not ledger["gsm"].eq(ledger["presumed_sample_library_key"]).all():
        raise ValueError(f"{path}: GSM and presumed sample/library key disagree")
    if ledger["gsm"].duplicated().any():
        raise ValueError(f"{path}: duplicate GSM rows")
    observed_group_map = dict(zip(ledger["gsm"], ledger["severity_group"], strict=True))
    if observed_group_map != EXPECTED_GROUPS:
        raise ValueError(f"{path}: expected exact selected GSE251686 keys/groups")
    if not normalized_bool(ledger["confirmatory_eligible"], "ledger confirmatory_eligible").eq(False).all():
        raise ValueError(f"{path}: exploratory ledger cannot contain a confirmatory-eligible row")
    for column in [
        "source_restricted_threshold_20_pass",
        "source_restricted_threshold_30_pass",
        "source_restricted_threshold_50_pass",
        "stream_integrity_pass",
        "identifier_audit_pass",
    ]:
        if not normalized_bool(ledger[column], column).all():
            raise ValueError(f"{path}: all selected keys must pass {column}")
    ledger["source_restricted_cells_numeric"] = pd.to_numeric(
        ledger["source_restricted_cells"], errors="coerce"
    )
    if not np.isfinite(ledger["source_restricted_cells_numeric"]).all() or (ledger["source_restricted_cells_numeric"] < 50).any():
        raise ValueError(f"{path}: selected keys must have at least 50 source-restricted cells")
    score_cells = scores.groupby("gsm", sort=True)["included_cells_numeric"].first().rename("score_cells")
    joined = ledger.set_index("gsm").join(score_cells, how="left")
    if joined["score_cells"].isna().any() or not np.array_equal(
        joined["source_restricted_cells_numeric"].to_numpy(dtype=int),
        joined["score_cells"].to_numpy(dtype=int),
    ):
        raise ValueError(f"{path}: score and library-ledger retained-cell counts disagree")
    return ledger


def build_group_rows(scores: pd.DataFrame, confidence_level: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for module_id in MODULE_ORDER:
        module = scores.loc[scores["module_id"].eq(module_id)]
        for group in ("mild", "severe"):
            values = module.loc[module["severity_group"].eq(group), "module_score"].to_numpy(dtype=float)
            lower, upper, status = mean_interval(values, confidence_level)
            rows.append(
                {
                    "dataset": "GSE251686",
                    "analysis_role": "isolated_incomplete_nonbalanced_exploratory_effect_display_only",
                    "compartment": "NP",
                    "module_id": module_id,
                    "severity_group": group,
                    "presumed_sample_library_key_n": len(values),
                    "mean_module_score_log1p_cpm": numeric_or_blank(float(np.mean(values))),
                    "sample_sd_module_score_log1p_cpm": numeric_or_blank(float(np.std(values, ddof=1))),
                    "minimum_module_score_log1p_cpm": numeric_or_blank(float(np.min(values))),
                    "maximum_module_score_log1p_cpm": numeric_or_blank(float(np.max(values))),
                    "mean_95_ci_lower": numeric_or_blank(lower),
                    "mean_95_ci_upper": numeric_or_blank(upper),
                    "mean_95_ci_status": status,
                    "inference_boundary": "Descriptive presumed sample/library-key summary only; cells are nested observations.",
                }
            )
    return rows


def build_effect_and_loko_rows(
    scores: pd.DataFrame,
    confidence_level: float,
    n_bootstrap: int,
    random_seed: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    effect_rows: list[dict[str, Any]] = []
    loko_rows: list[dict[str, Any]] = []
    for module_id in MODULE_ORDER:
        module = scores.loc[scores["module_id"].eq(module_id)].copy()
        mild = module.loc[module["severity_group"].eq("mild")].sort_values("gsm")
        severe = module.loc[module["severity_group"].eq("severe")].sort_values("gsm")
        comparison = mild["module_score"].to_numpy(dtype=float)
        target = severe["module_score"].to_numpy(dtype=float)
        effect = float(np.mean(target) - np.mean(comparison))
        lower, upper, standard_error, welch_df, welch_status = welch_interval(
            target, comparison, confidence_level
        )
        seed = stable_seed(random_seed, module_id)
        boot_lower, boot_upper, bootstrap_status = bootstrap_interval(
            target, comparison, confidence_level, n_bootstrap, seed
        )
        full_direction = direction(effect)
        matching = 0
        directional = 0
        available = 0
        for _, excluded in module.sort_values(["severity_group", "gsm"]).iterrows():
            retained_mild = mild.loc[mild["gsm"].ne(excluded["gsm"]), "module_score"].to_numpy(dtype=float)
            retained_severe = severe.loc[severe["gsm"].ne(excluded["gsm"]), "module_score"].to_numpy(dtype=float)
            if len(retained_mild) == 0 or len(retained_severe) == 0:
                loko_effect = None
                loko_status = "not_estimable_one_or_both_remaining_groups_empty"
            else:
                loko_effect = float(np.mean(retained_severe) - np.mean(retained_mild))
                loko_status = "available_unweighted_target_minus_comparison_mean_difference"
                available += 1
            loko_direction = direction(loko_effect)
            matches: bool | None = None
            if full_direction not in {"zero", "not_estimable"} and loko_direction not in {"zero", "not_estimable"}:
                directional += 1
                matches = loko_direction == full_direction
                matching += int(matches)
            loko_rows.append(
                {
                    "dataset": "GSE251686",
                    "analysis_role": "isolated_incomplete_nonbalanced_exploratory_effect_display_only",
                    "compartment": "NP",
                    "module_id": module_id,
                    "contrast_label": "severe_minus_mild",
                    "excluded_presumed_sample_library_key": excluded["gsm"],
                    "excluded_arm": excluded["severity_group"],
                    "remaining_mild_n": len(retained_mild),
                    "remaining_severe_n": len(retained_severe),
                    "leave_one_key_out_effect_severe_minus_mild": numeric_or_blank(loko_effect),
                    "leave_one_key_out_direction": loko_direction,
                    "full_effect_direction": full_direction,
                    "direction_agrees_with_full_effect": "" if matches is None else str(matches).lower(),
                    "leave_one_key_out_status": loko_status,
                    "inference_boundary": "Small-sample sensitivity only; not a validation or replication analysis.",
                }
            )
        if full_direction in {"zero", "not_estimable"}:
            retention_fraction = None
            retention_status = "not_assessable_full_effect_zero_or_nonfinite"
        elif directional == 0:
            retention_fraction = None
            retention_status = "not_assessable_no_nonzero_leave_one_key_out_effect"
        else:
            retention_fraction = matching / directional
            retention_status = (
                "direction_retained_at_or_above_0.80"
                if retention_fraction >= 0.80
                else "direction_retention_below_0.80"
            )
        effect_rows.append(
            {
                "dataset": "GSE251686",
                "cohort_role": "isolated_incomplete_nonbalanced_exploratory_check",
                "analysis_role": "isolated_incomplete_nonbalanced_exploratory_effect_display_only",
                "inference_key_type": "presumed_sample_library_key",
                "compartment": "NP",
                "module_id": module_id,
                "comparison_label": "mild",
                "target_label": "severe",
                "contrast_label": "severe_minus_mild",
                "severity_orientation": "higher_severity_minus_lower_severity",
                "comparison_n": len(comparison),
                "target_n": len(target),
                "comparison_mean_module_score_log1p_cpm": numeric_or_blank(float(np.mean(comparison))),
                "target_mean_module_score_log1p_cpm": numeric_or_blank(float(np.mean(target))),
                "mean_difference_target_minus_comparison": numeric_or_blank(effect),
                "effect_direction": full_direction,
                "effect_status": "available_unweighted_target_minus_comparison_mean_difference",
                "analytic_ci_method": "Welch_t_interval_for_unweighted_presumed_sample_library_key_mean_difference",
                "analytic_ci_lower": numeric_or_blank(lower),
                "analytic_ci_upper": numeric_or_blank(upper),
                "analytic_standard_error": numeric_or_blank(standard_error),
                "analytic_welch_df": numeric_or_blank(welch_df),
                "analytic_ci_status": welch_status,
                "bootstrap_ci_method": "independent_percentile_presumed_sample_library_key_bootstrap",
                "bootstrap_ci_lower": numeric_or_blank(boot_lower),
                "bootstrap_ci_upper": numeric_or_blank(boot_upper),
                "bootstrap_ci_status": bootstrap_status,
                "bootstrap_replicates": n_bootstrap,
                "bootstrap_seed": seed,
                "leave_one_key_out_total_keys": len(module),
                "leave_one_key_out_available_effects": available,
                "leave_one_key_out_nonzero_directional_effects": directional,
                "leave_one_key_out_matching_full_direction": matching,
                "leave_one_key_out_direction_retention_fraction": numeric_or_blank(retention_fraction),
                "leave_one_key_out_stability_status": retention_status,
                "confirmatory_eligible": "false",
                "excluded_record": EXCLUDED_GSM,
                "cohort_limitations": (
                    "Two mild and three severe presumed sample/library keys after permanent exclusion of "
                    "GSM7986002 for stream-integrity failure; patient identity, age, sex, individual grade, "
                    "and covariates are unavailable."
                ),
                "inference_boundary": (
                    "Isolated descriptive effect only; cells are nested, no p-value, meta-analysis, validation, "
                    "replication, causal, age-adjusted, biomarker, or therapeutic interpretation."
                ),
            }
        )
    return effect_rows, loko_rows


def write_readme(output_dir: Path, confidence_level: float, n_bootstrap: int, min_mapped_fraction: float) -> None:
    output_dir.joinpath("README.md").write_text(
        "# GSE251686 isolated exploratory effect summary\n\n"
        "This directory is intentionally separate from the default `data/derived/"
        "donor_module_effect_summary/` result. It contains a descriptive mild n=2 "
        "versus severe n=3 comparison of the five stream-integrity-passing presumed "
        "sample/library keys. `GSM7986002` is permanently excluded because its Matrix "
        "Market payload failed the independent stream-integrity audit.\n\n"
        "The summary uses unweighted severe-minus-mild differences, descriptive Welch "
        f"{confidence_level:.0%} intervals, {n_bootstrap} within-arm presumed-key bootstrap "
        "intervals, and leave-one-key-out direction checks. These are not hypothesis "
        "tests: no p-values, multiple-testing adjustment, formal meta-analysis, "
        "validation, replication adjudication, causal claim, biomarker claim, or "
        "therapeutic claim is produced.\n\n"
        f"Every score row must meet the locked mapped-gene fraction of {min_mapped_fraction:.3f}. "
        "The scorer records a 30-cell primary eligibility gate; this effect summary "
        "additionally requires every selected key to pass the 20-, 30-, and 50-cell "
        "source-restricted gates. "
        "Input hashes are recorded in `GSE251686_exploratory_effect_parameters.csv`; "
        "the run manifest records exact generated-artifact hashes. "
        "`run_artifacts.csv` and `run_manifest.json` are excluded from generated-artifact "
        "hashing because each would otherwise be self-referential.\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score-table", type=Path, required=True)
    parser.add_argument("--library-ledger", type=Path, required=True)
    parser.add_argument("--score-parameters", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--confidence-level", type=float, default=0.95)
    parser.add_argument("--bootstrap-replicates", type=int, default=10_000)
    parser.add_argument("--random-seed", type=int, default=20_260_814)
    parser.add_argument("--min-mapped-fraction", type=float, default=0.80)
    args = parser.parse_args()
    if not 0 < args.confidence_level < 1:
        raise ValueError("--confidence-level must be between 0 and 1")
    if args.bootstrap_replicates < 0:
        raise ValueError("--bootstrap-replicates cannot be negative")
    if not 0 < args.min_mapped_fraction <= 1:
        raise ValueError("--min-mapped-fraction must be in (0, 1]")

    score_table = args.score_table.resolve()
    library_ledger = args.library_ledger.resolve()
    score_parameters = args.score_parameters.resolve()
    parameters = load_score_parameters(score_parameters)
    if parameters.get("selected_gsms") != ";".join(EXPECTED_GROUPS):
        raise ValueError("Score parameter ledger does not contain the expected selected GSM order")
    if parameters.get("excluded_gsms") != EXCLUDED_GSM:
        raise ValueError(f"Score parameter ledger does not lock {EXCLUDED_GSM} as excluded")
    if float(parameters.get("min_mapped_fraction", "nan")) != args.min_mapped_fraction:
        raise ValueError("Score parameter mapped-gene threshold differs from effect-summary threshold")
    if parameters.get("inference_boundary", "").find("never confirmatory") < 0:
        raise ValueError("Score parameter ledger lacks the required non-confirmatory boundary")

    scores = load_and_validate_scores(score_table, args.min_mapped_fraction)
    load_and_validate_ledger(library_ledger, scores)
    group_rows = build_group_rows(scores, args.confidence_level)
    effect_rows, loko_rows = build_effect_and_loko_rows(
        scores, args.confidence_level, args.bootstrap_replicates, args.random_seed
    )

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    group_path = output_dir / "GSE251686_exploratory_group_descriptives.csv"
    effect_path = output_dir / "GSE251686_exploratory_module_effects.csv"
    loko_path = output_dir / "GSE251686_exploratory_module_leave_one_out.csv"
    parameter_path = output_dir / "GSE251686_exploratory_effect_parameters.csv"
    artifact_path = output_dir / "run_artifacts.csv"
    manifest_path = output_dir / "run_manifest.json"

    write_csv(
        group_path,
        group_rows,
        [
            "dataset", "analysis_role", "compartment", "module_id", "severity_group",
            "presumed_sample_library_key_n", "mean_module_score_log1p_cpm",
            "sample_sd_module_score_log1p_cpm", "minimum_module_score_log1p_cpm",
            "maximum_module_score_log1p_cpm", "mean_95_ci_lower", "mean_95_ci_upper",
            "mean_95_ci_status", "inference_boundary",
        ],
    )
    write_csv(
        effect_path,
        effect_rows,
        [
            "dataset", "cohort_role", "analysis_role", "inference_key_type", "compartment", "module_id",
            "comparison_label", "target_label", "contrast_label", "severity_orientation", "comparison_n",
            "target_n", "comparison_mean_module_score_log1p_cpm", "target_mean_module_score_log1p_cpm",
            "mean_difference_target_minus_comparison", "effect_direction", "effect_status", "analytic_ci_method",
            "analytic_ci_lower", "analytic_ci_upper", "analytic_standard_error", "analytic_welch_df",
            "analytic_ci_status", "bootstrap_ci_method", "bootstrap_ci_lower", "bootstrap_ci_upper",
            "bootstrap_ci_status", "bootstrap_replicates", "bootstrap_seed", "leave_one_key_out_total_keys",
            "leave_one_key_out_available_effects", "leave_one_key_out_nonzero_directional_effects",
            "leave_one_key_out_matching_full_direction", "leave_one_key_out_direction_retention_fraction",
            "leave_one_key_out_stability_status", "confirmatory_eligible", "excluded_record", "cohort_limitations",
            "inference_boundary",
        ],
    )
    write_csv(
        loko_path,
        loko_rows,
        [
            "dataset", "analysis_role", "compartment", "module_id", "contrast_label",
            "excluded_presumed_sample_library_key", "excluded_arm", "remaining_mild_n", "remaining_severe_n",
            "leave_one_key_out_effect_severe_minus_mild", "leave_one_key_out_direction", "full_effect_direction",
            "direction_agrees_with_full_effect", "leave_one_key_out_status", "inference_boundary",
        ],
    )
    parameter_rows = [
        {
            "parameter_class": "score_input",
            "parameter_name": "GSE251686_exploratory_module_scores.csv",
            "parameter_value": str(score_table),
            "sha256": sha256(score_table),
            "notes": "Exact separate exploratory score table; not a default-summary input.",
        },
        {
            "parameter_class": "score_input",
            "parameter_name": "GSE251686_exploratory_library_ledger.csv",
            "parameter_value": str(library_ledger),
            "sha256": sha256(library_ledger),
            "notes": "Hard-gated presumed sample/library key and retained-cell ledger.",
        },
        {
            "parameter_class": "score_input",
            "parameter_name": "GSE251686_exploratory_score_parameters.csv",
            "parameter_value": str(score_parameters),
            "sha256": sha256(score_parameters),
            "notes": "Locked scorer inputs, permanent exclusion, and scoring boundary.",
        },
        {
            "parameter_class": "analysis_script",
            "parameter_name": "summarize_gse251686_exploratory.py",
            "parameter_value": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
            "notes": "Exact implementation for this isolated descriptive effect display.",
        },
        {
            "parameter_class": "analysis_setting",
            "parameter_name": "confidence_level",
            "parameter_value": f"{args.confidence_level:.8f}",
            "sha256": "",
            "notes": "Descriptive Welch and percentile bootstrap interval level.",
        },
        {
            "parameter_class": "analysis_setting",
            "parameter_name": "bootstrap_replicates",
            "parameter_value": str(args.bootstrap_replicates),
            "sha256": "",
            "notes": "Independent resampling within mild and severe presumed-key arms.",
        },
        {
            "parameter_class": "analysis_setting",
            "parameter_name": "random_seed",
            "parameter_value": str(args.random_seed),
            "sha256": "",
            "notes": "Module-specific deterministic seeds are derived from this root seed.",
        },
        {
            "parameter_class": "analysis_setting",
            "parameter_name": "min_mapped_fraction",
            "parameter_value": f"{args.min_mapped_fraction:.8f}",
            "sha256": "",
            "notes": "Must match the locked scorer threshold.",
        },
        {
            "parameter_class": "analysis_boundary",
            "parameter_name": "default_summary_inclusion",
            "parameter_value": "false",
            "sha256": "",
            "notes": "This directory is excluded from the default 20-effect cross-cohort summary.",
        },
    ]
    write_csv(
        parameter_path,
        parameter_rows,
        ["parameter_class", "parameter_name", "parameter_value", "sha256", "notes"],
    )
    write_readme(output_dir, args.confidence_level, args.bootstrap_replicates, args.min_mapped_fraction)
    artifact_names = [
        "README.md",
        group_path.name,
        effect_path.name,
        loko_path.name,
        parameter_path.name,
        artifact_path.name,
        manifest_path.name,
    ]
    write_csv(
        artifact_path,
        [
            {
                "artifact": name,
                "generated_by_this_invocation": "true",
                "purpose": "See README.md and run_manifest.json for the isolated exploratory contract.",
            }
            for name in artifact_names
        ],
        ["artifact", "generated_by_this_invocation", "purpose"],
    )
    manifest = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_scope": "GSE251686 isolated mild n=2 versus severe n=3 exploratory display",
        "default_cross_cohort_summary_inclusion": False,
        "permanent_exclusion": EXCLUDED_GSM,
        "unit_of_observation": "presumed sample/library key; cells nested",
        "no_hypothesis_tests_or_p_values": True,
        "input_sha256": {
            "GSE251686_exploratory_module_scores.csv": sha256(score_table),
            "GSE251686_exploratory_library_ledger.csv": sha256(library_ledger),
            "GSE251686_exploratory_score_parameters.csv": sha256(score_parameters),
            "summarize_gse251686_exploratory.py": sha256(Path(__file__).resolve()),
        },
        "generated_artifact_sha256": {
            name: sha256(output_dir / name)
            for name in artifact_names
            if name not in {manifest_path.name, artifact_path.name}
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
