"""Summarize donor/library-level retained-cell threshold sensitivity runs.

This script does not ingest cell-level expression.  It only compares the
already aggregated module scores and donor/library effect summaries generated
for the pre-specified 20, 30, and 50 source-restricted-cell thresholds.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path
from typing import Iterable

import pandas as pd


THRESHOLDS = (20, 30, 50)
SCORE_KEY = ["dataset", "gsm", "donor_id", "compartment", "disease_state", "module_id"]
EFFECT_KEY = ["cohort_id", "compartment", "module_id"]


def write_csv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_columns(frame: pd.DataFrame, columns: list[str], path: Path) -> None:
    missing = set(columns) - set(frame.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")


def write_readme(output_dir: Path) -> None:
    """Document the threshold-analysis contract next to its generated tables."""
    text = """# Discovery Retained-Cell Threshold Sensitivity

This directory summarizes the pre-specified 20, 30, and 50 retained-cell
threshold analysis for the combined GSE230809 discovery parent project
(GSE229711 plus GSE230808). The threshold is a donor/library eligibility rule:
each source-restricted, QC-passing donor/library must contain at least the
specified number of cells before its full pseudobulk is scored. It is not a
random cell downsampling analysis and it does not make cells independent
replicates.

For each threshold, `score_module_pseudobulk.py` was rerun from the original
10x TAR archives using the same locked module configuration, QC ledger, and
annotation ledger. `summarize_donor_module_effects.py` then calculated only
unweighted donor/library target-minus-comparison score differences and
leave-one-donor/library-out displays for AF and NP separately.

Read `threshold_run_summary.csv` before interpreting the effect tables.
`threshold_score_identity.csv` compares every threshold-specific module score
to the primary 30-cell run. `threshold_effect_stability_vs_30.csv` compares
the resulting donor/library effects and directions. All effects remain
exploratory because GSE230809 is one parent project with three healthy donors
per compartment and complete age-disease confounding.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    input_root = args.input_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    score_frames: list[pd.DataFrame] = []
    effect_frames: list[pd.DataFrame] = []
    eligibility_rows: list[dict[str, object]] = []
    run_rows: list[dict[str, object]] = []
    artifact_rows: list[dict[str, object]] = []

    for threshold in THRESHOLDS:
        threshold_root = input_root / f"threshold_{threshold}"
        summary_path = threshold_root / "effect_summary" / "donor_module_effects.csv"
        effect = pd.read_csv(summary_path, dtype=str)
        require_columns(
            effect,
            EFFECT_KEY
            + [
                "comparison_n",
                "target_n",
                "mean_difference_target_minus_comparison",
                "effect_direction",
                "effect_status",
                "leave_one_out_direction_retention_fraction",
            ],
            summary_path,
        )
        effect["retained_cell_threshold"] = threshold
        effect_frames.append(effect)

        threshold_score_paths = []
        threshold_library_paths = []
        for child_series in ("GSE229711", "GSE230808"):
            series_dir = threshold_root / child_series
            score_path = series_dir / f"{child_series}_RAW_module_scores.csv"
            ledger_path = series_dir / f"{child_series}_RAW_library_pseudobulk_ledger.csv"
            score = pd.read_csv(score_path, dtype=str)
            require_columns(score, SCORE_KEY + ["module_score_log1p_cpm", "included_cells"], score_path)
            score["retained_cell_threshold"] = threshold
            score_frames.append(score)
            threshold_score_paths.append(score_path)
            threshold_library_paths.append(ledger_path)

            ledger = pd.read_csv(ledger_path, dtype=str)
            require_columns(
                ledger,
                ["dataset", "gsm", "donor_id", "compartment", "disease_state", "included_cells", "total_umi_included_cells"],
                ledger_path,
            )
            for row in ledger.to_dict("records"):
                included_cells = int(row["included_cells"])
                eligibility_rows.append(
                    {
                        "retained_cell_threshold": threshold,
                        "dataset": row["dataset"],
                        "gsm": row["gsm"],
                        "donor_id": row["donor_id"],
                        "compartment": row["compartment"],
                        "disease_state": row["disease_state"],
                        "source_restricted_qc_passing_cells": included_cells,
                        "total_umi_included_cells": row["total_umi_included_cells"],
                        "threshold_pass": str(included_cells >= threshold).lower(),
                        "inference_unit": "donor/library; cells are nested observations",
                    }
                )

        eligible = [row for row in eligibility_rows if row["retained_cell_threshold"] == threshold]
        failed = [row for row in eligible if row["threshold_pass"] != "true"]
        run_rows.append(
            {
                "retained_cell_threshold": threshold,
                "child_series": "GSE229711;GSE230808",
                "libraries_scored": len(eligible),
                "libraries_passing_threshold": len(eligible) - len(failed),
                "libraries_failing_threshold": len(failed),
                "minimum_observed_source_restricted_qc_passing_cells": min(
                    int(row["source_restricted_qc_passing_cells"]) for row in eligible
                ),
                "threshold_run_status": "all_libraries_scored" if not failed else "unexpected_below_threshold_library_in_output",
                "inference_unit": "donor/library; cells are nested observations",
                "effect_summary_path": str(summary_path),
            }
        )
        for path in threshold_score_paths + threshold_library_paths + [summary_path]:
            artifact_rows.append(
                {
                    "retained_cell_threshold": threshold,
                    "artifact": str(path),
                    "sha256": sha256(path),
                    "purpose": "Threshold-specific scorer or donor/library effect-summary input.",
                }
            )

    score_all = pd.concat(score_frames, ignore_index=True)
    for column in ["module_score_log1p_cpm", "included_cells"]:
        score_all[column] = pd.to_numeric(score_all[column], errors="raise")
    base_scores = score_all.loc[score_all["retained_cell_threshold"].eq(30), SCORE_KEY + ["module_score_log1p_cpm", "included_cells"]].copy()
    if base_scores.duplicated(SCORE_KEY).any():
        raise ValueError("Threshold 30 score table has duplicate library/module keys")
    identity_rows: list[dict[str, object]] = []
    for threshold in THRESHOLDS:
        current = score_all.loc[
            score_all["retained_cell_threshold"].eq(threshold), SCORE_KEY + ["module_score_log1p_cpm", "included_cells"]
        ].copy()
        if current.duplicated(SCORE_KEY).any():
            raise ValueError(f"Threshold {threshold} score table has duplicate library/module keys")
        merged = base_scores.merge(current, on=SCORE_KEY, how="outer", suffixes=("_threshold_30", "_current"), indicator=True)
        matched = merged.loc[merged["_merge"].eq("both")].copy()
        if matched.empty:
            max_delta = None
            included_cells_match = False
        else:
            deltas = (matched["module_score_log1p_cpm_threshold_30"] - matched["module_score_log1p_cpm_current"]).abs()
            max_delta = float(deltas.max())
            included_cells_match = bool(
                matched["included_cells_threshold_30"].eq(matched["included_cells_current"]).all()
            )
        identity_rows.append(
            {
                "retained_cell_threshold": threshold,
                "reference_threshold": 30,
                "score_rows_reference_30": len(base_scores),
                "score_rows_current": len(current),
                "matched_score_rows": len(matched),
                "reference_only_score_rows": int(merged["_merge"].eq("left_only").sum()),
                "current_only_score_rows": int(merged["_merge"].eq("right_only").sum()),
                "maximum_absolute_module_score_delta_vs_30": "" if max_delta is None else f"{max_delta:.12g}",
                "included_cell_counts_identical_vs_30": str(included_cells_match).lower(),
                "score_identity_status": (
                    "identical_score_keys_and_values" if len(matched) == len(base_scores) == len(current) and max_delta == 0.0 and included_cells_match
                    else "different_score_keys_or_values"
                ),
            }
        )

    effects_all = pd.concat(effect_frames, ignore_index=True)
    for column in ["comparison_n", "target_n", "mean_difference_target_minus_comparison"]:
        effects_all[column] = pd.to_numeric(effects_all[column], errors="coerce")
    base_effects = effects_all.loc[
        effects_all["retained_cell_threshold"].eq(30),
        EFFECT_KEY + ["comparison_n", "target_n", "mean_difference_target_minus_comparison", "effect_direction", "effect_status", "leave_one_out_direction_retention_fraction"],
    ].copy()
    if base_effects.duplicated(EFFECT_KEY).any():
        raise ValueError("Threshold 30 effect summary has duplicate cohort/compartment/module keys")
    effect_stability_rows: list[dict[str, object]] = []
    for threshold in THRESHOLDS:
        current = effects_all.loc[
            effects_all["retained_cell_threshold"].eq(threshold),
            EFFECT_KEY + ["comparison_n", "target_n", "mean_difference_target_minus_comparison", "effect_direction", "effect_status", "leave_one_out_direction_retention_fraction"],
        ].copy()
        merged = base_effects.merge(current, on=EFFECT_KEY, how="outer", suffixes=("_threshold_30", "_current"), indicator=True)
        for row in merged.to_dict("records"):
            if row["_merge"] != "both":
                effect_delta = ""
                directions_agree = ""
            else:
                reference_effect = row["mean_difference_target_minus_comparison_threshold_30"]
                current_effect = row["mean_difference_target_minus_comparison_current"]
                effect_delta = "" if pd.isna(reference_effect) or pd.isna(current_effect) else f"{current_effect - reference_effect:.12g}"
                directions_agree = str(row["effect_direction_threshold_30"] == row["effect_direction_current"]).lower()
            effect_stability_rows.append(
                {
                    "retained_cell_threshold": threshold,
                    "reference_threshold": 30,
                    "cohort_id": row["cohort_id"],
                    "compartment": row["compartment"],
                    "module_id": row["module_id"],
                    "row_match_status": row["_merge"],
                    "comparison_n_threshold_30": row.get("comparison_n_threshold_30", ""),
                    "target_n_threshold_30": row.get("target_n_threshold_30", ""),
                    "comparison_n_current": row.get("comparison_n_current", ""),
                    "target_n_current": row.get("target_n_current", ""),
                    "effect_threshold_30": row.get("mean_difference_target_minus_comparison_threshold_30", ""),
                    "effect_current": row.get("mean_difference_target_minus_comparison_current", ""),
                    "effect_delta_current_minus_30": effect_delta,
                    "effect_direction_threshold_30": row.get("effect_direction_threshold_30", ""),
                    "effect_direction_current": row.get("effect_direction_current", ""),
                    "effect_directions_agree": directions_agree,
                    "effect_status_threshold_30": row.get("effect_status_threshold_30", ""),
                    "effect_status_current": row.get("effect_status_current", ""),
                    "inference_boundary": "Descriptive donor/library-level comparison; cells are not independent replicates.",
                }
            )

    fields_eligibility = [
        "retained_cell_threshold", "dataset", "gsm", "donor_id", "compartment", "disease_state",
        "source_restricted_qc_passing_cells", "total_umi_included_cells", "threshold_pass", "inference_unit",
    ]
    fields_runs = [
        "retained_cell_threshold", "child_series", "libraries_scored", "libraries_passing_threshold",
        "libraries_failing_threshold", "minimum_observed_source_restricted_qc_passing_cells", "threshold_run_status",
        "inference_unit", "effect_summary_path",
    ]
    fields_identity = [
        "retained_cell_threshold", "reference_threshold", "score_rows_reference_30", "score_rows_current",
        "matched_score_rows", "reference_only_score_rows", "current_only_score_rows",
        "maximum_absolute_module_score_delta_vs_30", "included_cell_counts_identical_vs_30", "score_identity_status",
    ]
    fields_effects = list(effects_all.columns)
    fields_stability = [
        "retained_cell_threshold", "reference_threshold", "cohort_id", "compartment", "module_id", "row_match_status",
        "comparison_n_threshold_30", "target_n_threshold_30", "comparison_n_current", "target_n_current",
        "effect_threshold_30", "effect_current", "effect_delta_current_minus_30", "effect_direction_threshold_30",
        "effect_direction_current", "effect_directions_agree", "effect_status_threshold_30", "effect_status_current",
        "inference_boundary",
    ]
    write_csv(output_dir / "library_threshold_eligibility.csv", eligibility_rows, fields_eligibility)
    write_csv(output_dir / "threshold_run_summary.csv", run_rows, fields_runs)
    write_csv(output_dir / "threshold_score_identity.csv", identity_rows, fields_identity)
    write_csv(output_dir / "discovery_effects_by_retained_cell_threshold.csv", effects_all.to_dict("records"), fields_effects)
    write_csv(output_dir / "threshold_effect_stability_vs_30.csv", effect_stability_rows, fields_stability)
    write_csv(
        output_dir / "input_artifact_hashes.csv",
        artifact_rows,
        ["retained_cell_threshold", "artifact", "sha256", "purpose"],
    )
    write_readme(output_dir)
    write_csv(
        output_dir / "run_artifacts.csv",
        [
            {"artifact": name, "generated_by_this_invocation": "true", "purpose": "See sensitivity documentation."}
            for name in [
                "README.md",
                "library_threshold_eligibility.csv",
                "threshold_run_summary.csv",
                "threshold_score_identity.csv",
                "discovery_effects_by_retained_cell_threshold.csv",
                "threshold_effect_stability_vs_30.csv",
                "input_artifact_hashes.csv",
                "run_artifacts.csv",
            ]
        ],
        ["artifact", "generated_by_this_invocation", "purpose"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
