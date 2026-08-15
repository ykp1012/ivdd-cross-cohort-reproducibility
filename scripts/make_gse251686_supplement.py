"""Create isolated GSE251686 exploratory sensitivity deliverables.

This script never reads the default 20-effect summary. It accepts only the
separately audited GSE251686 package, verifies its provenance contract, and
writes a supplementary table and figure that retain the non-confirmatory,
non-balanced mild n=2 versus severe n=3 interpretation boundary.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


MODULE_ORDER = [
    "ecm_collagen_remodeling",
    "inflammatory_nfkb",
    "hypoxia_oxidative_stress",
    "disc_matrix_homeostasis",
]
MODULE_LABELS = {
    "ecm_collagen_remodeling": "ECM / collagen remodeling",
    "inflammatory_nfkb": "Inflammatory / NF-kB",
    "hypoxia_oxidative_stress": "Hypoxia / oxidative stress",
    "disc_matrix_homeostasis": "Disc matrix homeostasis",
}
PACKAGE_FILES = {
    "GSE251686_exploratory_module_scores.csv",
    "GSE251686_exploratory_library_ledger.csv",
    "GSE251686_exploratory_score_parameters.csv",
    "GSE251686_exploratory_effect_parameters.csv",
    "GSE251686_exploratory_group_descriptives.csv",
    "GSE251686_exploratory_module_effects.csv",
    "GSE251686_exploratory_module_leave_one_out.csv",
    "README.md",
    "run_manifest.json",
}
EXCLUDED_GSM = "GSM7986002"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def require_isolated_package(package_dir: Path) -> dict:
    missing = sorted(name for name in PACKAGE_FILES if not (package_dir / name).is_file())
    if missing:
        raise FileNotFoundError(f"Missing required GSE251686 package files: {missing}")
    manifest = json.loads((package_dir / "run_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("default_cross_cohort_summary_inclusion") is not False:
        raise ValueError("GSE251686 supplementary source must be excluded from the default summary")
    if manifest.get("permanent_exclusion") != EXCLUDED_GSM:
        raise ValueError("GSE251686 source manifest does not lock the malformed GSM exclusion")
    if manifest.get("unit_of_observation") != "presumed sample/library key; cells nested":
        raise ValueError("GSE251686 source manifest has an unexpected experimental unit")
    if manifest.get("no_hypothesis_tests_or_p_values") is not True:
        raise ValueError("GSE251686 source manifest lacks the no-hypothesis-test boundary")

    expected_inputs = {
        "GSE251686_exploratory_module_scores.csv",
        "GSE251686_exploratory_library_ledger.csv",
        "GSE251686_exploratory_score_parameters.csv",
        "summarize_gse251686_exploratory.py",
    }
    if set(manifest.get("input_sha256", {})) != expected_inputs:
        raise ValueError("GSE251686 source manifest does not contain the expected input hash contract")
    for name in expected_inputs - {"summarize_gse251686_exploratory.py"}:
        if sha256(package_dir / name) != manifest["input_sha256"][name]:
            raise ValueError(f"GSE251686 input hash mismatch: {name}")
    source_script = Path(__file__).with_name("summarize_gse251686_exploratory.py")
    if sha256(source_script) != manifest["input_sha256"]["summarize_gse251686_exploratory.py"]:
        raise ValueError("GSE251686 source summary script hash does not match the package manifest")

    generated = manifest.get("generated_artifact_sha256", {})
    if {"run_artifacts.csv", "run_manifest.json"} & set(generated):
        raise ValueError("Self-referential package ledgers must not be listed as generated-artifact hashes")
    for name, expected_hash in generated.items():
        if sha256(package_dir / name) != expected_hash:
            raise ValueError(f"GSE251686 generated-artifact hash mismatch: {name}")
    return manifest


def load_effects(package_dir: Path) -> pd.DataFrame:
    path = package_dir / "GSE251686_exploratory_module_effects.csv"
    effects = pd.read_csv(path, dtype=str)
    required = {
        "dataset",
        "compartment",
        "module_id",
        "comparison_label",
        "target_label",
        "comparison_n",
        "target_n",
        "mean_difference_target_minus_comparison",
        "analytic_ci_lower",
        "analytic_ci_upper",
        "bootstrap_ci_lower",
        "bootstrap_ci_upper",
        "leave_one_key_out_direction_retention_fraction",
        "confirmatory_eligible",
        "excluded_record",
        "cohort_limitations",
        "inference_boundary",
    }
    missing = required - set(effects.columns)
    if missing:
        raise ValueError(f"GSE251686 effect table is missing columns: {sorted(missing)}")
    if len(effects) != len(MODULE_ORDER):
        raise ValueError("GSE251686 effect table must contain exactly four module effects")
    if set(effects["module_id"]) != set(MODULE_ORDER):
        raise ValueError("GSE251686 effect table does not have the locked module set")
    if not effects["dataset"].eq("GSE251686").all() or not effects["compartment"].eq("NP").all():
        raise ValueError("GSE251686 effect table has an unexpected dataset or compartment")
    if not effects["comparison_label"].eq("mild").all() or not effects["target_label"].eq("severe").all():
        raise ValueError("GSE251686 effect orientation must be severe minus mild")
    if not effects["comparison_n"].eq("2").all() or not effects["target_n"].eq("3").all():
        raise ValueError("GSE251686 effect table must retain mild n=2 and severe n=3")
    if not effects["confirmatory_eligible"].eq("false").all():
        raise ValueError("GSE251686 supplementary results must remain non-confirmatory")
    if not effects["excluded_record"].eq(EXCLUDED_GSM).all():
        raise ValueError("GSE251686 effect table does not retain the permanent GSM exclusion")
    for column in [
        "mean_difference_target_minus_comparison",
        "analytic_ci_lower",
        "analytic_ci_upper",
        "bootstrap_ci_lower",
        "bootstrap_ci_upper",
        "leave_one_key_out_direction_retention_fraction",
    ]:
        effects[column] = pd.to_numeric(effects[column], errors="raise")
    return effects.set_index("module_id").loc[MODULE_ORDER].reset_index()


def write_supplementary_table(effects: pd.DataFrame, output: Path) -> None:
    rows: list[dict[str, str]] = []
    for _, row in effects.iterrows():
        rows.append(
            {
                "Dataset and analysis role": "GSE251686 isolated exploratory sensitivity display",
                "Pre-specified module": MODULE_LABELS[row["module_id"]],
                "Recorded lower vs higher group n": "2 vs 3",
                "Higher-minus-lower score difference": f"{row['mean_difference_target_minus_comparison']:.4f}",
                "Welch 95% interval": f"[{row['analytic_ci_lower']:.4f}, {row['analytic_ci_upper']:.4f}]",
                "Bootstrap 95% interval": f"[{row['bootstrap_ci_lower']:.4f}, {row['bootstrap_ci_upper']:.4f}]",
                "LOKO direction retention (fraction)": f"{row['leave_one_key_out_direction_retention_fraction']:.3f}",
                "Permanent exclusion": EXCLUDED_GSM,
                "Interpretation boundary": (
                    "Isolated descriptive sensitivity only; excluded from the default 20-effect summary, "
                    "sign-alignment count, and main figures/tables; not validation or replication."
                ),
            }
        )
    write_csv(
        output,
        rows,
        [
            "Dataset and analysis role",
            "Pre-specified module",
            "Recorded lower vs higher group n",
            "Higher-minus-lower score difference",
            "Welch 95% interval",
            "Bootstrap 95% interval",
            "LOKO direction retention (fraction)",
            "Permanent exclusion",
            "Interpretation boundary",
        ],
    )


def make_figure(effects: pd.DataFrame, output_dir: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(8.8, 4.2), gridspec_kw={"width_ratios": [1.35, 0.95]})
    y_positions = np.arange(len(MODULE_ORDER))[::-1]
    effect_values = effects["mean_difference_target_minus_comparison"].to_numpy(dtype=float)
    ci_lower = effects["analytic_ci_lower"].to_numpy(dtype=float)
    ci_upper = effects["analytic_ci_upper"].to_numpy(dtype=float)
    bound = max(1.0, float(np.ceil(np.abs(np.r_[ci_lower, ci_upper]).max() * 2.0) / 2.0))
    colors = np.where(effect_values >= 0, "#0072B2", "#D55E00")
    for ypos, value, lower, upper, color in zip(y_positions, effect_values, ci_lower, ci_upper, colors, strict=True):
        axes[0].errorbar(
            value,
            ypos,
            xerr=np.array([[value - lower], [upper - value]]),
            fmt="o",
            color=color,
            ecolor=color,
            markersize=5.5,
            elinewidth=1.35,
            capsize=3,
            zorder=3,
        )
        axes[0].text(
            bound * 0.98,
            ypos,
            f"{value:+.2f}",
            ha="right",
            va="center",
            fontsize=8,
        )
    axes[0].axvline(0, color="#4D4D4D", linewidth=0.8, zorder=0)
    axes[0].grid(axis="x", color="#D9D9D9", linewidth=0.5, zorder=0)
    axes[0].set_xlim(-bound, bound)
    axes[0].set_yticks(y_positions)
    axes[0].set_yticklabels([MODULE_LABELS[module] for module in MODULE_ORDER], fontsize=8)
    axes[0].set_xlabel("Severe minus mild module-score difference\n(descriptive Welch 95% interval)", fontsize=8)
    axes[0].set_title("A  Isolated GSE251686 effect display", loc="left", fontsize=10, fontweight="bold")

    retention = effects["leave_one_key_out_direction_retention_fraction"].to_numpy(dtype=float)
    for ypos, value, effect in zip(y_positions, retention, effect_values, strict=True):
        color = "#0072B2" if effect >= 0 else "#D55E00"
        axes[1].plot([0, value], [ypos, ypos], color="#BDBDBD", linewidth=1.0, zorder=1)
        axes[1].scatter(value, ypos, s=42, color=color, zorder=3)
        axes[1].text(min(value + 0.035, 1.03), ypos, f"{value:.2f}", va="center", fontsize=8)
    axes[1].axvline(0.80, color="#4D4D4D", linestyle="--", linewidth=0.8, zorder=0)
    axes[1].set_xlim(0, 1.08)
    axes[1].set_xticks([0, 0.5, 0.8, 1.0])
    axes[1].set_yticks(y_positions)
    axes[1].set_yticklabels([])
    axes[1].grid(axis="x", color="#D9D9D9", linewidth=0.5, zorder=0)
    axes[1].set_xlabel("Leave-one-key-out direction retention", fontsize=8)
    axes[1].set_title("B  Stability display", loc="left", fontsize=10, fontweight="bold")
    axes[1].text(0.80, -0.68, "0.80 display threshold", ha="center", va="top", fontsize=7, color="#4D4D4D")

    for axis in axes:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.tick_params(labelsize=8)
    figure.text(
        0.5,
        -0.055,
        "Mild n=2 and severe n=3 presumed sample/library keys after permanent exclusion of GSM7986002. "
        "This isolated exploratory sensitivity display is excluded from the default 20-effect summary and "
        "does not provide a p-value, validation, replication, pooled estimate, or causal interpretation.",
        ha="center",
        va="top",
        fontsize=7.2,
        wrap=True,
    )
    figure.tight_layout()
    figure.savefig(output_dir / "supplementary_figure_s1_gse251686_exploratory_sensitivity.pdf", bbox_inches="tight")
    figure.savefig(
        output_dir / "supplementary_figure_s1_gse251686_exploratory_sensitivity.png",
        dpi=600,
        bbox_inches="tight",
    )
    plt.close(figure)


def write_readme(output: Path) -> None:
    output.write_text(
        "# Supplementary IVDD Result Deliverables\n\n"
        "## Default descriptive analysis (unchanged)\n\n"
        "The four default score cohorts and their 20-effect descriptive summary remain the only authoritative "
        "main-analysis results. Supplementary files in this section document that analysis; they do not recalculate "
        "or extend it.\n\n"
        "- **S1:** isolated GSE251686 descriptive sensitivity. It remains outside the default summary; `GSM7986002` "
        "is permanently excluded after a Matrix Market stream-integrity failure.\n"
        "- **S2:** cohort disposition, observation keys, identity checks, and default-summary eligibility.\n"
        "- **S3:** locked program-module definitions and hashes.\n"
        "- **S4:** default leave-one-key-out stability results.\n"
        "- **S5a/S5b:** discovery retained-cell-threshold sensitivity.\n"
        "- **S6:** reproducibility contract index.\n\n"
        "S1--S6 do not provide confirmatory inference, replication, causal claims, biomarker claims, or therapeutic "
        "interpretation. The default summary contains no p-values or hypothesis tests.\n\n"
        "## Supplementary Table S7: exploratory four-cohort NP meta-analysis\n\n"
        "S7a--S7d are a separate cohort-level random-effects synthesis of the four default NP cohorts. They are "
        "exploratory supplementary analyses, not a replacement for the default descriptive analysis and not confirmation, "
        "replication, mechanism, biomarker, or therapy evidence. The four cohorts retain `confirmatory_eligible=false`; "
        "GSE251686 and `GSM7986002` remain outside S7.\n\n"
        "## Supplementary Table S8: post-hoc six-cohort NP expansion\n\n"
        "S8a--S8d are a separate post-hoc exploratory expansion of S7 that adds GSE186542 and the GSE167931 FPKM "
        "representation, for six cohort/module contrasts per module. The paired GSE167931 TPM representation is a "
        "same-sample processing sensitivity, not an additional cohort. S8 does not alter the unchanged 20-effect default "
        "analysis and is not independent patient-level validation.\n\n"
        "## Supplementary Table S9: source-family replacement sensitivity\n\n"
        "S9a--S9d replace the GSE167931 FPKM representation in S8 with the native GSE245147 Degenerated n=3 versus "
        "No-degenerated n=3 comparison. GSE167931 and GSE245147 are never pooled together because their source "
        "lab/author family overlaps and patient-level reuse cannot be excluded. S9 is a post-hoc source-family sensitivity "
        "only; it does not add a seventh independent cohort, alter the frozen default summary, or establish validation.\n\n"
        "S7-S9 report HKSJ and four-module BH p-values only for transparent description; none can establish confirmation, "
        "replication, biomarker, mechanistic, causal, or therapeutic evidence.\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-dir", type=Path, required=True)
    parser.add_argument("--table-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    args = parser.parse_args()

    package_dir = args.package_dir.resolve()
    table_dir = args.table_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    if table_dir.parent != figure_dir.parent:
        raise ValueError("Supplementary table and figure directories must share one results root")
    results_root = table_dir.parent
    manifest = require_isolated_package(package_dir)
    effects = load_effects(package_dir)
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    table_path = table_dir / "supplementary_table_s1_gse251686_exploratory_effects.csv"
    figure_pdf = figure_dir / "supplementary_figure_s1_gse251686_exploratory_sensitivity.pdf"
    figure_png = figure_dir / "supplementary_figure_s1_gse251686_exploratory_sensitivity.png"
    readme_path = table_dir / "README.md"
    artifact_path = table_dir / "run_artifacts.csv"
    manifest_path = table_dir / "run_manifest.json"

    write_supplementary_table(effects, table_path)
    make_figure(effects, figure_dir)
    write_readme(readme_path)
    artifact_names = [table_path.name, readme_path.name, artifact_path.name, manifest_path.name]
    write_csv(
        artifact_path,
        [
            {
                "artifact": name,
                "generated_by_this_invocation": "true",
                "purpose": "Isolated GSE251686 sensitivity artifact; not a default-summary input.",
            }
            for name in artifact_names
        ],
        ["artifact", "generated_by_this_invocation", "purpose"],
    )
    output_manifest = {
        "schema_version": 1,
        "source_package": str(package_dir),
        "source_package_run_manifest_sha256": sha256(package_dir / "run_manifest.json"),
        "source_default_cross_cohort_summary_inclusion": manifest["default_cross_cohort_summary_inclusion"],
        "source_permanent_exclusion": manifest["permanent_exclusion"],
        "default_cross_cohort_summary_inclusion": False,
        "main_figure_or_table_inclusion": False,
        "generated_artifact_path_root": str(results_root),
        "unit_of_observation": "presumed sample/library key; cells nested",
        "no_hypothesis_tests_or_p_values": True,
        "generated_artifact_sha256": {
            str(table_path.relative_to(results_root)): sha256(table_path),
            str(readme_path.relative_to(results_root)): sha256(readme_path),
            str(figure_pdf.relative_to(results_root)): sha256(figure_pdf),
            str(figure_png.relative_to(results_root)): sha256(figure_png),
        },
        "self_referential_hash_exclusions": [artifact_path.name, manifest_path.name],
    }
    manifest_path.write_text(json.dumps(output_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
