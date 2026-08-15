"""Create tables and figures for the audited descriptive IVDD score summary.

This script deliberately visualizes cohort-specific donor or presumed sample-key
effects. It never pools cohorts, computes p-values, or labels sign alignment as
replication. Inputs are limited to the current invocation artifacts listed by
the donor-module summary directory.
"""

from __future__ import annotations

import argparse
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
COHORT_ORDER = [
    "GSE230809_discovery",
    "GSE244889_directional",
    "GSE153066_support",
    "GSE165722_score_level",
]
COHORT_LABELS = {
    "GSE230809_discovery": "GSE230809 NP\nexploratory (3 vs 8)",
    "GSE244889_directional": "GSE244889 NP\ndirectional support (4 vs 3)",
    "GSE153066_support": "GSE153066 NP\ncount-level support (8 vs 8)",
    "GSE165722_score_level": "GSE165722 NP\nscore-level support (4 vs 4)",
}
# Labels used only in the formal manuscript figures. They describe the source
# and score representation without implying validation, replication, or support.
FORMAL_COHORT_LABELS = {
    "GSE230809_discovery": "GSE230809 NP\ndiscovery cohort (3 vs 8)",
    "GSE244889_directional": "GSE244889 NP\nexternal score-level cohort (4 vs 3)",
    "GSE153066_support": "GSE153066 NP\nexternal dense-count cohort (8 vs 8)",
    "GSE165722_score_level": "GSE165722 NP\nexternal normalized-count cohort (4 vs 4)",
}
FORMAL_DIRECTION_LABELS = {
    "GSE230809_discovery": "GSE230809\ndiscovery cohort\n3 vs 8",
    "GSE244889_directional": "GSE244889\nexternal\nscore-level cohort\n4 vs 3",
    "GSE153066_support": "GSE153066\nexternal\ndense-count cohort\n8 vs 8",
    "GSE165722_score_level": "GSE165722\nexternal\nnormalized-count cohort\n4 vs 4",
}
FORMAL_MODULE_LABELS = {
    **MODULE_LABELS,
    "inflammatory_nfkb": "Inflammatory / NF-κB",
}
# Table 1 describes the parent project's complete AF/NP score-to-ledger
# crosswalk, so its label must not imply that the 24 matched keys are NP-only.
TABLE_ONE_COHORT_LABELS = {
    **COHORT_LABELS,
    "GSE230809_discovery": "GSE230809 AF + NP\nexploratory (AF 3 vs 10; NP 3 vs 8)",
}
COHORT_COLORS = {
    "GSE230809_discovery": "#0072B2",
    "GSE244889_directional": "#D55E00",
    "GSE153066_support": "#009E73",
    "GSE165722_score_level": "#CC79A7",
}


def require_current_artifacts(summary_dir: Path) -> None:
    artifacts_path = summary_dir / "run_artifacts.csv"
    if not artifacts_path.is_file():
        raise FileNotFoundError(f"Missing current-run artifact ledger: {artifacts_path}")
    artifacts = pd.read_csv(artifacts_path, dtype=str)
    current = set(artifacts.loc[artifacts["generated_by_this_invocation"].eq("true"), "artifact"])
    required = {
        "donor_module_effects.csv",
        "cross_cohort_direction_consistency.csv",
        "cohort_score_availability.csv",
        "score_ledger_identity_crosswalk.csv",
    }
    missing = required - current
    if missing:
        raise ValueError(f"Current run artifact ledger is missing required entries: {sorted(missing)}")
    missing_files = [name for name in required if not (summary_dir / name).is_file()]
    if missing_files:
        raise FileNotFoundError(f"Current-run files are missing: {missing_files}")


def clean_effects(summary_dir: Path) -> pd.DataFrame:
    effects = pd.read_csv(summary_dir / "donor_module_effects.csv", dtype=str)
    needed = {
        "cohort_id",
        "compartment",
        "module_id",
        "comparison_n",
        "target_n",
        "mean_difference_target_minus_comparison",
        "analytic_ci_lower",
        "analytic_ci_upper",
        "bootstrap_ci_lower",
        "bootstrap_ci_upper",
        "effect_direction",
        "leave_one_out_direction_retention_fraction",
        "confirmatory_eligible",
        "cohort_limitations",
    }
    missing = needed - set(effects.columns)
    if missing:
        raise ValueError(f"Effect table is missing required columns: {sorted(missing)}")
    effects = effects.loc[
        effects["cohort_id"].isin(COHORT_ORDER) & effects["compartment"].eq("NP")
    ].copy()
    for column in [
        "comparison_n",
        "target_n",
        "mean_difference_target_minus_comparison",
        "analytic_ci_lower",
        "analytic_ci_upper",
        "bootstrap_ci_lower",
        "bootstrap_ci_upper",
        "leave_one_out_direction_retention_fraction",
    ]:
        effects[column] = pd.to_numeric(effects[column], errors="coerce")
    expected = {(cohort, module) for cohort in COHORT_ORDER for module in MODULE_ORDER}
    observed = set(zip(effects["cohort_id"], effects["module_id"], strict=True))
    if observed != expected:
        raise ValueError(
            "NP effect table does not contain exactly the expected current cohort/module grid; "
            f"missing={sorted(expected - observed)}, extra={sorted(observed - expected)}"
        )
    if not effects["confirmatory_eligible"].eq("false").all():
        raise ValueError("This deliverable script is restricted to the current non-confirmatory summary.")
    return effects


def write_table_one(summary_dir: Path, output: Path) -> None:
    availability = pd.read_csv(summary_dir / "cohort_score_availability.csv", dtype=str)
    crosswalk = pd.read_csv(summary_dir / "score_ledger_identity_crosswalk.csv", dtype=str)
    availability = availability.loc[availability["cohort_id"].isin(COHORT_ORDER)].copy()
    if set(availability["cohort_id"]) != set(COHORT_ORDER):
        raise ValueError("Cohort availability table lacks a current NP score cohort.")
    crosswalk = crosswalk.loc[crosswalk["cohort_id"].isin(COHORT_ORDER)].copy()
    if not crosswalk["crosswalk_status"].eq("pass_exact_identity_crosswalk").all():
        raise ValueError("A score-to-ledger identity crosswalk did not pass exact matching.")
    counts = crosswalk.set_index("cohort_id")[
        ["unique_retained_score_sample_keys", "matched_keys", "crosswalk_status"]
    ].rename(
        columns={
            "unique_retained_score_sample_keys": "sample_keys_after_score_filter",
            "matched_keys": "exactly_matched_ledger_keys",
            "crosswalk_status": "identity_crosswalk",
        }
    )
    table = availability.set_index("cohort_id").join(counts, how="left").reset_index()
    table["cohort_display"] = table["cohort_id"].map(TABLE_ONE_COHORT_LABELS)
    table["confirmatory_status"] = "No current cohort is confirmatory eligible"
    table["interpretation"] = table["analysis_boundary"]
    table = table[
        [
            "cohort_display",
            "audited_group_structure",
            "sample_key_type",
            "sample_keys_after_score_filter",
            "exactly_matched_ledger_keys",
            "identity_crosswalk",
            "interpretation",
            "confirmatory_status",
        ]
    ]
    table.columns = [
        "Cohort and role",
        "Recorded group structure",
        "Observation key",
        "Keys after score filter",
        "Exact ledger matches",
        "Identity check",
        "Interpretation boundary",
        "Confirmatory status",
    ]
    table.to_csv(output, index=False, encoding="utf-8")


def write_table_two(effects: pd.DataFrame, output: Path) -> None:
    table = effects.copy()
    table["_cohort_order"] = pd.Categorical(table["cohort_id"], categories=COHORT_ORDER, ordered=True)
    table["_module_order"] = pd.Categorical(table["module_id"], categories=MODULE_ORDER, ordered=True)
    table = table.sort_values(["_cohort_order", "_module_order"], kind="stable").reset_index(drop=True)
    table["cohort_display"] = table["cohort_id"].map(COHORT_LABELS)
    table["module_display"] = table["module_id"].map(MODULE_LABELS)
    table["recorded_group_n"] = (
        table["comparison_n"].astype("Int64").astype(str)
        + " vs "
        + table["target_n"].astype("Int64").astype(str)
    )
    table["mean_difference"] = table["mean_difference_target_minus_comparison"].map("{:.4f}".format)
    table["welch_95_ci"] = table.apply(
        lambda row: f"[{row['analytic_ci_lower']:.4f}, {row['analytic_ci_upper']:.4f}]", axis=1
    )
    table["bootstrap_95_ci"] = table.apply(
        lambda row: f"[{row['bootstrap_ci_lower']:.4f}, {row['bootstrap_ci_upper']:.4f}]", axis=1
    )
    # Keep three decimals so proportions such as 0.625 are not silently
    # rounded to a materially different two-decimal display (0.62).
    table["lodo_direction_retention"] = table["leave_one_out_direction_retention_fraction"].map(
        "{:.3f}".format
    )
    table["direction"] = np.where(
        table["mean_difference_target_minus_comparison"] > 0,
        "positive higher recorded severity",
        "negative lower recorded severity",
    )
    table = table[
        [
            "cohort_display",
            "module_display",
            "recorded_group_n",
            "mean_difference",
            "welch_95_ci",
            "bootstrap_95_ci",
            "direction",
            "lodo_direction_retention",
            "cohort_limitations",
        ]
    ]
    table.columns = [
        "Cohort and role",
        "Pre-specified module",
        "Lower vs higher group n",
        "Higher-minus-lower score difference",
        "Welch 95% interval",
        "Bootstrap 95% interval",
        "Observed direction",
        "LODO direction retention (fraction)",
        "Cohort limitation",
    ]
    table.to_csv(output, index=False, encoding="utf-8")


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.axvline(0, color="#4D4D4D", linewidth=0.8, zorder=0)
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.5, zorder=0)
    ax.tick_params(axis="both", labelsize=7)


def make_effect_figure(
    effects: pd.DataFrame,
    output_dir: Path,
    cohort_labels: dict[str, str] = COHORT_LABELS,
    module_labels: dict[str, str] = MODULE_LABELS,
    formal: bool = False,
) -> None:
    # Keep Figure 1 compact enough for its caption and image to stay on one
    # landscape manuscript page at the fixed 9.2-inch insertion width.
    figure_size = (9.4, 5.7) if formal else (7.2, 5.6)
    fig, axes = plt.subplots(2, 2, figsize=figure_size, sharex=True, constrained_layout=True)
    all_boundaries = pd.concat(
        [effects["analytic_ci_lower"], effects["analytic_ci_upper"]], ignore_index=True
    ).dropna()
    bound = max(1.0, float(np.ceil(all_boundaries.abs().max() * 2.0) / 2.0))
    xlim = (-bound, bound)
    y_positions = np.arange(len(COHORT_ORDER))[::-1]
    for panel, (axis, module_id) in enumerate(zip(axes.flat, MODULE_ORDER, strict=True)):
        sub = effects.loc[effects["module_id"].eq(module_id)].set_index("cohort_id").loc[COHORT_ORDER]
        estimate = sub["mean_difference_target_minus_comparison"].to_numpy(dtype=float)
        lower = sub["analytic_ci_lower"].to_numpy(dtype=float)
        upper = sub["analytic_ci_upper"].to_numpy(dtype=float)
        for ypos, cohort, value, lo, hi in zip(y_positions, COHORT_ORDER, estimate, lower, upper, strict=True):
            color = COHORT_COLORS[cohort]
            axis.errorbar(
                value,
                ypos,
                xerr=np.array([[value - lo], [hi - value]]),
                fmt="o",
                markersize=5,
                color=color,
                ecolor=color,
                elinewidth=1.25,
                capsize=2.5,
                zorder=3,
            )
        axis.set_title(module_labels[module_id], fontsize=9, fontweight="bold", loc="left")
        axis.text(-0.13, 1.07, chr(65 + panel), transform=axis.transAxes, fontsize=10, fontweight="bold")
        axis.set_yticks(y_positions)
        axis.set_yticklabels([cohort_labels[cohort] for cohort in COHORT_ORDER], fontsize=6.4)
        axis.set_xlim(xlim)
        style_axis(axis)
    if formal:
        fig.supxlabel(
            "Higher-minus-lower module-score difference (unitless; cohort-specific; not pooled)",
            fontsize=8,
        )
    else:
        for axis in axes[1, :]:
            axis.set_xlabel(
                "Higher recorded severity minus lower recorded severity\n"
                "unitless module-score difference (cohort-specific; not pooled)",
                fontsize=8,
            )
    if not formal:
        fig.text(
            0.5,
            -0.035,
            "Points are unweighted donor or presumed sample-key mean differences; bars are descriptive Welch 95% intervals. "
            "No p-values, pooled estimate, or replication decision is shown. Magnitudes are not comparable across processing scales; GSE165722 uses supplied normalized-count values, not raw CPM.",
            ha="center",
            va="top",
            fontsize=7,
        )
    save_kwargs = {"bbox_inches": "tight", "pad_inches": 0.14} if formal else {"bbox_inches": "tight"}
    fig.savefig(output_dir / "figure_1_np_cohort_module_effects.pdf", **save_kwargs)
    fig.savefig(output_dir / "figure_1_np_cohort_module_effects.png", dpi=600, **save_kwargs)
    plt.close(fig)


def make_direction_figure(
    summary_dir: Path,
    effects: pd.DataFrame,
    output_dir: Path,
    cohort_labels: dict[str, str] | None = None,
    module_labels: dict[str, str] = MODULE_LABELS,
    include_panel_label: bool = True,
    formal: bool = False,
) -> None:
    directions = pd.read_csv(summary_dir / "cross_cohort_direction_consistency.csv", dtype=str)
    directions = directions.loc[
        directions["compartment"].eq("NP") & directions["module_id"].isin(MODULE_ORDER)
    ].copy()
    directions = directions.set_index("module_id").loc[MODULE_ORDER].reset_index()
    if len(directions) != len(MODULE_ORDER):
        raise ValueError("Direction-consistency table lacks a pre-specified NP module.")
    matrix = effects.pivot(index="module_id", columns="cohort_id", values="mean_difference_target_minus_comparison")
    matrix = matrix.loc[MODULE_ORDER, COHORT_ORDER]
    signs = np.sign(matrix.to_numpy(dtype=float))
    figure_size = (10.4, 4.8) if formal else (8.6, 3.45)
    fig, ax = plt.subplots(figsize=figure_size, constrained_layout=True)
    cmap = matplotlib.colors.ListedColormap(["#D55E00", "#F7F7F7", "#0072B2"])
    image = ax.imshow(signs, vmin=-1, vmax=1, cmap=cmap, aspect="auto")
    ax.set_xticks(np.arange(len(COHORT_ORDER)))
    if cohort_labels is None:
        cohort_labels = {
            "GSE230809_discovery": "GSE230809\nexploratory\n3 vs 8",
            "GSE244889_directional": "GSE244889\ndirectional support\n4 vs 3",
            "GSE153066_support": "GSE153066\ncount-level support\n8 vs 8",
            "GSE165722_score_level": "GSE165722\nscore-level support\n4 vs 4",
        }
    ax.set_xticklabels([cohort_labels[cohort] for cohort in COHORT_ORDER], fontsize=6.4 if formal else 6.7)
    if formal:
        ax.tick_params(axis="x", pad=4)
    ax.set_yticks(np.arange(len(MODULE_ORDER)))
    ax.set_yticklabels([module_labels[item] for item in MODULE_ORDER], fontsize=8)
    for row_idx, module_id in enumerate(MODULE_ORDER):
        direction_row = directions.loc[directions["module_id"].eq(module_id)].iloc[0]
        for col_idx, cohort_id in enumerate(COHORT_ORDER):
            value = float(matrix.loc[module_id, cohort_id])
            ax.text(col_idx, row_idx, f"{value:+.2f}", ha="center", va="center", fontsize=8, color="#111111")
    ax.set_title("NP cohort-specific directions and descriptive sign alignment", fontsize=10, fontweight="bold", loc="left")
    if include_panel_label:
        ax.text(-0.07, 1.06, "A", transform=ax.transAxes, fontsize=10, fontweight="bold")
    ax.set_xticks(np.arange(-0.5, len(COHORT_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(MODULE_ORDER), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)
    colorbar = fig.colorbar(image, ax=ax, pad=0.02, shrink=0.78, ticks=[-1, 0, 1])
    colorbar.ax.set_yticklabels(["negative", "zero", "positive"])
    colorbar.ax.tick_params(labelsize=7)
    if not formal:
        fig.text(
            0.5,
            -0.06,
            "Values are cohort-specific, unitless higher-minus-lower score differences; colors encode sign only. "
            "Hypoxia/oxidative stress: 4 positive and 0 negative cohort directions. Alignment is descriptive and does not establish replication or a universal program.",
            ha="center",
            va="top",
            fontsize=7,
        )
    save_kwargs = {"bbox_inches": "tight", "pad_inches": 0.14} if formal else {"bbox_inches": "tight"}
    fig.savefig(output_dir / "figure_2_np_direction_alignment.pdf", **save_kwargs)
    fig.savefig(output_dir / "figure_2_np_direction_alignment.png", dpi=600, **save_kwargs)
    plt.close(fig)


def write_readme(output: Path) -> None:
    text = """# Current Descriptive IVDD Deliverables

Generated by `scripts/make_current_summary_deliverables.py` from the files
listed in `data/derived/donor_module_effect_summary/run_artifacts.csv`.

- `table_1_current_cohort_roles.csv`: cohort structure, identity checks, and boundaries.
- `table_2_np_module_effects.csv`: all 16 current NP module effects, descriptive intervals, and LODO retention.
- `figure_1_np_cohort_module_effects.pdf/png`: cohort-specific effect display using Welch 95% intervals.
- `figure_2_np_direction_alignment.pdf/png`: sign and effect matrix without pooling.

All results are donor or presumed donor/sample-key-level descriptive summaries.
Cells are nested observations. Score magnitudes are cohort-specific and are not
pooled across processing scales. No p-values, formal meta-analysis, replication
adjudication, causal inference, biomarker claim, or therapeutic interpretation
is included.
"""
    output.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-dir", type=Path, required=True)
    parser.add_argument("--table-dir", type=Path)
    parser.add_argument("--figure-dir", type=Path)
    parser.add_argument(
        "--formal-figure-dir",
        type=Path,
        help=(
            "Optional separate output directory for manuscript figures. The formal versions "
            "use neutral cohort labels and omit the single-panel Figure 2 label."
        ),
    )
    parser.add_argument(
        "--formal-only",
        action="store_true",
        help="Generate only the separate formal manuscript figures; requires --formal-figure-dir.",
    )
    args = parser.parse_args()
    summary_dir = args.summary_dir.resolve()
    require_current_artifacts(summary_dir)
    effects = clean_effects(summary_dir)
    if args.formal_only:
        if args.formal_figure_dir is None:
            parser.error("--formal-only requires --formal-figure-dir")
        formal_figure_dir = args.formal_figure_dir.resolve()
        formal_figure_dir.mkdir(parents=True, exist_ok=True)
        make_effect_figure(
            effects,
            formal_figure_dir,
            cohort_labels=FORMAL_COHORT_LABELS,
            module_labels=FORMAL_MODULE_LABELS,
            formal=True,
        )
        make_direction_figure(
            summary_dir,
            effects,
            formal_figure_dir,
            cohort_labels=FORMAL_DIRECTION_LABELS,
            module_labels=FORMAL_MODULE_LABELS,
            include_panel_label=False,
            formal=True,
        )
        return 0
    if args.table_dir is None or args.figure_dir is None:
        parser.error("--table-dir and --figure-dir are required unless --formal-only is used")
    table_dir = args.table_dir.resolve()
    figure_dir = args.figure_dir.resolve()
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    write_table_one(summary_dir, table_dir / "table_1_current_cohort_roles.csv")
    write_table_two(effects, table_dir / "table_2_np_module_effects.csv")
    make_effect_figure(effects, figure_dir)
    make_direction_figure(summary_dir, effects, figure_dir)
    if args.formal_figure_dir is not None:
        formal_figure_dir = args.formal_figure_dir.resolve()
        formal_figure_dir.mkdir(parents=True, exist_ok=True)
        make_effect_figure(
            effects,
            formal_figure_dir,
            cohort_labels=FORMAL_COHORT_LABELS,
            module_labels=FORMAL_MODULE_LABELS,
            formal=True,
        )
        make_direction_figure(
            summary_dir,
            effects,
            formal_figure_dir,
            cohort_labels=FORMAL_DIRECTION_LABELS,
            module_labels=FORMAL_MODULE_LABELS,
            include_panel_label=False,
            formal=True,
        )
    write_readme(table_dir / "README.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
