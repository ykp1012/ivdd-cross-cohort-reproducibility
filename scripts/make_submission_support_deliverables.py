"""Create auditable supplementary tables, a workflow figure, and a graphical abstract.

The deliverables are derived from the frozen default descriptive IVDD summary,
the locked module ledger, the discovery retained-cell sensitivity package, the
isolated GSE251686 exploratory package, and three separately labeled NP
random-effects packages (S7-S9). The script rejects any attempt to turn
GSE251686 into a default-summary input or to promote any current contrast to
confirmatory inference.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import pandas as pd


DEFAULT_COHORT_ORDER = [
    "GSE230809_discovery",
    "GSE244889_directional",
    "GSE153066_support",
    "GSE165722_score_level",
]
DEFAULT_REQUIRED_ARTIFACTS = {
    "cohort_score_availability.csv",
    "cross_cohort_direction_consistency.csv",
    "donor_module_effects.csv",
    "donor_module_leave_one_out.csv",
    "run_parameters.csv",
    "score_ledger_identity_crosswalk.csv",
}
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
GSE251686_ID = "GSE251686_exploratory"
GSE251686_BAD_GSM = "GSM7986002"
META_ANALYSIS_SPECS = {
    "s7_exploratory_four_cohort": {
        "directory": "np_exploratory_meta_analysis",
        "expected_k": 4,
        "expected_cohorts": set(DEFAULT_COHORT_ORDER),
        "primary_table": "results/supplementary_tables/supplementary_table_s7b_np_meta_analysis_primary_results.csv",
        "study_table": "results/supplementary_tables/supplementary_table_s7a_np_meta_analysis_study_effects.csv",
        "purpose": "Separate four-cohort exploratory SMDH random-effects synthesis; non-confirmatory.",
    },
    "s8_post_hoc_external_expansion": {
        "directory": "np_post_hoc_external_expansion_meta_analysis",
        "expected_k": 6,
        "expected_cohorts": {
            *DEFAULT_COHORT_ORDER,
            "GSE186542_external_count_support",
            "GSE167931_external_fpkm_support",
        },
        "primary_table": "results/supplementary_tables/supplementary_table_s8b_np_post_hoc_external_expansion_primary_results.csv",
        "study_table": "results/supplementary_tables/supplementary_table_s8a_np_post_hoc_external_expansion_study_effects.csv",
        "purpose": "Separate post hoc six-cohort SMDH random-effects expansion; non-confirmatory.",
    },
    "s9_source_family_replacement": {
        "directory": "np_source_family_replacement_meta_analysis",
        "expected_k": 6,
        "expected_cohorts": {
            *DEFAULT_COHORT_ORDER,
            "GSE186542_external_count_support",
            "GSE245147_external_native_comparison_support",
        },
        "primary_table": "results/supplementary_tables/supplementary_table_s9b_np_source_family_replacement_sensitivity_primary_results.csv",
        "study_table": "results/supplementary_tables/supplementary_table_s9a_np_source_family_replacement_sensitivity_study_effects.csv",
        "purpose": "Separate source-family replacement SMDH random-effects sensitivity; non-confirmatory.",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Missing expected file: {path}")
    return pd.read_csv(path, dtype=str)


def verify_manifest_artifacts(directory: Path, manifest: dict[str, Any]) -> None:
    for name, expected_hash in manifest.get("generated_artifact_sha256", {}).items():
        path = directory / name
        if not path.is_file():
            raise FileNotFoundError(f"Manifest-listed artifact is missing: {path}")
        if sha256(path) != expected_hash:
            raise ValueError(f"Manifest hash mismatch: {path.name}")


def project_relative_path(project_root: Path, relative: str) -> Path:
    """Resolve a manifest path while rejecting paths outside the project root."""
    path = (project_root / relative).resolve()
    try:
        path.relative_to(project_root)
    except ValueError as exc:
        raise ValueError(f"Manifest path escapes the project root: {relative}") from exc
    return path


def require_meta_analysis_package(
    project_root: Path,
    package_dir: Path,
    spec: dict[str, Any],
) -> dict[str, Any]:
    """Validate an isolated S7/S8/S9 package before it is indexed for submission."""
    manifest_path = package_dir / "meta_analysis_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing meta-analysis manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    scope = manifest.get("scope", {})
    if set(scope.get("cohorts", [])) != spec["expected_cohorts"]:
        raise ValueError(f"{package_dir.name} has an unexpected cohort set")
    if set(scope.get("modules", [])) != set(MODULE_ORDER):
        raise ValueError(f"{package_dir.name} has an unexpected module set")
    if "SMDH" not in scope.get("primary_analysis", "") or "REML" not in scope.get("primary_analysis", ""):
        raise ValueError(f"{package_dir.name} does not retain the SMDH/REML primary-analysis contract")
    if "Knapp-Hartung" not in scope.get("primary_analysis", ""):
        raise ValueError(f"{package_dir.name} does not retain the Knapp-Hartung interval contract")

    generated = manifest.get("generated_artifact_sha256", {})
    if not generated:
        raise ValueError(f"{package_dir.name} has no generated artifact hashes")
    artifact_paths: list[Path] = []
    for relative, expected_hash in generated.items():
        path = project_relative_path(project_root, relative)
        if not path.is_file():
            raise FileNotFoundError(f"Meta-analysis artifact is missing: {path}")
        if sha256(path) != expected_hash:
            raise ValueError(f"Meta-analysis artifact hash mismatch: {path}")
        artifact_paths.append(path)

    primary_table = project_relative_path(project_root, spec["primary_table"])
    study_table = project_relative_path(project_root, spec["study_table"])
    if primary_table not in artifact_paths or study_table not in artifact_paths:
        raise ValueError(f"{package_dir.name} manifest does not hash both required result tables")
    primary = read_csv(primary_table)
    if len(primary) != len(MODULE_ORDER) or set(primary["module_id"]) != set(MODULE_ORDER):
        raise ValueError(f"{package_dir.name} primary table must contain the four locked modules")
    if not primary["k"].astype(str).eq(str(spec["expected_k"])).all():
        raise ValueError(f"{package_dir.name} primary table has an unexpected cohort count")
    if not primary["effect_measure"].str.contains("SMDH", regex=False).all():
        raise ValueError(f"{package_dir.name} primary table has an unexpected effect measure")
    if not primary["tau_squared_method"].eq("REML").all() or not primary["confidence_interval_method"].eq("Knapp-Hartung").all():
        raise ValueError(f"{package_dir.name} primary table has an unexpected model contract")
    if not primary["fit_control_note"].str.contains("maxiter=10000", regex=False).all():
        raise ValueError(f"{package_dir.name} primary table does not record the REML iteration control")
    required_numeric_columns = [
        "pooled_standardized_mean_difference",
        "ci_lower",
        "ci_upper",
        "hksj_p_value",
        "hksj_p_value_BH_four_modules",
    ]
    for column in required_numeric_columns:
        values = pd.to_numeric(primary[column], errors="coerce")
        if not values.notna().all():
            raise ValueError(f"{package_dir.name} primary table has non-finite {column} values")

    study = read_csv(study_table)
    if set(study["cohort_id"]) != spec["expected_cohorts"]:
        raise ValueError(f"{package_dir.name} study table has an unexpected cohort set")
    if not study["all_cohorts_confirmatory_eligible"].astype(str).str.lower().eq("false").all():
        raise ValueError(f"{package_dir.name} contains a confirmatory-eligible study effect")
    return {
        "manifest_path": manifest_path,
        "manifest": manifest,
        "artifact_paths": artifact_paths,
        "primary_table": primary_table,
        "study_table": study_table,
    }


def require_default_summary(summary_dir: Path) -> dict[str, Any]:
    manifest_path = summary_dir / "run_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing default-summary manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("result_contract") != "default_current_project_cross_cohort_summary_only":
        raise ValueError("Unexpected default-summary result contract")
    if manifest.get("default_summary_effect_count") != 20:
        raise ValueError("Default summary must contain exactly 20 effects")
    if manifest.get("all_effects_confirmatory_eligible") is not False:
        raise ValueError("Default summary may not contain confirmatory-eligible effects")
    if manifest.get("excluded_separate_exploratory_package", {}).get("cohort_id") != GSE251686_ID:
        raise ValueError("Default summary does not explicitly exclude the GSE251686 package")
    verify_manifest_artifacts(summary_dir, manifest)

    artifacts = read_csv(summary_dir / "run_artifacts.csv")
    present = set(artifacts.loc[artifacts["generated_by_this_invocation"].eq("true"), "artifact"])
    missing = DEFAULT_REQUIRED_ARTIFACTS - present
    if missing:
        raise ValueError(f"Default current-run ledger lacks required artifacts: {sorted(missing)}")

    effects = read_csv(summary_dir / "donor_module_effects.csv")
    if len(effects) != 20 or not effects["confirmatory_eligible"].eq("false").all():
        raise ValueError("Default effects do not satisfy the non-confirmatory 20-effect contract")
    np_effects = effects.loc[effects["compartment"].eq("NP")].copy()
    if len(np_effects) != 16:
        raise ValueError("Default summary must contain 16 NP effects")
    if set(np_effects["cohort_id"]) != set(DEFAULT_COHORT_ORDER):
        raise ValueError("Default NP effects have an unexpected cohort set")
    if set(np_effects["module_id"]) != set(MODULE_ORDER):
        raise ValueError("Default NP effects have an unexpected locked module set")

    availability = read_csv(summary_dir / "cohort_score_availability.csv")
    gse_row = availability.loc[availability["cohort_id"].eq(GSE251686_ID)]
    if len(gse_row) != 1 or gse_row.iloc[0]["effect_summary_included"] != "false":
        raise ValueError("GSE251686 must remain unavailable to the default effect summary")

    crosswalk = read_csv(summary_dir / "score_ledger_identity_crosswalk.csv")
    if set(crosswalk["cohort_id"]) != set(DEFAULT_COHORT_ORDER):
        raise ValueError("Default score-to-ledger identity crosswalk has an unexpected cohort set")
    if not crosswalk["crosswalk_status"].eq("pass_exact_identity_crosswalk").all():
        raise ValueError("A default score-to-ledger identity crosswalk did not pass")
    if int(pd.to_numeric(crosswalk["matched_keys"]).sum()) != 55:
        raise ValueError("Default score-to-ledger identity contract is not 55 matched keys")

    direction = read_csv(summary_dir / "cross_cohort_direction_consistency.csv")
    hypoxia = direction.loc[
        direction["compartment"].eq("NP") & direction["module_id"].eq("hypoxia_oxidative_stress")
    ]
    if len(hypoxia) != 1 or hypoxia.iloc[0]["n_positive"] != "4" or hypoxia.iloc[0]["n_negative"] != "0":
        raise ValueError("Default hypoxia sign-alignment record does not match the audited four-cohort display")
    return {
        "manifest": manifest,
        "effects": effects,
        "np_effects": np_effects,
        "availability": availability,
        "crosswalk": crosswalk,
        "direction": direction,
    }


def require_gse251686_package(package_dir: Path) -> dict[str, Any]:
    manifest_path = package_dir / "run_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing GSE251686 package manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("default_cross_cohort_summary_inclusion") is not False:
        raise ValueError("GSE251686 package is not isolated from the default summary")
    if manifest.get("permanent_exclusion") != GSE251686_BAD_GSM:
        raise ValueError("GSE251686 package does not retain the permanent malformed-GSM exclusion")
    if manifest.get("no_hypothesis_tests_or_p_values") is not True:
        raise ValueError("GSE251686 package has an unexpected inferential contract")
    if manifest.get("unit_of_observation") != "presumed sample/library key; cells nested":
        raise ValueError("GSE251686 package has an unexpected experimental unit")
    verify_manifest_artifacts(package_dir, manifest)
    effects = read_csv(package_dir / "GSE251686_exploratory_module_effects.csv")
    ledger = read_csv(package_dir / "GSE251686_exploratory_library_ledger.csv")
    if len(effects) != 4 or set(effects["module_id"]) != set(MODULE_ORDER):
        raise ValueError("GSE251686 package must contain four locked module effects")
    if not effects["comparison_n"].eq("2").all() or not effects["target_n"].eq("3").all():
        raise ValueError("GSE251686 package must retain mild n=2 and severe n=3")
    if not effects["excluded_record"].eq(GSE251686_BAD_GSM).all():
        raise ValueError("GSE251686 effect package does not retain the permanent exclusion")
    if not effects["confirmatory_eligible"].eq("false").all():
        raise ValueError("GSE251686 package contains an invalid confirmatory status")
    if len(ledger) != 5 or GSE251686_BAD_GSM in set(ledger["gsm"]):
        raise ValueError("GSE251686 library ledger does not contain the expected five included keys")
    for column in [
        "source_restricted_threshold_20_pass",
        "source_restricted_threshold_30_pass",
        "source_restricted_threshold_50_pass",
    ]:
        if not ledger[column].astype(str).str.lower().eq("true").all():
            raise ValueError(f"GSE251686 package does not pass its {column} gate")
    return {"manifest": manifest, "effects": effects, "ledger": ledger}


def require_gse251686_supplement(table_dir: Path, results_root: Path) -> dict[str, Any]:
    manifest_path = table_dir / "run_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing GSE251686 supplementary manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("default_cross_cohort_summary_inclusion") is not False:
        raise ValueError("GSE251686 supplementary outputs must be excluded from the default summary")
    if manifest.get("main_figure_or_table_inclusion") is not False:
        raise ValueError("GSE251686 supplementary outputs must not be main-figure or main-table inputs")
    if manifest.get("source_permanent_exclusion") != GSE251686_BAD_GSM:
        raise ValueError("GSE251686 supplementary outputs do not retain the permanent exclusion")
    if Path(manifest.get("generated_artifact_path_root", "")).resolve() != results_root:
        raise ValueError("GSE251686 supplementary manifest has an unexpected artifact root")
    expected_artifacts = {
        "supplementary_tables\\README.md",
        "supplementary_tables\\supplementary_table_s1_gse251686_exploratory_effects.csv",
        "supplementary_figures\\supplementary_figure_s1_gse251686_exploratory_sensitivity.pdf",
        "supplementary_figures\\supplementary_figure_s1_gse251686_exploratory_sensitivity.png",
    }
    observed_artifacts = set(manifest.get("generated_artifact_sha256", {}))
    if observed_artifacts != expected_artifacts:
        raise ValueError("GSE251686 supplementary manifest has an unexpected artifact set")
    for relative, expected_hash in manifest["generated_artifact_sha256"].items():
        path = results_root / relative
        if not path.is_file() or sha256(path) != expected_hash:
            raise ValueError(f"GSE251686 supplementary artifact hash mismatch: {relative}")
    return manifest


def require_discovery_sensitivity(sensitivity_dir: Path) -> dict[str, pd.DataFrame]:
    run_summary = read_csv(sensitivity_dir / "threshold_run_summary.csv")
    stability = read_csv(sensitivity_dir / "threshold_effect_stability_vs_30.csv")
    identity = read_csv(sensitivity_dir / "threshold_score_identity.csv")
    if set(run_summary["retained_cell_threshold"]) != {"20", "30", "50"}:
        raise ValueError("Discovery sensitivity package must include 20, 30, and 50 cell gates")
    if not run_summary["libraries_scored"].eq("24").all() or not run_summary["libraries_passing_threshold"].eq("24").all():
        raise ValueError("Discovery sensitivity package does not retain all 24 libraries")
    if not run_summary["minimum_observed_source_restricted_qc_passing_cells"].eq("471").all():
        raise ValueError("Discovery sensitivity minimum retained-cell count differs from the audited result")
    if not stability["effect_delta_current_minus_30"].eq("0").all() or not stability["effect_directions_agree"].eq("true").all():
        raise ValueError("Discovery threshold sensitivity is not identical to the 30-cell reference")
    required_identity_columns = {
        "matched_score_rows",
        "reference_only_score_rows",
        "current_only_score_rows",
        "maximum_absolute_module_score_delta_vs_30",
        "included_cell_counts_identical_vs_30",
        "score_identity_status",
    }
    missing_identity = required_identity_columns - set(identity.columns)
    if missing_identity:
        raise ValueError(f"Discovery score-identity audit is missing columns: {sorted(missing_identity)}")
    if (
        not identity["reference_only_score_rows"].eq("0").all()
        or not identity["current_only_score_rows"].eq("0").all()
        or not identity["maximum_absolute_module_score_delta_vs_30"].eq("0").all()
        or not identity["included_cell_counts_identical_vs_30"].eq("true").all()
        or not identity["score_identity_status"].eq("identical_score_keys_and_values").all()
    ):
        raise ValueError("Discovery sensitivity score identity does not match the 30-cell reference")
    return {"run_summary": run_summary, "stability": stability, "identity": identity}


def make_cohort_disposition_table(
    availability: pd.DataFrame,
    crosswalk: pd.DataFrame,
    gse_ledger: pd.DataFrame,
) -> list[dict[str, Any]]:
    crosswalk_by_cohort = crosswalk.set_index("cohort_id")
    rows: list[dict[str, Any]] = []
    ordered = [*DEFAULT_COHORT_ORDER, GSE251686_ID]
    for cohort_id in ordered:
        row = availability.loc[availability["cohort_id"].eq(cohort_id)].iloc[0]
        if cohort_id in crosswalk_by_cohort.index:
            identity = crosswalk_by_cohort.loc[cohort_id]
            keys_after_score = identity["unique_retained_score_sample_keys"]
            identity_check = f"{identity['crosswalk_status']}; {identity['matched_keys']}/{identity['unique_retained_score_sample_keys']} matched"
        else:
            keys_after_score = str(len(gse_ledger))
            identity_check = (
                "separate audit-gated score ledger; 5 selected presumed sample/library keys; "
                f"{GSE251686_BAD_GSM} permanently excluded"
            )
        rows.append(
            {
                "Analysis stream": "Default descriptive summary" if cohort_id != GSE251686_ID else "Separate exploratory sensitivity",
                "Cohort or parent project": row["dataset_or_parent_project"],
                "Compartment scope": row["compartment_scope"],
                "Observation key": row["sample_key_type"],
                "Recorded group structure": row["audited_group_structure"],
                "Keys after score filter": keys_after_score,
                "Identity or eligibility audit": identity_check,
                "Included in default 20-effect summary": row["effect_summary_included"],
                "Interpretation boundary": row["analysis_boundary"],
            }
        )
    return rows


def make_program_table(ledger: pd.DataFrame) -> list[dict[str, Any]]:
    required = {
        "module_id",
        "label",
        "source_class",
        "source_ids",
        "gene_count",
        "gene_symbols_sorted",
        "gene_list_sha256",
        "score_direction",
        "locked_at_utc",
    }
    missing = required - set(ledger.columns)
    if missing:
        raise ValueError(f"Program ledger is missing columns: {sorted(missing)}")
    if set(ledger["module_id"]) != set(MODULE_ORDER):
        raise ValueError("Program ledger does not contain the locked module set")
    ledger = ledger.set_index("module_id").loc[MODULE_ORDER].reset_index()
    return [
        {
            "Module ID": row["module_id"],
            "Module label": row["label"],
            "Source class": row["source_class"],
            "Source identifiers": row["source_ids"],
            "Gene count": row["gene_count"],
            "Locked gene symbols": row["gene_symbols_sorted"],
            "Gene-list SHA-256": row["gene_list_sha256"],
            "Score direction": row["score_direction"],
            "Locked at UTC": row["locked_at_utc"],
        }
        for _, row in ledger.iterrows()
    ]


def draw_box(ax: plt.Axes, x: float, y: float, width: float, height: float, text: str, color: str, edge: str = "#333333", fontsize: float = 8.0, linestyle: str = "-") -> None:
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.015",
        linewidth=1.15,
        edgecolor=edge,
        facecolor=color,
        linestyle=linestyle,
        transform=ax.transAxes,
        zorder=2,
    )
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height / 2, text, ha="center", va="center", fontsize=fontsize, transform=ax.transAxes, wrap=True, zorder=3)


def draw_arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str = "#4D4D4D", linestyle: str = "-") -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=1.15,
            color=color,
            linestyle=linestyle,
            transform=ax.transAxes,
            zorder=1,
        )
    )


def make_workflow_figure(output_dir: Path) -> tuple[Path, Path]:
    figure, ax = plt.subplots(figsize=(10.8, 6.0))
    ax.set_axis_off()
    draw_box(
        ax,
        0.035,
        0.72,
        0.19,
        0.16,
        "Public human IVDD\nGEO records and matrices",
        "#E6F0FA",
        fontsize=9,
    )
    draw_box(
        ax,
        0.285,
        0.72,
        0.20,
        0.16,
        "Provenance, matrix,\nidentity, QC, annotation,\nand mapping audits",
        "#F0F0F0",
        fontsize=8.5,
    )
    draw_arrow(ax, (0.225, 0.80), (0.285, 0.80))
    draw_box(
        ax,
        0.56,
        0.62,
        0.39,
        0.26,
        "Default descriptive summary\nGSE230809 parent project, GSE244889, GSE153066, GSE165722\n55 exact score-to-ledger matches; 20 effects\n(all confirmatory_eligible = false)",
        "#E7F4EA",
        fontsize=8.3,
    )
    draw_arrow(ax, (0.485, 0.80), (0.56, 0.76))
    draw_box(
        ax,
        0.56,
        0.24,
        0.39,
        0.23,
        "Cohort-specific differences, Welch/bootstrap intervals,\nleave-one-key-out stability, and descriptive sign alignment\nNo p-values, pooling, replication adjudication, or causal interpretation",
        "#FFF4D6",
        fontsize=8.1,
    )
    draw_arrow(ax, (0.755, 0.62), (0.755, 0.47))
    draw_box(
        ax,
        0.035,
        0.25,
        0.45,
        0.25,
        "Separate GSE251686 exploratory sensitivity\nmild n=2; severe n=3; five presumed sample/library keys\nGSM7986002 permanently excluded after stream-integrity failure\nNot a default-summary, sign-alignment, or main-figure input",
        "#FDE9E7",
        edge="#A33A2B",
        fontsize=8.2,
        linestyle="--",
    )
    draw_arrow(ax, (0.385, 0.72), (0.385, 0.50), linestyle="--")
    ax.text(0.035, 0.96, "Supplementary Figure S2. Cohort disposition and analysis boundary", fontsize=12, fontweight="bold", transform=ax.transAxes)
    ax.text(
        0.035,
        0.08,
        "Cells were nested observations throughout. GSE251686 was independently scored for transparent sensitivity reporting, but remained isolated because its usable records were incomplete and non-balanced.",
        fontsize=8.5,
        transform=ax.transAxes,
        wrap=True,
    )
    figure.tight_layout()
    pdf = output_dir / "supplementary_figure_s2_cohort_disposition_and_analysis_boundary.pdf"
    png = output_dir / "supplementary_figure_s2_cohort_disposition_and_analysis_boundary.png"
    figure.savefig(pdf, bbox_inches="tight")
    figure.savefig(png, dpi=600, bbox_inches="tight")
    plt.close(figure)
    return pdf, png


def make_graphical_abstract(output_dir: Path) -> tuple[Path, Path]:
    figure, ax = plt.subplots(figsize=(12.0, 5.2))
    ax.set_axis_off()
    ax.text(
        0.03,
        0.93,
        "Cohort-aware IVDD program audit: descriptive alignment and heterogeneity",
        fontsize=16,
        fontweight="bold",
        transform=ax.transAxes,
    )
    draw_box(
        ax,
        0.04,
        0.36,
        0.22,
        0.36,
        "Human public IVDD cohorts\n\nDonor or presumed\nsample/library key\nis the observation unit\n\nCells remain nested",
        "#E6F0FA",
        fontsize=10,
    )
    draw_arrow(ax, (0.26, 0.54), (0.34, 0.54))
    draw_box(
        ax,
        0.34,
        0.36,
        0.25,
        0.36,
        "Locked before external scoring\n\nECM/collagen\nInflammatory/NF-kB\nHypoxia/oxidative stress\nDisc matrix homeostasis",
        "#E7F4EA",
        fontsize=9.5,
    )
    draw_arrow(ax, (0.59, 0.54), (0.67, 0.54))
    draw_box(
        ax,
        0.67,
        0.36,
        0.29,
        0.36,
        "Four default NP contrasts\n\nHypoxia/oxidative-stress\npoint estimates positive in 4/4\n(all Welch intervals include zero)\n\nOther three modules: discordant directions",
        "#FFF4D6",
        fontsize=9.2,
    )
    draw_box(
        ax,
        0.68,
        0.11,
        0.28,
        0.18,
        "Cohort-specific descriptive evidence\nnot a universal program, mechanism, biomarker, or target",
        "#FDE9E7",
        edge="#A33A2B",
        fontsize=7.8,
    )
    draw_arrow(ax, (0.815, 0.36), (0.815, 0.29), color="#A33A2B")
    ax.text(
        0.04,
        0.09,
        "Separate sensitivity: GSE251686 mild n=2 versus severe n=3 after permanent exclusion of GSM7986002; not a default-summary input.",
        fontsize=7.6,
        transform=ax.transAxes,
        color="#6B2D25",
        wrap=True,
    )
    figure.tight_layout()
    pdf = output_dir / "graphical_abstract_cohort_aware_ivdd.pdf"
    png = output_dir / "graphical_abstract_cohort_aware_ivdd.png"
    figure.savefig(pdf, bbox_inches="tight")
    figure.savefig(png, dpi=600, bbox_inches="tight")
    plt.close(figure)
    return pdf, png


def write_support_readme(path: Path) -> None:
    path.write_text(
        "# Submission-Support Deliverables\n\n"
        "This directory contains publication-support tables and figures generated from audited, existing analysis outputs. "
        "The frozen default descriptive summary retains 20 effects across four score cohorts and remains the only authoritative "
        "main-analysis result. It excludes the separately audited GSE251686 package, which is represented only by "
        "Supplementary Table S1 and Supplementary Figure S1.\n\n"
        "Supplementary Table S2 records cohort disposition, observation keys, identity checks, and default-summary eligibility. "
        "Supplementary Table S3 records the locked program definitions and hashes. Supplementary Table S4 contains all "
        "default leave-one-key-out effects. Supplementary Tables S5a and S5b record the discovery retained-cell threshold "
        "sensitivity. Supplementary Table S6 indexes the reproducibility contract. Together, S1--S6 document the unchanged "
        "descriptive analysis, which has no p-values or hypothesis tests.\n\n"
        "Supplementary Tables S7a--S7d and Supplementary Figure S3 are a separate exploratory random-effects meta-analysis "
        "of the four default NP cohorts. Supplementary Tables S8a--S8d and Supplementary Figure S4 are a separate post hoc "
        "six-cohort expansion that adds GSE186542 and the GSE167931 FPKM representation. Supplementary Tables S9a--S9d and "
        "Supplementary Figure S5 are a source-family replacement sensitivity in which the native GSE245147 3-versus-3 clinical "
        "comparison replaces GSE167931. GSE167931 and GSE245147 are never pooled together. S7-S9 report HKSJ and four-module "
        "BH p-values only for transparent description; none is confirmatory, patient-level validation, replication, mechanism, "
        "biomarker, causal, or therapy evidence.\n\n"
        "Supplementary Figure S2 describes the default analysis streams and their separation. The graphical abstract is an "
        "evidence-bound visual summary, not a mechanistic model. The top-level submission-support manifest hashes the artifacts "
        "it lists and preserves the separate default and supplementary analysis contracts.\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--default-summary-dir", type=Path, required=True)
    parser.add_argument("--gse251686-package-dir", type=Path, required=True)
    parser.add_argument("--program-ledger", type=Path, required=True)
    parser.add_argument("--discovery-sensitivity-dir", type=Path, required=True)
    parser.add_argument("--results-root", type=Path, required=True)
    args = parser.parse_args()

    summary_dir = args.default_summary_dir.resolve()
    gse_dir = args.gse251686_package_dir.resolve()
    ledger_path = args.program_ledger.resolve()
    sensitivity_dir = args.discovery_sensitivity_dir.resolve()
    results_root = args.results_root.resolve()
    project_root = results_root.parent
    table_dir = results_root / "supplementary_tables"
    figure_dir = results_root / "supplementary_figures"
    graphical_dir = results_root / "graphical_abstract"
    main_table_dir = results_root / "tables"
    main_figure_dir = results_root / "figures"
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    graphical_dir.mkdir(parents=True, exist_ok=True)

    main_artifact_paths = [
        main_table_dir / "table_1_current_cohort_roles.csv",
        main_table_dir / "table_2_np_module_effects.csv",
        main_figure_dir / "figure_1_np_cohort_module_effects.pdf",
        main_figure_dir / "figure_1_np_cohort_module_effects.png",
        main_figure_dir / "figure_2_np_direction_alignment.pdf",
        main_figure_dir / "figure_2_np_direction_alignment.png",
    ]
    missing_main_artifacts = [path for path in main_artifact_paths if not path.is_file()]
    if missing_main_artifacts:
        raise FileNotFoundError(
            "Missing main deliverables; run make_current_summary_deliverables.py first: "
            + ", ".join(str(path) for path in missing_main_artifacts)
        )

    default = require_default_summary(summary_dir)
    gse = require_gse251686_package(gse_dir)
    gse_supplement = require_gse251686_supplement(table_dir, results_root)
    sensitivity = require_discovery_sensitivity(sensitivity_dir)
    meta_packages = {
        key: require_meta_analysis_package(
            project_root,
            project_root / "data" / "derived" / spec["directory"],
            spec,
        )
        for key, spec in META_ANALYSIS_SPECS.items()
    }
    program_ledger = read_csv(ledger_path)

    cohort_table = table_dir / "supplementary_table_s2_cohort_disposition_and_identity.csv"
    program_table = table_dir / "supplementary_table_s3_locked_program_modules.csv"
    loko_table = table_dir / "supplementary_table_s4_default_leave_one_key_out.csv"
    sensitivity_summary_table = table_dir / "supplementary_table_s5a_discovery_threshold_summary.csv"
    sensitivity_stability_table = table_dir / "supplementary_table_s5b_discovery_threshold_effect_stability.csv"
    reproducibility_table = table_dir / "supplementary_table_s6_reproducibility_contract.csv"
    readme = table_dir / "README_submission_support.md"
    workflow_pdf, workflow_png = make_workflow_figure(figure_dir)
    abstract_pdf, abstract_png = make_graphical_abstract(graphical_dir)

    cohort_rows = make_cohort_disposition_table(default["availability"], default["crosswalk"], gse["ledger"])
    write_csv(cohort_table, cohort_rows, list(cohort_rows[0]))
    program_rows = make_program_table(program_ledger)
    write_csv(program_table, program_rows, list(program_rows[0]))
    default_loko = read_csv(summary_dir / "donor_module_leave_one_out.csv")
    write_csv(loko_table, default_loko.to_dict("records"), list(default_loko.columns))
    write_csv(
        sensitivity_summary_table,
        sensitivity["run_summary"].to_dict("records"),
        list(sensitivity["run_summary"].columns),
    )
    write_csv(
        sensitivity_stability_table,
        sensitivity["stability"].to_dict("records"),
        list(sensitivity["stability"].columns),
    )

    source_paths = [
        ("default_summary_manifest", summary_dir / "run_manifest.json", "Default 20-effect contract and input/output hashes."),
        ("default_summary_artifacts", summary_dir / "run_artifacts.csv", "Default current-run artifact ledger."),
        ("default_score_identity_crosswalk", summary_dir / "score_ledger_identity_crosswalk.csv", "55/55 exact score-to-ledger matching audit."),
        ("locked_program_ledger", ledger_path, "Gene lists, source identifiers, locking time, and SHA-256 hashes."),
        ("discovery_threshold_input_hashes", sensitivity_dir / "input_artifact_hashes.csv", "Threshold sensitivity input hashes."),
        ("gse251686_exploratory_manifest", gse_dir / "run_manifest.json", "Separate exploratory package; excluded from default summary."),
        ("gse251686_supplement_manifest", table_dir / "run_manifest.json", "Separate GSE251686 sensitivity table/figure contract."),
        ("s7_exploratory_meta_analysis_manifest", meta_packages["s7_exploratory_four_cohort"]["manifest_path"], META_ANALYSIS_SPECS["s7_exploratory_four_cohort"]["purpose"]),
        ("s8_post_hoc_external_expansion_manifest", meta_packages["s8_post_hoc_external_expansion"]["manifest_path"], META_ANALYSIS_SPECS["s8_post_hoc_external_expansion"]["purpose"]),
        ("s9_source_family_replacement_manifest", meta_packages["s9_source_family_replacement"]["manifest_path"], META_ANALYSIS_SPECS["s9_source_family_replacement"]["purpose"]),
        ("main_deliverables_script", Path(__file__).with_name("make_current_summary_deliverables.py"), "Generator for the two main tables and four main figure files."),
        ("submission_support_script", Path(__file__).resolve(), "Generator for the submission-support tables, workflow figure, and graphical abstract."),
        ("gse251686_supplement_script", Path(__file__).with_name("make_gse251686_supplement.py"), "Generator for the isolated GSE251686 supplementary sensitivity table and figure."),
        ("s7_meta_analysis_script", Path(__file__).with_name("run_np_exploratory_meta_analysis.R"), "Generator for the separate S7 exploratory random-effects synthesis."),
        ("s8_s9_meta_analysis_script", Path(__file__).with_name("run_np_post_hoc_external_expansion_meta_analysis.R"), "Generator for the separate S8 expansion and S9 replacement sensitivity."),
        ("python_environment_lock", results_root.parent / "tools" / "python" / "requirements-lock.txt", "Project-local Python package lock."),
    ]
    reproducibility_rows = [
        {
            "Artifact role": role,
            "Project-relative path": str(path.relative_to(results_root.parent)),
            "SHA-256": sha256(path),
            "Purpose": purpose,
        }
        for role, path, purpose in source_paths
        if path.is_file()
    ]
    write_csv(reproducibility_table, reproducibility_rows, list(reproducibility_rows[0]))
    write_support_readme(readme)

    artifact_paths = [
        *main_artifact_paths,
        *[results_root / relative for relative in gse_supplement["generated_artifact_sha256"]],
        cohort_table,
        program_table,
        loko_table,
        sensitivity_summary_table,
        sensitivity_stability_table,
        reproducibility_table,
        readme,
        workflow_pdf,
        workflow_png,
        abstract_pdf,
        abstract_png,
    ]
    meta_artifact_paths = {
        key: package["artifact_paths"]
        for key, package in meta_packages.items()
    }
    for paths in meta_artifact_paths.values():
        artifact_paths.extend(path for path in paths if path.is_relative_to(results_root))
    artifact_paths = list(dict.fromkeys(artifact_paths))
    manifest_path = results_root / "submission_support_manifest.json"
    manifest = {
        "schema_version": 2,
        "default_summary_effect_count": 20,
        "default_summary_confirmatory_eligible": False,
        "default_score_ledger_exact_matches": "55/55",
        "gse251686_default_summary_inclusion": False,
        "gse251686_permanent_exclusion": GSE251686_BAD_GSM,
        "analysis_contract": {
            "default_summary": "Frozen 20-effect descriptive analysis; no p-values, formal meta-analysis, or replication adjudication.",
            "supplementary_s7_s9": "Separate non-confirmatory SMDH random-effects syntheses; HKSJ and four-module BH p-values are transparent descriptive outputs only.",
        },
        "default_no_hypothesis_tests_or_p_values": True,
        "supplementary_meta_analysis_p_values_transparency_only": True,
        "generated_artifact_sha256": {
            str(path.relative_to(results_root)): sha256(path)
            for path in artifact_paths
        },
        "meta_analysis_artifact_sha256": {
            key: {
                str(path.relative_to(project_root)): sha256(path)
                for path in paths
            }
            for key, paths in meta_artifact_paths.items()
        },
        "source_manifest_sha256": {
            "default_summary": sha256(summary_dir / "run_manifest.json"),
            "gse251686_exploratory": sha256(gse_dir / "run_manifest.json"),
            **{
                key: sha256(package["manifest_path"])
                for key, package in meta_packages.items()
            },
        },
        "self_referential_hash_exclusions": [str(manifest_path.relative_to(results_root))],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
