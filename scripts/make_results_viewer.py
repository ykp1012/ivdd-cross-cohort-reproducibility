"""Build a local HTML viewer for the audited IVDD figures and tables.

The viewer is intentionally a derived presentation artifact. It reads the
already-generated CSV/PNG/PDF files and does not recalculate any result.
PNG files are embedded so the page remains viewable when opened directly.
"""

from __future__ import annotations

import argparse
import base64
import csv
import html
from pathlib import Path


FIGURES = [
    ("Figure 1 - NP cohort module effects", "results/figures/figure_1_np_cohort_module_effects.png", "results/figures/figure_1_np_cohort_module_effects.pdf"),
    ("Figure 2 - NP direction alignment", "results/figures/figure_2_np_direction_alignment.png", "results/figures/figure_2_np_direction_alignment.pdf"),
    ("Supplementary Figure S1 - GSE251686 sensitivity", "results/supplementary_figures/supplementary_figure_s1_gse251686_exploratory_sensitivity.png", "results/supplementary_figures/supplementary_figure_s1_gse251686_exploratory_sensitivity.pdf"),
    ("Supplementary Figure S2 - cohort disposition", "results/supplementary_figures/supplementary_figure_s2_cohort_disposition_and_analysis_boundary.png", "results/supplementary_figures/supplementary_figure_s2_cohort_disposition_and_analysis_boundary.pdf"),
    ("Supplementary Figure S3 - four-cohort exploratory NP meta-analysis", "results/supplementary_figures/supplementary_figure_s3_np_exploratory_random_effects_meta_analysis.png", "results/supplementary_figures/supplementary_figure_s3_np_exploratory_random_effects_meta_analysis.pdf"),
    ("Supplementary Figure S4 - post-hoc six-cohort exploratory NP meta-analysis", "results/supplementary_figures/supplementary_figure_s4_np_post_hoc_external_expansion_meta_analysis.png", "results/supplementary_figures/supplementary_figure_s4_np_post_hoc_external_expansion_meta_analysis.pdf"),
    ("Supplementary Figure S5 - NP source-family replacement sensitivity", "results/supplementary_figures/supplementary_figure_s5_np_source_family_replacement_sensitivity.png", "results/supplementary_figures/supplementary_figure_s5_np_source_family_replacement_sensitivity.pdf"),
]


def esc(value: object) -> str:
    text = "" if value is None else str(value)
    return html.escape(text).replace("\n", "<br>")


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def relative_href(output: Path, target: Path) -> str:
    return target.resolve().relative_to(output.resolve().parent).as_posix()


def render_table(output: Path, title: str, source: Path, columns: list[str] | None = None) -> str:
    fieldnames, rows = read_csv(source)
    selected = [name for name in (columns or fieldnames) if name in fieldnames]
    header = "".join(f"<th scope='col'>{esc(name)}</th>" for name in selected)
    body = []
    for row in rows:
        cells = "".join(f"<td>{esc(row.get(name, ''))}</td>" for name in selected)
        body.append(f"<tr>{cells}</tr>")
    download = relative_href(output, source)
    return (
        f"<section class='table-section' id='{html.escape(title.lower().replace(' ', '-'))}'>"
        f"<div class='section-heading'><h3>{esc(title)}</h3>"
        f"<a class='file-link' href='{html.escape(download)}' download>Download CSV</a></div>"
        f"<div class='table-wrap'><table><thead><tr>{header}</tr></thead>"
        f"<tbody>{''.join(body)}</tbody></table></div></section>"
    )


def embedded_png(path: Path) -> str:
    payload = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{payload}"


def render_figure(output: Path, title: str, png: Path, pdf: Path) -> str:
    png_src = embedded_png(png)
    pdf_href = relative_href(output, pdf)
    png_href = relative_href(output, png)
    return (
        "<figure class='figure-block'>"
        f"<img src='{png_src}' alt='{esc(title)}' loading='lazy'>"
        f"<figcaption><strong>{esc(title)}</strong>"
        f"<span class='figure-links'><a href='{html.escape(png_href)}' download>PNG</a>"
        f"<a href='{html.escape(pdf_href)}' download>PDF</a></span></figcaption>"
        "</figure>"
    )


def build_page(output: Path) -> str:
    root = Path(__file__).resolve().parents[1]
    result_root = root / "results"
    main_table_1 = result_root / "tables/table_1_current_cohort_roles.csv"
    main_table_2 = result_root / "tables/table_2_np_module_effects.csv"
    s1 = result_root / "supplementary_tables/supplementary_table_s1_gse251686_exploratory_effects.csv"
    s2 = result_root / "supplementary_tables/supplementary_table_s2_cohort_disposition_and_identity.csv"
    s3 = result_root / "supplementary_tables/supplementary_table_s3_locked_program_modules.csv"
    s4 = result_root / "supplementary_tables/supplementary_table_s4_default_leave_one_key_out.csv"
    s5a = result_root / "supplementary_tables/supplementary_table_s5a_discovery_threshold_summary.csv"
    s7b = result_root / "supplementary_tables/supplementary_table_s7b_np_meta_analysis_primary_results.csv"
    s8b = result_root / "supplementary_tables/supplementary_table_s8b_np_post_hoc_external_expansion_primary_results.csv"
    s9b = result_root / "supplementary_tables/supplementary_table_s9b_np_source_family_replacement_sensitivity_primary_results.csv"

    table_1_columns = [
        "Cohort and role", "Recorded group structure", "Observation key",
        "Keys after score filter", "Exact ledger matches", "Identity check",
        "Interpretation boundary", "Confirmatory status",
    ]
    table_2_columns = [
        "Cohort and role", "Pre-specified module", "Lower vs higher group n",
        "Higher-minus-lower score difference", "Welch 95% interval",
        "Bootstrap 95% interval", "Observed direction",
        "LODO direction retention (fraction)",
    ]
    s1_columns = [
        "Dataset and analysis role", "Pre-specified module",
        "Recorded lower vs higher group n", "Higher-minus-lower score difference",
        "Welch 95% interval", "LOKO direction retention (fraction)",
        "Permanent exclusion", "Interpretation boundary",
    ]
    s2_columns = [
        "Analysis stream", "Cohort or parent project", "Compartment scope",
        "Observation key", "Recorded group structure", "Keys after score filter",
        "Included in default 20-effect summary", "Interpretation boundary",
    ]
    s3_columns = [
        "Module ID", "Module label", "Source class", "Source identifiers",
        "Gene count", "Gene-list SHA-256", "Score direction", "Locked at UTC",
    ]
    s4_columns = [
        "cohort_id", "compartment", "module_id", "contrast_label",
        "excluded_donor_or_library_key", "excluded_arm",
        "remaining_comparison_n", "remaining_target_n",
        "leave_one_out_effect_target_minus_comparison",
        "direction_agrees_with_full_effect", "leave_one_out_status",
    ]
    s5a_columns = [
        "retained_cell_threshold", "libraries_scored", "libraries_passing_threshold",
        "minimum_observed_source_restricted_qc_passing_cells", "threshold_run_status",
    ]
    meta_columns = [
        "module", "effect_measure", "analysis_scope", "k",
        "pooled_standardized_mean_difference", "ci_lower", "ci_upper",
        "I_squared_percent", "prediction_interval_lower", "prediction_interval_upper",
    ]

    table_html = "".join(
        [
            "<div class='table-group-heading'><h3>Default descriptive analysis (unchanged)</h3>"
            "<p>These are the only authoritative main-analysis results. The 20-effect descriptive summary and its four default score cohorts are unchanged.</p></div>",
            render_table(output, "Table 1 - cohort roles", main_table_1, table_1_columns),
            render_table(output, "Table 2 - NP module effects", main_table_2, table_2_columns),
            render_table(output, "Supplementary Table S1 - GSE251686 sensitivity", s1, s1_columns),
            render_table(output, "Supplementary Table S2 - cohort disposition", s2, s2_columns),
            render_table(output, "Supplementary Table S3 - locked program modules", s3, s3_columns),
            render_table(output, "Supplementary Table S4 - default leave-one-key-out", s4, s4_columns),
            render_table(output, "Supplementary Table S5a - threshold summary", s5a, s5a_columns),
            "<div class='table-group-heading exploratory'><h3>Separate exploratory meta-analysis packages</h3>"
            "<p>These supplementary results do not alter, replace, or validate the unchanged 20-effect descriptive main analysis.</p></div>",
            render_table(output, "Supplementary Table S7b - exploratory four-cohort NP meta-analysis", s7b, meta_columns),
            render_table(output, "Supplementary Table S8b - post-hoc six-cohort NP expansion", s8b, meta_columns),
            render_table(output, "Supplementary Table S9b - NP source-family replacement sensitivity", s9b, meta_columns),
        ]
    )

    figure_html = "".join(
        render_figure(output, title, root / png, root / pdf)
        for title, png, pdf in FIGURES
    )
    support_links = []
    for label, rel in [
        ("Supplementary Table S5b - discovery threshold effect stability", "results/supplementary_tables/supplementary_table_s5b_discovery_threshold_effect_stability.csv"),
        ("Supplementary Table S6 - reproducibility contract", "results/supplementary_tables/supplementary_table_s6_reproducibility_contract.csv"),
        ("Supplementary Table S7a - four-cohort study effects", "results/supplementary_tables/supplementary_table_s7a_np_meta_analysis_study_effects.csv"),
        ("Supplementary Table S7c - four-cohort model sensitivity", "results/supplementary_tables/supplementary_table_s7c_np_meta_analysis_model_sensitivity.csv"),
        ("Supplementary Table S7d - four-cohort leave-one-cohort-out", "results/supplementary_tables/supplementary_table_s7d_np_meta_analysis_leave_one_cohort_out.csv"),
        ("Supplementary Table S8a - post-hoc six-cohort study effects", "results/supplementary_tables/supplementary_table_s8a_np_post_hoc_external_expansion_study_effects.csv"),
        ("Supplementary Table S8c - post-hoc six-cohort model sensitivity", "results/supplementary_tables/supplementary_table_s8c_np_post_hoc_external_expansion_model_sensitivity.csv"),
        ("Supplementary Table S8d - post-hoc six-cohort leave-one-cohort-out", "results/supplementary_tables/supplementary_table_s8d_np_post_hoc_external_expansion_leave_one_cohort_out.csv"),
        ("Supplementary Table S9a - source-family replacement study effects", "results/supplementary_tables/supplementary_table_s9a_np_source_family_replacement_sensitivity_study_effects.csv"),
        ("Supplementary Table S9c - source-family replacement model sensitivity", "results/supplementary_tables/supplementary_table_s9c_np_source_family_replacement_sensitivity_model_sensitivity.csv"),
        ("Supplementary Table S9d - source-family replacement leave-one-cohort-out", "results/supplementary_tables/supplementary_table_s9d_np_source_family_replacement_sensitivity_leave_one_cohort_out.csv"),
        ("S7 meta-analysis manifest", "data/derived/np_exploratory_meta_analysis/meta_analysis_manifest.json"),
        ("S8 meta-analysis manifest", "data/derived/np_post_hoc_external_expansion_meta_analysis/meta_analysis_manifest.json"),
        ("S9 source-family replacement manifest", "data/derived/np_source_family_replacement_meta_analysis/meta_analysis_manifest.json"),
        ("Supplementary-table index and analysis boundaries", "results/supplementary_tables/README.md"),
        ("Submission manifest", "results/submission_support_manifest.json"),
        ("Full manuscript draft", "manuscript/04_manuscript_draft.md"),
        ("Submission checklist", "manuscript/05_submission_readiness_checklist.md"),
    ]:
        support_links.append(
            f"<li><a href='{html.escape(relative_href(output, root / rel))}'>{esc(label)}</a></li>"
        )

    css = """
    :root { --ink:#172033; --muted:#5c677d; --line:#dfe5ef; --soft:#f5f7fb; --accent:#1f5fae; --accent2:#0f766e; }
    * { box-sizing:border-box; }
    body { margin:0; color:var(--ink); background:#fff; font:15px/1.6 Arial, sans-serif; }
    header { background:#17324d; color:#fff; padding:32px clamp(18px,4vw,56px); }
    header h1 { max-width:1180px; margin:0 auto 8px; font-size:clamp(24px,3vw,38px); }
    header p { max-width:1180px; margin:0 auto; color:#d7e3ef; }
    nav { position:sticky; top:0; z-index:5; background:#fff; border-bottom:1px solid var(--line); padding:10px clamp(18px,4vw,56px); }
    nav div { max-width:1180px; margin:auto; display:flex; gap:18px; flex-wrap:wrap; }
    nav a, .file-link, .figure-links a { color:var(--accent); text-decoration:none; font-weight:600; }
    nav a:hover, .file-link:hover, .figure-links a:hover { text-decoration:underline; }
    main { max-width:1180px; margin:auto; padding:28px clamp(18px,4vw,56px) 64px; }
    section { margin:0 0 34px; }
    h2 { margin:0 0 12px; font-size:24px; border-bottom:2px solid var(--line); padding-bottom:7px; }
    h3 { margin:0; font-size:18px; }
    .summary { display:grid; grid-template-columns:repeat(auto-fit,minmax(180px,1fr)); gap:12px; margin:18px 0 24px; }
    .metric { border:1px solid var(--line); background:var(--soft); padding:14px; }
    .metric strong { display:block; font-size:22px; color:var(--accent2); }
    .figure-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(320px,1fr)); gap:20px; }
    .figure-block { margin:0; border:1px solid var(--line); padding:10px; background:#fff; }
    .figure-block img { display:block; width:100%; height:auto; max-height:640px; object-fit:contain; background:#fff; }
    figcaption { padding:10px 4px 2px; display:flex; justify-content:space-between; gap:12px; flex-wrap:wrap; }
    .figure-links { display:flex; gap:12px; font-size:13px; }
    .table-group-heading { margin:28px 0 12px; padding:12px 14px; border-left:4px solid var(--accent); background:var(--soft); }
    .table-group-heading h3 { margin:0 0 4px; }
    .table-group-heading p { margin:0; color:var(--muted); }
    .table-group-heading.exploratory { border-left-color:#b45309; background:#fff8eb; }
    .table-section { border-top:1px solid var(--line); padding-top:16px; }
    .section-heading { display:flex; justify-content:space-between; gap:14px; align-items:baseline; flex-wrap:wrap; margin-bottom:10px; }
    .table-wrap { overflow-x:auto; border:1px solid var(--line); }
    table { border-collapse:collapse; width:100%; min-width:760px; font-size:13px; }
    th, td { border-bottom:1px solid var(--line); padding:8px 10px; text-align:left; vertical-align:top; }
    th { background:#edf3f9; font-weight:700; white-space:nowrap; }
    tr:nth-child(even) td { background:#fbfcfe; }
    .note { border-left:4px solid var(--accent2); background:#eef8f5; padding:12px 14px; }
    footer { border-top:1px solid var(--line); color:var(--muted); max-width:1180px; margin:auto; padding:18px clamp(18px,4vw,56px) 40px; }
    @media (max-width:700px) { nav div { gap:10px; } .figure-grid { grid-template-columns:1fr; } }
    """
    return f"""<!doctype html>
<html lang='en'>
<head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'>
<title>IVDD Results Viewer</title><style>{css}</style></head>
<body>
<header><h1>IVDD Cross-Cohort Results Viewer</h1><p>Audited default descriptive results and separately labeled exploratory supplementary analyses.</p></header>
<nav><div><a href='#overview'>Overview</a><a href='#analysis-boundary'>Analysis boundary</a><a href='#figures'>Figures</a><a href='#tables'>Tables</a><a href='#files'>Files</a></div></nav>
<main>
<section id='overview'><h2>Overview</h2>
<p class='note'>The Markdown preview may show only text. This page embeds the PNG figures and renders the main CSV tables as HTML so they remain visible when the file is opened directly.</p>
<div class='summary'><div class='metric'><strong>20</strong>default descriptive effects</div><div class='metric'><strong>16</strong>default NP effects</div><div class='metric'><strong>55/55</strong>default identity matches</div><div class='metric'><strong>0</strong>default p-values or hypothesis tests</div></div>
<p>The four default NP hypoxia/oxidative-stress point estimates are positive, but every descriptive Welch interval crosses zero. GSE251686 is isolated; GSM7986002 remains permanently excluded.</p>
</section>
<section id='analysis-boundary'><h2>Analysis Boundary</h2>
<p class='note'><strong>Default 20-effect descriptive analysis unchanged.</strong> The four default score cohorts remain the only authoritative main-analysis evidence on this page; Supplementary Tables S3 and S4 document its locked modules and leave-one-key-out stability.</p>
<p><strong>Supplementary Table S7:</strong> a separate exploratory random-effects meta-analysis of the four default NP cohorts. It is a quantitative synthesis only and does not replace the descriptive main analysis or establish confirmation, replication, mechanism, biomarker, or therapy evidence.</p>
<p><strong>Supplementary Table S8:</strong> a separate post-hoc exploratory six-cohort expansion that adds GSE186542 and the GSE167931 FPKM representation. It is not an independent patient-level validation and does not change the default analysis.</p>
<p><strong>Supplementary Table S9:</strong> a source-family replacement sensitivity that substitutes the native GSE245147 3-versus-3 comparison for GSE167931. The related source-family datasets are never pooled together, and S9 is not an additional independent validation.</p>
</section>
<section id='figures'><h2>Figures</h2><div class='figure-grid'>{figure_html}</div></section>
<section id='tables'><h2>Tables</h2>{table_html}</section>
<section id='files'><h2>More Files</h2><ul>{''.join(support_links)}</ul></section>
</main><footer>Generated from existing IVDD result artifacts. No analysis was recalculated by this viewer.</footer>
</body></html>"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = (args.output or root / "ivdd_results_viewer.html").resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(build_page(output), encoding="utf-8")
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
