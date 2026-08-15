#!/usr/bin/env Rscript

# Run a deliberately bounded, cohort-level exploratory meta-analysis for the
# four default NP cohorts. The existing donor/library-level Welch mean
# differences remain the primary results. This script standardizes within each
# cohort only because the supplied score scales differ between cohorts.

required_packages <- c("metafor", "digest", "jsonlite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(metafor)
  library(digest)
  library(jsonlite)
})

MODULE_ORDER <- c(
  "ecm_collagen_remodeling",
  "inflammatory_nfkb",
  "hypoxia_oxidative_stress",
  "disc_matrix_homeostasis"
)

MODULE_LABELS <- c(
  ecm_collagen_remodeling = "ECM / collagen remodeling",
  inflammatory_nfkb = "Inflammatory / NF-kB",
  hypoxia_oxidative_stress = "Hypoxia / oxidative stress",
  disc_matrix_homeostasis = "Disc matrix homeostasis"
)

COHORT_ORDER <- c(
  "GSE230809_discovery",
  "GSE244889_directional",
  "GSE153066_support",
  "GSE165722_score_level"
)

COHORT_LABELS <- c(
  GSE230809_discovery = "GSE230809 parent project (NP 3 vs 8)",
  GSE244889_directional = "GSE244889 (NP 4 vs 3)",
  GSE153066_support = "GSE153066 (NP 8 vs 8)",
  GSE165722_score_level = "GSE165722 (NP 4 vs 4)"
)

usage <- function() {
  paste(
    "Usage:",
    "Rscript scripts/run_np_exploratory_meta_analysis.R",
    "--summary-dir data/derived/donor_module_effect_summary",
    "--table-dir results/supplementary_tables",
    "--figure-dir results/supplementary_figures",
    "--output-dir data/derived/np_exploratory_meta_analysis",
    sep = "\n"
  )
}

parse_args <- function(args) {
  expected <- c("summary-dir", "table-dir", "figure-dir", "output-dir")
  values <- setNames(rep(NA_character_, length(expected)), expected)
  index <- 1L
  while (index <= length(args)) {
    token <- args[[index]]
    if (identical(token, "--help") || identical(token, "-h")) {
      cat(usage(), "\n")
      quit(status = 0L)
    }
    if (!startsWith(token, "--")) {
      stop("Unexpected argument: ", token, "\n", usage(), call. = FALSE)
    }
    key <- substring(token, 3L)
    if (!key %in% expected) {
      stop("Unknown option: ", token, "\n", usage(), call. = FALSE)
    }
    if (index == length(args)) {
      stop("Missing value for option: ", token, call. = FALSE)
    }
    values[[key]] <- args[[index + 1L]]
    index <- index + 2L
  }
  if (anyNA(values) || any(!nzchar(values))) {
    stop("All four output/input options are required.\n", usage(), call. = FALSE)
  }
  values
}

assert_columns <- function(data, expected, source_name) {
  missing <- setdiff(expected, names(data))
  if (length(missing) > 0L) {
    stop(source_name, " is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

as_number <- function(data, columns, source_name) {
  for (column in columns) {
    data[[column]] <- suppressWarnings(as.numeric(data[[column]]))
    if (any(!is.finite(data[[column]]))) {
      stop(source_name, " has a non-finite numeric value in ", column, call. = FALSE)
    }
  }
  data
}

sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

manifest_path <- function(path) {
  normalized <- gsub("\\\\", "/", path)
  sub("^\\./", "", normalized)
}

script_path <- function() {
  command_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(command_file) == 1L) {
    return(sub("^--file=", "", command_file))
  }
  NA_character_
}

make_wide_groups <- function(group_file, effect_file) {
  groups <- read.csv(group_file, stringsAsFactors = FALSE, check.names = FALSE)
  effects <- read.csv(effect_file, stringsAsFactors = FALSE, check.names = FALSE)

  assert_columns(
    groups,
    c(
      "cohort_id", "compartment", "module_id", "contrast_label", "contrast_arm",
      "severity_group", "n_donor_or_library_keys", "mean_module_score_log1p_cpm",
      "sd_module_score_log1p_cpm"
    ),
    basename(group_file)
  )
  assert_columns(
    effects,
    c("cohort_id", "compartment", "module_id", "confirmatory_eligible", "cohort_limitations"),
    basename(effect_file)
  )

  groups <- groups[
    groups$cohort_id %in% COHORT_ORDER &
      groups$compartment == "NP" &
      groups$module_id %in% MODULE_ORDER,
    ,
    drop = FALSE
  ]
  effects <- effects[
    effects$cohort_id %in% COHORT_ORDER &
      effects$compartment == "NP" &
      effects$module_id %in% MODULE_ORDER,
    ,
    drop = FALSE
  ]

  expected_groups <- length(COHORT_ORDER) * length(MODULE_ORDER) * 2L
  expected_effects <- length(COHORT_ORDER) * length(MODULE_ORDER)
  if (nrow(groups) != expected_groups) {
    stop("Expected ", expected_groups, " default NP group rows, found ", nrow(groups), call. = FALSE)
  }
  if (nrow(effects) != expected_effects) {
    stop("Expected ", expected_effects, " default NP effect rows, found ", nrow(effects), call. = FALSE)
  }
  if (!all(effects$confirmatory_eligible == "false")) {
    stop("The current project boundary requires every pooled cohort to remain non-confirmatory.", call. = FALSE)
  }

  groups <- as_number(
    groups,
    c("n_donor_or_library_keys", "mean_module_score_log1p_cpm", "sd_module_score_log1p_cpm"),
    basename(group_file)
  )
  if (any(groups$n_donor_or_library_keys < 2L) || any(groups$sd_module_score_log1p_cpm <= 0)) {
    stop("Every meta-analysis arm must have n >= 2 and a positive sample SD.", call. = FALSE)
  }

  key_columns <- c("cohort_id", "module_id", "contrast_label")
  target <- groups[groups$contrast_arm == "target", c(key_columns, "severity_group", "n_donor_or_library_keys", "mean_module_score_log1p_cpm", "sd_module_score_log1p_cpm")]
  comparison <- groups[groups$contrast_arm == "comparison", c(key_columns, "severity_group", "n_donor_or_library_keys", "mean_module_score_log1p_cpm", "sd_module_score_log1p_cpm")]
  if (nrow(target) != expected_effects || nrow(comparison) != expected_effects) {
    stop("Every default NP cohort/module contrast must have exactly one target and one comparison arm.", call. = FALSE)
  }

  names(target)[(length(key_columns) + 1L):ncol(target)] <- c("target_group", "target_n", "target_mean", "target_sd")
  names(comparison)[(length(key_columns) + 1L):ncol(comparison)] <- c("comparison_group", "comparison_n", "comparison_mean", "comparison_sd")
  wide <- merge(target, comparison, by = key_columns, all = TRUE, sort = FALSE)
  if (nrow(wide) != expected_effects || any(!complete.cases(wide))) {
    stop("Failed to form complete target/comparison arms for every default NP contrast.", call. = FALSE)
  }

  effect_meta <- effects[, c("cohort_id", "module_id", "cohort_limitations")]
  wide <- merge(wide, effect_meta, by = c("cohort_id", "module_id"), all.x = TRUE, sort = FALSE)
  if (any(is.na(wide$cohort_limitations))) {
    stop("A default NP contrast lacks its cohort limitation text.", call. = FALSE)
  }

  wide$cohort_id <- factor(wide$cohort_id, levels = COHORT_ORDER)
  wide$module_id <- factor(wide$module_id, levels = MODULE_ORDER)
  wide <- wide[order(wide$module_id, wide$cohort_id), , drop = FALSE]
  wide$cohort_id <- as.character(wide$cohort_id)
  wide$module_id <- as.character(wide$module_id)
  wide$raw_mean_difference_target_minus_comparison <- wide$target_mean - wide$comparison_mean
  wide
}

calculate_effects <- function(wide, measure) {
  calculated <- metafor::escalc(
    measure = measure,
    m1i = wide$target_mean,
    sd1i = wide$target_sd,
    n1i = wide$target_n,
    m2i = wide$comparison_mean,
    sd2i = wide$comparison_sd,
    n2i = wide$comparison_n
  )
  if (any(!is.finite(calculated$yi)) || any(!is.finite(calculated$vi)) || any(calculated$vi <= 0)) {
    stop("metafor::escalc returned an invalid standardized effect or variance for ", measure, call. = FALSE)
  }
  data.frame(effect = calculated$yi, variance = calculated$vi, standard_error = sqrt(calculated$vi))
}

fit_random_effects <- function(effect_data, method) {
  arguments <- list(
    yi = effect_data$effect,
    vi = effect_data$variance,
    data = effect_data,
    method = method,
    test = "knha",
    level = 95
  )
  if (identical(method, "REML")) {
    arguments$control <- list(maxiter = 10000)
  }
  fit <- do.call(metafor::rma.uni, arguments)
  attr(fit, "fit_control_note") <- if (identical(method, "REML")) {
    "metafor REML maxiter=10000; all other control values default"
  } else {
    "default_metafor_control"
  }
  fit
}

model_row <- function(fit, module_id, measure_label, tau_method, scope_label, omitted_cohort = NA_character_) {
  prediction <- predict(fit)
  fit_control_note <- attr(fit, "fit_control_note", exact = TRUE)
  if (is.null(fit_control_note)) {
    fit_control_note <- "default_metafor_control"
  }
  data.frame(
    module_id = module_id,
    module = unname(MODULE_LABELS[[module_id]]),
    effect_measure = measure_label,
    tau_squared_method = tau_method,
    confidence_interval_method = "Knapp-Hartung",
    analysis_scope = scope_label,
    omitted_cohort = omitted_cohort,
    k = fit$k,
    pooled_standardized_mean_difference = as.numeric(fit$b[1L]),
    ci_lower = fit$ci.lb,
    ci_upper = fit$ci.ub,
    hksj_test_statistic = fit$zval,
    hksj_degrees_freedom = fit$dfs,
    hksj_p_value = fit$pval,
    tau_squared = fit$tau2,
    I_squared_percent = fit$I2,
    Cochran_Q = fit$QE,
    Cochran_Q_degrees_freedom = fit$k - 1L,
    Cochran_Q_p_value = fit$QEp,
    prediction_interval_lower = prediction$pi.lb,
    prediction_interval_upper = prediction$pi.ub,
    fit_control_note = fit_control_note,
    stringsAsFactors = FALSE
  )
}

run_models <- function(study_effects, effect_column, variance_column, measure_label, tau_method, scope_label) {
  rows <- list()
  for (module_id in MODULE_ORDER) {
    data <- study_effects[study_effects$module_id == module_id, , drop = FALSE]
    data$effect <- data[[effect_column]]
    data$variance <- data[[variance_column]]
    fit <- fit_random_effects(data, tau_method)
    tau_label <- if (identical(tau_method, "PM")) "Paule-Mandel" else tau_method
    rows[[module_id]] <- model_row(fit, module_id, measure_label, tau_label, scope_label)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

run_leave_one_out <- function(study_effects) {
  rows <- list()
  index <- 1L
  for (module_id in MODULE_ORDER) {
    module_data <- study_effects[study_effects$module_id == module_id, , drop = FALSE]
    for (omitted in COHORT_ORDER) {
      data <- module_data[module_data$cohort_id != omitted, , drop = FALSE]
      data$effect <- data$smdh
      data$variance <- data$smdh_variance
      fit <- fit_random_effects(data, "REML")
      rows[[index]] <- model_row(
        fit,
        module_id,
        "SMDH (heteroscedastic standardized mean difference)",
        "REML",
        "leave-one-cohort-out exploratory sensitivity",
        omitted
      )
      index <- index + 1L
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

write_forest_plot <- function(study_effects, figure_dir) {
  make_plot <- function(device) {
    device()
    on.exit(dev.off(), add = TRUE)
    par(mfrow = c(2, 2), oma = c(5.2, 0.5, 3.2, 0.5), mar = c(3.7, 4.3, 2.5, 1.0))
    for (module_id in MODULE_ORDER) {
      data <- study_effects[study_effects$module_id == module_id, , drop = FALSE]
      data <- data[match(COHORT_ORDER, data$cohort_id), , drop = FALSE]
      data$effect <- data$smdh
      data$variance <- data$smdh_variance
      fit <- fit_random_effects(data, "REML")
      labels <- unname(COHORT_LABELS[data$cohort_id])
      forest(
        fit,
        slab = labels,
        rows = c(9, 7, 5, 3),
        ylim = c(-1, 11),
        xlim = c(-8.0, 8.0),
        alim = c(-4.0, 4.0),
        at = seq(-4, 4, by = 2),
        cex = 0.69,
        addfit = FALSE,
        header = FALSE,
        xlab = "Higher recorded severity minus lower recorded severity (standardized)",
        main = unname(MODULE_LABELS[[module_id]])
      )
      text(-4.0, 10.15, "Cohort", pos = 4, cex = 0.69, font = 2)
      text(4.0, 10.15, "SMDH [95% CI]", pos = 2, cex = 0.69, font = 2)
      addpoly(
        fit,
        row = 1,
        mlab = "RE pooled",
        cex = 0.69
      )
      abline(v = 0, lty = 3, col = "#4D4D4D")
    }
    mtext(
      "Supplementary Figure S3. Exploratory cohort-level random-effects meta-analysis of NP module scores.",
      outer = TRUE,
      side = 3,
      line = 1.2,
      cex = 0.95,
      font = 2
    )
    mtext(
      "Primary standardized effect: SMDH; REML tau^2; Knapp-Hartung intervals; k = 4 per module.",
      outer = TRUE,
      side = 1,
      line = 3.3,
      cex = 0.60
    )
    mtext(
      "All cohorts remain non-confirmatory; pooling does not remove age, source, platform, or presumed-key limitations.",
      outer = TRUE,
      side = 1,
      line = 2.1,
      cex = 0.55
    )
  }
  make_plot(function() pdf(file.path(figure_dir, "supplementary_figure_s3_np_exploratory_random_effects_meta_analysis.pdf"), width = 8.6, height = 8.2))
  make_plot(function() png(file.path(figure_dir, "supplementary_figure_s3_np_exploratory_random_effects_meta_analysis.png"), width = 5160, height = 4920, res = 600))
}

write_readme <- function(path, primary_results) {
  lines <- c(
    "# Exploratory NP Random-Effects Meta-Analysis",
    "",
    "This package adds a cohort-level exploratory quantitative synthesis to the existing IVDD project.",
    "The pre-existing donor/library-level Welch mean differences remain the primary analysis.",
    "",
    "## Scope and boundary",
    "",
    "- Four default NP cohorts are included: GSE230809 parent project, GSE244889, GSE153066, and GSE165722.",
    "- GSE229711 and GSE230808 remain one GSE230809 parent project, never two independent studies.",
    "- GSE251686 and GSM7986002 remain outside this primary meta-analysis package.",
    "- One effect is used per cohort and module; AF is not pooled with NP.",
    "- All included cohorts retain confirmatory_eligible=false.",
    "",
    "## Methods",
    "",
    "The primary standardized effect is metafor SMDH, which allows unequal group variances.",
    "A conventional pooled-SD Hedges g analysis and Paule-Mandel tau-squared estimate are sensitivity analyses.",
    "All models use REML or Paule-Mandel random effects with Knapp-Hartung intervals.",
    "REML Fisher scoring uses maxiter=10000 with all other metafor control values left at default; this changes only the iteration ceiling, not the estimator or interval method.",
    "With k=4 per module, heterogeneity estimates, Q tests, prediction intervals, and leave-one-cohort-out results are descriptive and unstable.",
    "No meta-regression, funnel plot, Egger test, causal claim, validation claim, or therapeutic inference is permitted.",
    "",
    "## Primary SMDH results"
  )
  for (index in seq_len(nrow(primary_results))) {
    row <- primary_results[index, ]
    lines <- c(
      lines,
      sprintf(
        "- %s: pooled SMDH %.4f, 95%% CI [%.4f, %.4f], prediction interval [%.4f, %.4f], I2 %.1f%%.",
        row$module,
        row$pooled_standardized_mean_difference,
        row$ci_lower,
        row$ci_upper,
        row$prediction_interval_lower,
        row$prediction_interval_upper,
        row$I_squared_percent
      )
    )
  }
  lines <- c(
    lines,
    "",
    "These results quantify the current cross-cohort pattern; they do not establish a universal IVDD program or a confirmatory biological conclusion."
  )
  writeLines(lines, path, useBytes = TRUE)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  # Keep caller-provided relative paths. Under the installed R 4.4.1 build,
  # normalizePath() can make absolute paths under this Unicode workspace fail
  # file.exists(), while the original relative paths are valid.
  summary_dir <- args[["summary-dir"]]
  table_dir <- args[["table-dir"]]
  figure_dir <- args[["figure-dir"]]
  output_dir <- args[["output-dir"]]
  if (!dir.exists(summary_dir)) {
    stop("The summary directory does not exist: ", summary_dir, call. = FALSE)
  }
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  group_file <- file.path(summary_dir, "donor_module_group_descriptives.csv")
  effect_file <- file.path(summary_dir, "donor_module_effects.csv")
  if (!file.exists(group_file) || !file.exists(effect_file)) {
    stop(
      "The current donor-level summary directory must contain group descriptives and effects. ",
      "group_file=", group_file, " (exists=", file.exists(group_file), "); ",
      "effect_file=", effect_file, " (exists=", file.exists(effect_file), ")",
      call. = FALSE
    )
  }

  wide <- make_wide_groups(group_file, effect_file)
  hedges <- calculate_effects(wide, "SMD")
  smdh <- calculate_effects(wide, "SMDH")
  study_effects <- cbind(
    wide,
    hedges_g = hedges$effect,
    hedges_g_variance = hedges$variance,
    hedges_g_standard_error = hedges$standard_error,
    smdh = smdh$effect,
    smdh_variance = smdh$variance,
    smdh_standard_error = smdh$standard_error
  )
  study_effects$cohort <- unname(COHORT_LABELS[study_effects$cohort_id])
  study_effects$module <- unname(MODULE_LABELS[study_effects$module_id])
  study_effects$all_cohorts_confirmatory_eligible <- FALSE
  study_effects$analysis_boundary <- "Exploratory cohort-level standardization only; cells remain nested observations; not confirmation, replication, mechanism, biomarker, or therapy evidence."

  primary_results <- run_models(
    study_effects,
    "smdh",
    "smdh_variance",
    "SMDH (heteroscedastic standardized mean difference)",
    "REML",
    "primary exploratory random-effects meta-analysis"
  )
  primary_results$hksj_p_value_BH_four_modules <- p.adjust(primary_results$hksj_p_value, method = "BH")
  hedges_sensitivity <- run_models(
    study_effects,
    "hedges_g",
    "hedges_g_variance",
    "Hedges g (pooled-SD standardized mean difference)",
    "REML",
    "pooled-SD standardized-effect sensitivity"
  )
  pm_sensitivity <- run_models(
    study_effects,
    "smdh",
    "smdh_variance",
    "SMDH (heteroscedastic standardized mean difference)",
    "PM",
    "tau-squared estimator sensitivity"
  )
  leave_one_out <- run_leave_one_out(study_effects)

  study_path <- file.path(table_dir, "supplementary_table_s7a_np_meta_analysis_study_effects.csv")
  primary_path <- file.path(table_dir, "supplementary_table_s7b_np_meta_analysis_primary_results.csv")
  sensitivity_path <- file.path(table_dir, "supplementary_table_s7c_np_meta_analysis_model_sensitivity.csv")
  leave_one_out_path <- file.path(table_dir, "supplementary_table_s7d_np_meta_analysis_leave_one_cohort_out.csv")
  params_path <- file.path(output_dir, "run_parameters.csv")
  readme_path <- file.path(output_dir, "README.md")

  write.csv(study_effects, study_path, row.names = FALSE, na = "")
  write.csv(primary_results, primary_path, row.names = FALSE, na = "")
  write.csv(rbind(hedges_sensitivity, pm_sensitivity), sensitivity_path, row.names = FALSE, na = "")
  write.csv(leave_one_out, leave_one_out_path, row.names = FALSE, na = "")
  write.csv(
    data.frame(
      parameter = c(
        "scope",
        "primary_effect_measure",
        "primary_tau_squared_estimator",
        "primary_interval_method",
        "primary_cohort_count_per_module",
        "excluded_cohort",
        "unit_of_analysis",
        "REML_convergence_control",
        "pooled_SD_sensitivity",
        "tau_squared_sensitivity",
        "multiplicity_adjustment",
        "prohibited_analyses"
      ),
      value = c(
        "four default NP cohort/module contrasts only",
        "SMDH via metafor::escalc(measure = 'SMDH')",
        "REML",
        "Knapp-Hartung",
        "4",
        "GSE251686; GSM7986002 remains permanently excluded",
        "donor or presumed donor/sample/library key; cells nested",
        "metafor REML maxiter=10000 with all other control values default; numerical iteration ceiling only, with the control used recorded per model",
        "Hedges g via metafor::escalc(measure = 'SMD')",
        "Paule-Mandel",
        "Benjamini-Hochberg across four primary module-level HKSJ p-values",
        "meta-regression; funnel plot; Egger test; pooled AF/NP analysis; confirmation or causal inference"
      ),
      stringsAsFactors = FALSE
    ),
    params_path,
    row.names = FALSE
  )
  write_readme(readme_path, primary_results)
  write_forest_plot(study_effects, figure_dir)

  output_paths <- c(
    study_path,
    primary_path,
    sensitivity_path,
    leave_one_out_path,
    params_path,
    readme_path,
    file.path(figure_dir, "supplementary_figure_s3_np_exploratory_random_effects_meta_analysis.pdf"),
    file.path(figure_dir, "supplementary_figure_s3_np_exploratory_random_effects_meta_analysis.png")
  )
  if (any(!file.exists(output_paths))) {
    stop("A declared meta-analysis output was not generated.", call. = FALSE)
  }
  generated_hashes <- setNames(vapply(output_paths, sha256, character(1)), manifest_path(output_paths))
  input_paths <- c(group_file, effect_file, script_path())
  input_paths <- input_paths[!is.na(input_paths)]
  input_hashes <- setNames(vapply(input_paths, sha256, character(1)), manifest_path(input_paths))
  manifest <- list(
    schema_version = 1L,
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    purpose = "Exploratory cohort-level random-effects meta-analysis of four default human NP score cohorts.",
    scope = list(
      cohorts = unname(COHORT_ORDER),
      modules = unname(MODULE_ORDER),
      excluded = "GSE251686 is excluded; GSM7986002 remains permanently excluded.",
      primary_analysis = "SMDH, REML, Knapp-Hartung; all included cohorts are non-confirmatory.",
      primary_result_table = gsub("\\\\", "/", primary_path)
    ),
    software = list(
      R = R.version.string,
      metafor = as.character(packageVersion("metafor")),
      digest = as.character(packageVersion("digest")),
      jsonlite = as.character(packageVersion("jsonlite"))
    ),
    input_sha256 = as.list(input_hashes),
    generated_artifact_sha256 = as.list(generated_hashes)
  )
  jsonlite::write_json(
    manifest,
    file.path(output_dir, "meta_analysis_manifest.json"),
    pretty = TRUE,
    auto_unbox = TRUE,
    na = "null"
  )

  cat("Exploratory NP random-effects meta-analysis completed.\n")
  cat("Primary results: ", primary_path, "\n", sep = "")
  cat("Forest plots: ", figure_dir, "\n", sep = "")
  cat("Manifest: ", file.path(output_dir, "meta_analysis_manifest.json"), "\n", sep = "")
}

main()
