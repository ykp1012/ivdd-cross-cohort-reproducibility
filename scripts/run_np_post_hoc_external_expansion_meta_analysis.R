#!/usr/bin/env Rscript

# Run an isolated, post hoc, source-family-aware random-effects synthesis of human NP
# module scores. This script never overwrites the frozen four-cohort S7/S3
# package. All cohorts remain exploratory, and accession-level independence is
# not interpreted as patient-level confirmation.

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

SIX_COHORT_ORDER <- c(
  "GSE230809_discovery",
  "GSE244889_directional",
  "GSE153066_support",
  "GSE165722_score_level",
  "GSE186542_external_count_support",
  "GSE167931_external_fpkm_support"
)

REPLACEMENT_COHORT_ORDER <- c(
  "GSE230809_discovery",
  "GSE244889_directional",
  "GSE153066_support",
  "GSE165722_score_level",
  "GSE186542_external_count_support",
  "GSE245147_external_native_comparison_support"
)

COHORT_LABELS <- c(
  GSE230809_discovery = "GSE230809 parent (3 vs 8)",
  GSE244889_directional = "GSE244889 (NP 4 vs 3)",
  GSE153066_support = "GSE153066 (NP 8 vs 8)",
  GSE165722_score_level = "GSE165722 (NP 4 vs 4)",
  GSE186542_external_count_support = "GSE186542 (3 vs 3)",
  GSE167931_external_fpkm_support = "GSE167931 FPKM (4 vs 5)",
  GSE245147_external_native_comparison_support = "GSE245147 native (3 vs 3)"
)

usage <- function() {
  paste(
    "Usage:",
    "Rscript scripts/run_np_post_hoc_external_expansion_meta_analysis.R",
    "--summary-dir data/derived/donor_module_effect_summary_external_expansion",
    "--contrast-spec config/np_post_hoc_external_expansion_contrast_spec.csv",
    "--table-dir results/supplementary_tables",
    "--figure-dir results/supplementary_figures",
    "--output-dir data/derived/np_post_hoc_external_expansion_meta_analysis",
    "[--analysis-set six|replace167931]",
    sep = "\n"
  )
}

parse_args <- function(args) {
  expected <- c("summary-dir", "contrast-spec", "table-dir", "figure-dir", "output-dir", "analysis-set")
  values <- setNames(rep(NA_character_, length(expected)), expected)
  values[["analysis-set"]] <- "six"
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
  required <- c("summary-dir", "contrast-spec", "table-dir", "figure-dir", "output-dir")
  if (anyNA(values[required]) || any(!nzchar(values[required]))) {
    stop("The five input/output options are required.\n", usage(), call. = FALSE)
  }
  if (!values[["analysis-set"]] %in% c("six", "replace167931")) {
    stop("--analysis-set must be 'six' or 'replace167931'.", call. = FALSE)
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

read_cohort_order <- function(config_path, analysis_set) {
  config <- read.csv(config_path, stringsAsFactors = FALSE, check.names = FALSE)
  assert_columns(config, c("cohort_id", "compartment", "confirmatory_eligible"), basename(config_path))
  config <- config[config$compartment == "NP", , drop = FALSE]
  expected_order <- if (identical(analysis_set, "six")) SIX_COHORT_ORDER else REPLACEMENT_COHORT_ORDER
  if (!identical(config$cohort_id, expected_order)) {
    stop(
      "The post hoc expansion cohort configuration does not match --analysis-set. Found: ",
      paste(config$cohort_id, collapse = "; "),
      call. = FALSE
    )
  }
  if (any(tolower(config$confirmatory_eligible) != "false")) {
    stop("All expansion cohorts must remain non-confirmatory.", call. = FALSE)
  }
  config$cohort_id
}

make_wide_groups <- function(group_file, effect_file, cohort_order) {
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
    groups$cohort_id %in% cohort_order &
      groups$compartment == "NP" &
      groups$module_id %in% MODULE_ORDER,
    ,
    drop = FALSE
  ]
  effects <- effects[
    effects$cohort_id %in% cohort_order &
      effects$compartment == "NP" &
      effects$module_id %in% MODULE_ORDER,
    ,
    drop = FALSE
  ]

  expected_groups <- length(cohort_order) * length(MODULE_ORDER) * 2L
  expected_effects <- length(cohort_order) * length(MODULE_ORDER)
  if (nrow(groups) != expected_groups) {
    stop("Expected ", expected_groups, " expansion NP group rows, found ", nrow(groups), call. = FALSE)
  }
  if (nrow(effects) != expected_effects) {
    stop("Expected ", expected_effects, " expansion NP effect rows, found ", nrow(effects), call. = FALSE)
  }
  if (!all(tolower(effects$confirmatory_eligible) == "false")) {
    stop("The project boundary requires every pooled cohort to remain non-confirmatory.", call. = FALSE)
  }

  groups <- as_number(
    groups,
    c("n_donor_or_library_keys", "mean_module_score_log1p_cpm", "sd_module_score_log1p_cpm"),
    basename(group_file)
  )
  if (any(groups$n_donor_or_library_keys < 2L) || any(groups$sd_module_score_log1p_cpm <= 0)) {
    stop("Every expansion meta-analysis arm must have n >= 2 and a positive sample SD.", call. = FALSE)
  }

  key_columns <- c("cohort_id", "module_id", "contrast_label")
  target <- groups[groups$contrast_arm == "target", c(key_columns, "severity_group", "n_donor_or_library_keys", "mean_module_score_log1p_cpm", "sd_module_score_log1p_cpm")]
  comparison <- groups[groups$contrast_arm == "comparison", c(key_columns, "severity_group", "n_donor_or_library_keys", "mean_module_score_log1p_cpm", "sd_module_score_log1p_cpm")]
  if (nrow(target) != expected_effects || nrow(comparison) != expected_effects) {
    stop("Every expansion NP cohort/module contrast must have exactly one target and one comparison arm.", call. = FALSE)
  }

  names(target)[(length(key_columns) + 1L):ncol(target)] <- c("target_group", "target_n", "target_mean", "target_sd")
  names(comparison)[(length(key_columns) + 1L):ncol(comparison)] <- c("comparison_group", "comparison_n", "comparison_mean", "comparison_sd")
  wide <- merge(target, comparison, by = key_columns, all = TRUE, sort = FALSE)
  if (nrow(wide) != expected_effects || any(!complete.cases(wide))) {
    stop("Failed to form complete target/comparison arms for every expansion NP contrast.", call. = FALSE)
  }

  effect_meta <- effects[, c("cohort_id", "module_id", "cohort_limitations")]
  wide <- merge(wide, effect_meta, by = c("cohort_id", "module_id"), all.x = TRUE, sort = FALSE)
  if (any(is.na(wide$cohort_limitations))) {
    stop("An expansion NP contrast lacks its cohort limitation text.", call. = FALSE)
  }

  wide$cohort_id <- factor(wide$cohort_id, levels = cohort_order)
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

run_leave_one_out <- function(study_effects, cohort_order, scope_label) {
  rows <- list()
  index <- 1L
  for (module_id in MODULE_ORDER) {
    module_data <- study_effects[study_effects$module_id == module_id, , drop = FALSE]
    for (omitted in cohort_order) {
      data <- module_data[module_data$cohort_id != omitted, , drop = FALSE]
      data$effect <- data$smdh
      data$variance <- data$smdh_variance
      fit <- fit_random_effects(data, "REML")
      rows[[index]] <- model_row(
        fit,
        module_id,
        "SMDH (heteroscedastic standardized mean difference)",
        "REML",
        scope_label,
        omitted
      )
      index <- index + 1L
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

draw_forest_panel <- function(data, fit, module_id, cohort_order) {
  data <- data[match(cohort_order, data$cohort_id), , drop = FALSE]
  study_lower <- data$smdh - qnorm(0.975) * data$smdh_standard_error
  study_upper <- data$smdh + qnorm(0.975) * data$smdh_standard_error
  limit <- ceiling(max(abs(c(study_lower, study_upper, fit$ci.lb, fit$ci.ub, 0))) + 0.3)
  limit <- max(limit, 3)
  rows <- rev(seq_len(nrow(data))) + 1L

  par(xpd = FALSE)
  plot(
    NA,
    xlim = c(-limit, limit),
    ylim = c(0.3, nrow(data) + 2.3),
    xlab = "Higher recorded severity minus lower recorded severity (SMDH)",
    ylab = "",
    yaxt = "n",
    bty = "n",
    main = unname(MODULE_LABELS[[module_id]])
  )
  axis(1, at = pretty(c(-limit, limit), n = 5))
  abline(v = 0, lty = 3, col = "#4D4D4D")
  par(xpd = NA)
  label_left <- -limit - (0.04 * limit)
  label_right <- limit + (0.04 * limit)
  text(label_left, nrow(data) + 1.7, "Cohort", pos = 2, font = 2, cex = 0.60)
  text(label_right, nrow(data) + 1.7, "SMDH [95% CI]", pos = 4, font = 2, cex = 0.60)
  for (index in seq_len(nrow(data))) {
    segments(study_lower[index], rows[index], study_upper[index], rows[index], lwd = 1.1, col = "#1A1A1A")
    points(data$smdh[index], rows[index], pch = 15, cex = 0.72, col = "#0072B2")
    text(label_left, rows[index], unname(COHORT_LABELS[[data$cohort_id[index]]]), pos = 2, cex = 0.52)
    text(
      label_right,
      rows[index],
      sprintf("%.2f [%.2f, %.2f]", data$smdh[index], study_lower[index], study_upper[index]),
      pos = 4,
      cex = 0.50
    )
  }
  diamond_y <- 1.1
  polygon(
    x = c(fit$ci.lb, as.numeric(fit$b[1L]), fit$ci.ub, as.numeric(fit$b[1L])),
    y = c(diamond_y, diamond_y + 0.22, diamond_y, diamond_y - 0.22),
    border = "#D55E00",
    col = "#D55E00"
  )
  text(label_left, diamond_y, "RE pooled (REML, Knapp-Hartung)", pos = 2, font = 2, cex = 0.53)
  text(label_right, diamond_y, sprintf("%.2f [%.2f, %.2f]", fit$b[1L], fit$ci.lb, fit$ci.ub), pos = 4, font = 2, cex = 0.51)
  par(xpd = FALSE)
}

write_forest_plot <- function(study_effects, figure_dir, cohort_order, analysis_set) {
  if (identical(analysis_set, "six")) {
    figure_basename <- "supplementary_figure_s4_np_post_hoc_external_expansion_meta_analysis"
    figure_title <- "Supplementary Figure S4. Post hoc exploratory six-cohort random-effects meta-analysis of human NP module scores."
    figure_note <- "GSE186542/GSE167931 are accession-level additions only; all cohorts remain non-confirmatory score-level analyses."
  } else {
    figure_basename <- "supplementary_figure_s5_np_source_family_replacement_sensitivity"
    figure_title <- "Supplementary Figure S5. Source-family replacement sensitivity for the post hoc NP random-effects synthesis."
    figure_note <- "GSE245147 replaces GSE167931; the two same-family cohorts are never pooled as independent studies."
  }
  k_label <- length(cohort_order)
  make_plot <- function(device) {
    device()
    on.exit(dev.off(), add = TRUE)
    par(mfrow = c(2, 2), oma = c(5.0, 0.5, 3.2, 0.5), mar = c(3.9, 8.5, 2.3, 8.1), cex.main = 0.83)
    for (module_id in MODULE_ORDER) {
      data <- study_effects[study_effects$module_id == module_id, , drop = FALSE]
      data$effect <- data$smdh
      data$variance <- data$smdh_variance
      fit <- fit_random_effects(data, "REML")
      draw_forest_panel(data, fit, module_id, cohort_order)
    }
    mtext(
      figure_title,
      outer = TRUE,
      side = 3,
      line = 1.2,
      cex = 0.94,
      font = 2
    )
    mtext(
      sprintf("Primary effect: SMDH; REML tau^2; Knapp-Hartung intervals; k = %d per module.", k_label),
      outer = TRUE,
      side = 1,
      line = 3.2,
      cex = 0.63
    )
    mtext(
      figure_note,
      outer = TRUE,
      side = 1,
      line = 2.0,
      cex = 0.54
    )
  }
  make_plot(function() pdf(file.path(figure_dir, paste0(figure_basename, ".pdf")), width = 11, height = 10))
  make_plot(function() png(file.path(figure_dir, paste0(figure_basename, ".png")), width = 6600, height = 6000, res = 600))
}

write_readme <- function(path, primary_results, cohort_order, analysis_set) {
  if (identical(analysis_set, "six")) {
    heading <- "# Post Hoc External Expansion: Exploratory NP Random-Effects Meta-Analysis"
    opening <- "It adds GSE186542 and the GSE167931 FPKM representation after candidate GEO auditing."
    boundary <- c(
      "- Six human NP cohort/module contrasts are included once each; k = 6 per module.",
      "- GSE167931 TPM is a same-sample processing sensitivity, not an additional study.",
      "- GSE186542 and GSE167931 are accession/BioProject-level additions; patient-level overlap cannot be excluded from public metadata."
    )
    conclusion <- "The expansion quantifies how two newly audited public datasets shift the cross-cohort pattern. It does not establish a universal IVDD program or an independent patient-level validation."
  } else {
    heading <- "# Source-Family Replacement Sensitivity: Exploratory NP Random-Effects Meta-Analysis"
    opening <- "It replaces the GSE167931 FPKM representation with the native clinical-comparison subset of GSE245147."
    boundary <- c(
      "- Six human NP cohort/module contrasts are included once each; k = 6 per module.",
      "- GSE245147 includes only Degenerated_1-3 versus NO_Degenerated_1-3; P2/P8 and DMSO/H-151 arms are excluded.",
      "- GSE167931 and GSE245147 are not pooled together because their source lab/author family overlaps and patient-level reuse cannot be excluded."
    )
    conclusion <- "This is a source-family replacement sensitivity, not an independent validation analysis. It tests whether the pooled pattern depends on the choice of one potentially related external cohort."
  }
  lines <- c(
    heading,
    "",
    "This package is separate from the frozen four-cohort S7/S3 analysis.",
    opening,
    "",
    "## Scope and boundary",
    "",
    "- GSE230809 remains one parent project and is not split into GSE229711 and GSE230808 studies.",
    "- GSE251686 and GSM7986002 remain outside this package.",
    boundary,
    "- All included contrasts remain post hoc, score-level, and confirmatory_eligible=false.",
    "",
    "## Methods",
    "",
    "The primary standardized effect is metafor SMDH, which permits unequal group variances.",
    "A conventional pooled-SD Hedges g analysis and Paule-Mandel tau-squared estimate are sensitivity analyses.",
    "All models use random effects with Knapp-Hartung intervals.",
    "REML Fisher scoring uses maxiter=10000 with all other metafor control values left at default; this changes only the iteration ceiling, not the estimator or interval method.",
    sprintf("With k=%d per module and heterogeneous source material, prediction intervals, Q tests, I2, and leave-one-cohort-out results are descriptive and unstable.", length(cohort_order)),
    "No meta-regression, funnel plot, Egger test, causal claim, patient-level replication claim, biomarker claim, or therapeutic inference is permitted.",
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
    conclusion
  )
  writeLines(lines, path, useBytes = TRUE)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  # Preserve caller-provided relative paths in this Unicode workspace because
  # normalizing to an absolute path is unreliable in the installed R build.
  summary_dir <- args[["summary-dir"]]
  contrast_spec <- args[["contrast-spec"]]
  table_dir <- args[["table-dir"]]
  figure_dir <- args[["figure-dir"]]
  output_dir <- args[["output-dir"]]
  analysis_set <- args[["analysis-set"]]
  for (input_path in c(summary_dir, contrast_spec)) {
    if (!file.exists(input_path) && !dir.exists(input_path)) {
      stop("Required input does not exist: ", input_path, call. = FALSE)
    }
  }
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cohort_order <- read_cohort_order(contrast_spec, analysis_set)
  group_file <- file.path(summary_dir, "donor_module_group_descriptives.csv")
  effect_file <- file.path(summary_dir, "donor_module_effects.csv")
  if (!file.exists(group_file) || !file.exists(effect_file)) {
    stop("The expansion summary directory must contain group descriptives and effects.", call. = FALSE)
  }

  wide <- make_wide_groups(group_file, effect_file, cohort_order)
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
  study_effects$analysis_boundary <- "Post hoc exploratory accession-level external expansion; score-level standardization only; not patient-level confirmation, replication, mechanism, biomarker, or therapy evidence."

  scope_label <- if (identical(analysis_set, "six")) {
    "post hoc exploratory six-cohort random-effects meta-analysis"
  } else {
    "source-family replacement sensitivity: six-cohort random-effects meta-analysis"
  }
  primary_results <- run_models(
    study_effects,
    "smdh",
    "smdh_variance",
    "SMDH (heteroscedastic standardized mean difference)",
    "REML",
    scope_label
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
  loo_scope_label <- if (identical(analysis_set, "six")) {
    "leave-one-cohort-out post hoc exploratory sensitivity"
  } else {
    "leave-one-cohort-out source-family replacement sensitivity"
  }
  leave_one_out <- run_leave_one_out(study_effects, cohort_order, loo_scope_label)

  if (identical(analysis_set, "six")) {
    table_stem <- "np_post_hoc_external_expansion"
    table_prefix <- "s8"
  } else {
    table_stem <- "np_source_family_replacement_sensitivity"
    table_prefix <- "s9"
  }
  study_path <- file.path(table_dir, sprintf("supplementary_table_%sa_%s_study_effects.csv", table_prefix, table_stem))
  primary_path <- file.path(table_dir, sprintf("supplementary_table_%sb_%s_primary_results.csv", table_prefix, table_stem))
  sensitivity_path <- file.path(table_dir, sprintf("supplementary_table_%sc_%s_model_sensitivity.csv", table_prefix, table_stem))
  leave_one_out_path <- file.path(table_dir, sprintf("supplementary_table_%sd_%s_leave_one_cohort_out.csv", table_prefix, table_stem))
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
        "newly_added_external_cohorts",
        "independence_boundary",
        "unit_of_analysis",
        "paired_processing_representation_excluded",
        "excluded_cohort",
        "REML_convergence_control",
        "pooled_SD_sensitivity",
        "tau_squared_sensitivity",
        "multiplicity_adjustment",
        "prohibited_analyses"
      ),
      value = c(
        if (identical(analysis_set, "six")) "post hoc six-cohort human NP expansion; frozen four-cohort S7/S3 package unchanged" else "source-family replacement sensitivity; GSE245147 replaces GSE167931 and frozen four-cohort S7/S3 package remains unchanged",
        "SMDH via metafor::escalc(measure = 'SMDH')",
        "REML",
        "Knapp-Hartung",
        as.character(length(cohort_order)),
        if (identical(analysis_set, "six")) "GSE186542 external count-score support; GSE167931 FPKM external processed-score support" else "GSE186542 external count-score support; GSE245147 native clinical-comparison RPKM score support",
        if (identical(analysis_set, "six")) "BioProject/accession-level additions only; patient-level independence remains unverified from public metadata" else "GSE167931 and GSE245147 are source-family related; do not interpret replacement as independent validation",
        "donor or presumed donor/sample/library key; cells remain nested",
        if (identical(analysis_set, "six")) "GSE167931 TPM paired representation" else "GSE167931 FPKM is excluded by replacement design; P2/P8 and DMSO/H-151 GSE245147 arms are excluded",
        "GSE251686; GSM7986002 remains permanently excluded",
        "metafor REML maxiter=10000 with all other control values default; numerical iteration ceiling only, with the control recorded per model",
        "Hedges g via metafor::escalc(measure = 'SMD')",
        "Paule-Mandel",
        "Benjamini-Hochberg across four primary module-level HKSJ p-values",
        "meta-regression; funnel plot; Egger test; pooled AF/NP analysis; patient-level confirmation; causal or therapeutic inference"
      ),
      stringsAsFactors = FALSE
    ),
    params_path,
    row.names = FALSE
  )
  write_readme(readme_path, primary_results, cohort_order, analysis_set)
  write_forest_plot(study_effects, figure_dir, cohort_order, analysis_set)

  if (identical(analysis_set, "six")) {
    figure_basename <- "supplementary_figure_s4_np_post_hoc_external_expansion_meta_analysis"
  } else {
    figure_basename <- "supplementary_figure_s5_np_source_family_replacement_sensitivity"
  }

  output_paths <- c(
    study_path,
    primary_path,
    sensitivity_path,
    leave_one_out_path,
    params_path,
    readme_path,
    file.path(figure_dir, paste0(figure_basename, ".pdf")),
    file.path(figure_dir, paste0(figure_basename, ".png"))
  )
  if (any(!file.exists(output_paths))) {
    stop("A declared post hoc expansion output was not generated.", call. = FALSE)
  }
  generated_hashes <- setNames(vapply(output_paths, sha256, character(1)), manifest_path(output_paths))
  input_paths <- c(group_file, effect_file, contrast_spec, script_path())
  input_paths <- input_paths[!is.na(input_paths)]
  input_hashes <- setNames(vapply(input_paths, sha256, character(1)), manifest_path(input_paths))
  manifest <- list(
    schema_version = 1L,
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    purpose = if (identical(analysis_set, "six")) "Post hoc exploratory six-cohort random-effects synthesis of human NP module scores." else "Source-family replacement sensitivity for a six-cohort exploratory random-effects synthesis of human NP module scores.",
    scope = list(
      cohorts = unname(cohort_order),
      modules = unname(MODULE_ORDER),
      newly_added = if (identical(analysis_set, "six")) c("GSE186542", "GSE167931 FPKM") else c("GSE186542", "GSE245147 native subset"),
      excluded = if (identical(analysis_set, "six")) "GSE167931 TPM is a paired processing sensitivity; GSE251686 and GSM7986002 remain excluded." else "GSE167931 is excluded by source-family replacement; GSE245147 passage/treatment arms, GSE251686, and GSM7986002 remain excluded.",
      independence_boundary = if (identical(analysis_set, "six")) "GSE186542 and GSE167931 are accession/BioProject-level additions; patient-level independence is not established." else "GSE245147 is a source-family replacement for GSE167931; patient-level independence from that family is not established.",
      primary_analysis = if (identical(analysis_set, "six")) "SMDH, REML, Knapp-Hartung; all included cohorts are post hoc and non-confirmatory." else "SMDH, REML, Knapp-Hartung; replacement is sensitivity-only and non-confirmatory.",
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

  cat("Post hoc external-expansion NP random-effects meta-analysis completed.\n")
  cat("Primary results: ", primary_path, "\n", sep = "")
  cat("Forest plots: ", figure_dir, "\n", sep = "")
  cat("Manifest: ", file.path(output_dir, "meta_analysis_manifest.json"), "\n", sep = "")
}

main()
