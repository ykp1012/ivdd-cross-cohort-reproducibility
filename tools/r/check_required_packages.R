#!/usr/bin/env Rscript

requirements_path <- file.path("tools", "r", "requirements.tsv")
if (!file.exists(requirements_path)) {
  stop("Run this command from the repository root.", call. = FALSE)
}

requirements <- read.delim(
  requirements_path,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

runtime_row <- requirements[requirements$package == "R", , drop = FALSE]
if (nrow(runtime_row) != 1L) {
  stop("requirements.tsv must contain exactly one R runtime row.", call. = FALSE)
}

issues <- character()
if (getRversion() != runtime_row$version[[1L]]) {
  issues <- c(
    issues,
    sprintf("R %s is installed; expected R %s.", getRversion(), runtime_row$version[[1L]])
  )
}

package_rows <- requirements[requirements$package != "R", , drop = FALSE]
for (index in seq_len(nrow(package_rows))) {
  package_name <- package_rows$package[[index]]
  expected_version <- package_rows$version[[index]]
  if (!requireNamespace(package_name, quietly = TRUE)) {
    issues <- c(issues, sprintf("%s is not installed; expected %s.", package_name, expected_version))
    next
  }

  installed_version <- as.character(utils::packageVersion(package_name))
  if (installed_version != expected_version) {
    issues <- c(
      issues,
      sprintf("%s %s is installed; expected %s.", package_name, installed_version, expected_version)
    )
  }
}

if (length(issues) > 0L) {
  writeLines(c("R environment check failed:", paste0("- ", issues)))
  quit(status = 1L)
}

writeLines("R environment check passed.")
