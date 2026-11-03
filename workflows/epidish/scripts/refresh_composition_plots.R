#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: refresh_composition_plots.R COMPLETED_RESULT_DIR NEW_OUTPUT_DIR", call. = FALSE)
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
source(file.path(dirname(script_path), "epidish_common.R"))

result_dir <- assert_dir(args[[1L]], "completed comparison result")
completion_path <- assert_file(file.path(result_dir, "EPIDISHR_COMPLETED"))
manifest_path <- assert_file(file.path(result_dir, "provenance", "run_manifest.json"))
input_manifest_path <- assert_file(file.path(result_dir, "provenance", "input_manifest.json"))
summary_dir <- comparison_summary_dir(result_dir)
commands_path <- assert_file(file.path(result_dir, "provenance", "commands.log"))

atomic_ggsave <- function(path, plot, width, height) {
  temporary <- paste0(path, ".tmp.pdf")
  on.exit(unlink(temporary), add = TRUE)
  ggplot2::ggsave(temporary, plot, width = width, height = height)
  if (!file.rename(temporary, path)) stopf("Could not publish refreshed plot: %s", path)
  invisible(path)
}

atomic_json <- function(value, path) {
  temporary <- paste0(path, ".tmp")
  on.exit(unlink(temporary), add = TRUE)
  jsonlite::write_json(
    value, temporary, auto_unbox = TRUE, pretty = TRUE, null = "null"
  )
  if (!file.rename(temporary, path)) stopf("Could not publish refreshed manifest: %s", path)
  invisible(path)
}

samples <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)$samples
if (!is.character(samples) || !length(samples)) stopf("Run manifest lacks sample order")

input_manifest <- jsonlite::read_json(input_manifest_path, simplifyVector = FALSE)
thresholds <- sort(as.integer(unlist(input_manifest$coverage_thresholds)))
output_dir <- args[[2L]]
reserve_output_bundle(output_dir, script_path, c(manifest_path, input_manifest_path))
started_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
plot_records <- list()
options(error = function() {
  try(finish_output_bundle(output_dir, script_path, FALSE), silent = TRUE)
  quit(save = "no", status = 1L, runLast = FALSE)
})

for (profile in c("full", "hematopoietic")) {
  for (threshold in thresholds) {
    comparison_dir <- assert_dir(file.path(
      result_dir, "comparisons", profile, threshold_label(threshold)
    ))
    table_files <- c(epidish_rpc = "epidish_rpc_fractions.tsv")
    if (file.exists(file.path(comparison_dir, "nnls_unweighted_fractions.tsv"))) {
      table_files <- c(
        table_files, nnls_unweighted = "nnls_unweighted_fractions.tsv"
      )
    }
    table_files <- c(
      table_files,
      uxm_weighted_matched = "uxm_weighted_matched_fractions.tsv",
      uxm_native_context = "uxm_native_fractions.tsv"
    )
    tables <- lapply(
      table_files,
      function(name) read_tsv(file.path(comparison_dir, name))
    )
    composition <- do.call(rbind, Map(
      function(table, method) fractions_long(table, profile, threshold, method),
      tables, names(tables)
    ))
    composition$sample <- factor(composition$sample, levels = samples)
    composition$method <- factor(composition$method, levels = names(tables))
    plot <- ggplot2::ggplot(
      composition, ggplot2::aes(x = sample, y = fraction, fill = cell_type)
    ) +
      ggplot2::geom_col(width = 0.82) +
      ggplot2::facet_wrap(~method, ncol = 1L) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        title = sprintf("%s, coverage >= %d", profile, threshold),
        x = NULL, y = "Estimated fraction", fill = "Cell type"
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    destination <- file.path(output_dir, "comparisons", profile, threshold_label(threshold))
    dir.create(destination, recursive = TRUE)
    atomic_ggsave(
      file.path(destination, "composition_comparison.pdf"),
      plot, width = 11, height = 10
    )

    plot_records[[length(plot_records) + 1L]] <- list(
      profile = profile, coverage = threshold,
      input_tables = lapply(file.path(comparison_dir, table_files), artifact_record),
      output = artifact_record(file.path(destination, "composition_comparison.pdf")))
  }
}

manifest <- list(started_utc = started_utc,
  completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  source_result = result_dir, source_manifest = artifact_record(manifest_path),
  script = artifact_record(script_path), plots = plot_records,
  purpose = "Replot unchanged saved estimates; source completed bundle remains immutable")
atomic_json(manifest, file.path(output_dir, "plot_manifest.json"))
finish_output_bundle(output_dir, script_path)
message("Composition plots exported: ", normalizePath(output_dir))
