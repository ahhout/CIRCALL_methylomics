#!/usr/bin/env Rscript
# Display saved cohort estimates together; never refit or pool marker masks.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) {
  stop("Usage: export_saved_deconvolution.R NEW_OUTPUT_DIR METADATA_TSV COVERAGE COMPLETED_EPIDISH_RUN...", call. = FALSE)
}
started_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]]))
source(file.path(dirname(script_path), "epidish_common.R"))
source(file.path(dirname(script_path), "plot_epidish_standalone.R"))
output_dir <- args[[1L]]
name <- basename(sub("/+$", "", output_dir))
stamp <- substr(name, 1L, 16L)
parsed <- as.POSIXct(strptime(stamp, "%Y-%m-%d_%H-%M", tz = "UTC"))
if (!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}_([01][0-9]|2[0-3])-[0-5][0-9]-[A-Za-z0-9][A-Za-z0-9._-]*$", name) ||
    is.na(parsed) || format(parsed, "%Y-%m-%d_%H-%M", tz = "UTC") != stamp) {
  stopf("Output name requires a valid UTC YYYY-MM-DD_HH-MM timestamp")
}
if (file.exists(output_dir) || (!is.na(Sys.readlink(output_dir)) && nzchar(Sys.readlink(output_dir)))) stopf("Output already exists: %s", output_dir)
metadata_path <- assert_file(args[[2L]])
metadata <- read_tsv(metadata_path)
if (!all(c("sample", "patient", "compartment", "blast_pct") %in% names(metadata)) ||
    !nrow(metadata) || anyDuplicated(metadata$sample)) stopf("Invalid expected sample metadata")
if (!grepl("^[1-9][0-9]*$", args[[3L]])) stopf("Coverage must be a positive integer")
coverage <- as.integer(args[[3L]])
sources <- vapply(args[-(1:3)], assert_dir, character(1))
if (anyDuplicated(sources)) stopf("Duplicate source runs")
fractions <- warnings <- native <- records <- settings <- list()
seen <- character()
saved_metadata_rows <- metadata_differences <- list()
for (i in seq_along(sources)) {
  run <- sources[[i]]
  assert_file(file.path(run, "EPIDISHR_COMPLETED"))
  input_path <- file.path(run, "provenance", "input_manifest.json")
  manifest_path <- file.path(run, "provenance", "run_manifest.json")
  input <- jsonlite::read_json(assert_file(input_path), simplifyVector = FALSE)
  manifest <- jsonlite::read_json(assert_file(manifest_path), simplifyVector = FALSE)
  if (sha256_file(input_path) != manifest$input_manifest_sha256) stopf("Input manifest checksum differs: %s", run)
  samples <- unlist(input$samples, use.names = FALSE)
  if (!length(samples) || anyDuplicated(samples) || length(intersect(seen, samples)) ||
      !setequal(samples, unlist(manifest$samples))) stopf("Invalid or overlapping source cohort: %s", run)
  seen <- c(seen, samples)
  cohort <- paste0("cohort-", i, "-", length(samples), "samples")
  summary_dir <- comparison_summary_dir(run)
  checked <- function(name) {
    path <- assert_file(file.path(summary_dir, name))
    matches <- Filter(function(r) identical(r$name, name), manifest$summary_files)
    if (length(matches) != 1L || sha256_file(path) != matches[[1L]]$sha256) {
      stopf("Saved summary checksum differs: %s", path)
    }
    path
  }
  fraction_path <- checked("fractions_long.tsv.gz")
  warning_path <- checked("epidish_warnings_all.tsv")
  x <- read_tsv(fraction_path)
  w <- read_tsv(warning_path)
  thresholds <- as.integer(unlist(input$coverage_thresholds))
  if (!coverage %in% thresholds || !setequal(unique(x$sample), samples)) stopf("Source lacks requested coverage or samples")
  methods <- c("epidish_rpc", "nnls_unweighted", "uxm_weighted_matched", "uxm_native_context")
  if (!setequal(unique(x$method), methods)) stopf("Unexpected source methods")
  keys <- x[, c("sample", "profile", "coverage_threshold", "method", "cell_type")]
  if (anyDuplicated(keys) || any(!is.finite(x$fraction)) || any(x$fraction < 0 | x$fraction > 1)) stopf("Invalid saved fractions")
  for (profile in c("full", "hematopoietic")) {
    for (threshold in thresholds) for (method in methods) {
      d <- x[x$profile == profile & x$coverage_threshold == threshold & x$method == method, ]
      cells <- if (profile == "full") 40L else 7L
      if (nrow(d) != length(samples) * cells || !setequal(unique(d$sample), samples) ||
          length(unique(d$cell_type)) != cells || any(table(d$sample) != cells) ||
          any(abs(tapply(d$fraction, d$sample, sum) - 1) > 1e-5)) stopf("Incomplete fraction grid: %s", run)
    }
    key <- if (profile == "full") "full_uxm_run" else "hematopoietic_uxm_run"
    native_path <- assert_file(file.path(input[[key]]$path, "results", "lineage_fractions.csv"))
    if (sha256_file(native_path) != input[[key]]$native_result_sha256) stopf("Native UXM checksum differs")
    wide <- read.csv(native_path, check.names = FALSE, stringsAsFactors = FALSE)
    if (!identical(names(wide), c("CellType", samples))) stopf("Native UXM sample order differs")
    validate_fraction_matrix(as.matrix(wide[, samples, drop = FALSE]), samples, wide$CellType,
                             "native UXM", cells_by_samples = TRUE)
    n <- fractions_long(wide, profile, 0L, "uxm_native_context")
    n$source_run <- normalizePath(input[[key]]$path)
    n$source_cohort <- cohort
    native[[length(native) + 1L]] <- n
    records[[length(records) + 1L]] <- artifact_record(native_path)
  }
  saved_metadata_path <- assert_file(file.path(run, "provenance", "sample_metadata.tsv"))
  saved_metadata <- read_tsv(saved_metadata_path)
  for (column in c("patient", "compartment", "blast_pct")) {
    a <- as.character(saved_metadata[[column]][match(samples, saved_metadata$sample)])
    b <- as.character(metadata[[column]][match(samples, metadata$sample)])
    if (column != "blast_pct" && !identical(a, b)) stopf("Saved metadata differs for %s", column)
    different <- which((is.na(a) != is.na(b)) | (!is.na(a) & !is.na(b) & a != b))
    if (length(different)) metadata_differences[[length(metadata_differences) + 1L]] <- data.frame(
      sample = samples[different], field = column, saved_value = a[different], current_value = b[different])
  }
  saved_metadata_rows[[i]] <- saved_metadata[match(samples, saved_metadata$sample), , drop = FALSE]
  x$source_run <- run
  x$source_cohort <- cohort
  x$rpc_maxit <- input$epidish_maxit
  w$source_run <- rep(run, nrow(w))
  w$source_cohort <- rep(cohort, nrow(w))
  fractions[[i]] <- x
  warnings[[i]] <- w
  settings[[i]] <- data.frame(source_cohort = cohort, source_run = run,
    sample_count = length(samples), coverage_thresholds = paste(thresholds, collapse = ","),
    rpc_maxit = input$epidish_maxit, epidish_version = input$epidish_version)
  for (path in c(input_path, manifest_path, fraction_path, warning_path, saved_metadata_path)) {
    records[[length(records) + 1L]] <- artifact_record(path)
  }
}
if (!setequal(seen, metadata$sample)) stopf("Source cohort union differs from expected metadata samples")
plot_metadata <- do.call(rbind, saved_metadata_rows)
plot_metadata <- plot_metadata[match(metadata$sample, plot_metadata$sample), , drop = FALSE]
all_fractions <- do.call(rbind, fractions)
all_warnings <- do.call(rbind, warnings)
all_native <- do.call(rbind, native)
source_settings <- do.call(rbind, settings)
reserve_output_bundle(output_dir, script_path, c(metadata_path, file.path(sources, "provenance", "run_manifest.json")))
manifest_path <- file.path(output_dir, "display_manifest.json")
provenance <- list(started_utc = started_utc, completed_utc = NULL, status = "running",
  purpose = "Combine saved estimates for display; no refitting, sample merging or new cohort-wide mask",
  samples = metadata$sample, plotted_coverage = coverage,
  sources = source_settings, inputs = records, metadata = artifact_record(metadata_path),
  scripts = lapply(c(script_path, file.path(dirname(script_path), "plot_epidish_standalone.R"),
                    file.path(dirname(script_path), "epidish_common.R")), artifact_record),
  command = commandArgs(), r_version = as.character(getRversion()))
jsonlite::write_json(provenance, manifest_path, pretty = TRUE, auto_unbox = TRUE, null = "null")
options(error = function() {
  try(finish_output_bundle(output_dir, script_path, FALSE), silent = TRUE)
  provenance$status <- "failed"
  provenance$completed_utc <- NULL
  provenance$ended_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  try(jsonlite::write_json(provenance, manifest_path, pretty = TRUE, auto_unbox = TRUE, null = "null"), silent = TRUE)
  quit(save = "no", status = 1L, runLast = FALSE)
})
write_tsv(source_settings, file.path(output_dir, "source_settings.tsv"))
write_tsv(plot_metadata, file.path(output_dir, "sample_metadata.tsv"))
plot_metadata_path <- file.path(output_dir, "sample_metadata.tsv")
if (length(metadata_differences)) write_tsv(do.call(rbind, metadata_differences),
                                           file.path(output_dir, "metadata_differences.tsv"))
write_tsv(data.frame(sample = metadata$sample,
  source_cohort = source_settings$source_cohort[match(vapply(metadata$sample, function(s)
    unique(all_fractions$source_run[all_fractions$sample == s]), character(1)), source_settings$source_run)],
  epidish = TRUE, uxm_full = TRUE, uxm_hematopoietic = TRUE), file.path(output_dir, "sample_coverage.tsv"))
epidish_dir <- file.path(output_dir, "epidish_r")
uxm_dir <- file.path(output_dir, "uxm")
comparison_dir <- file.path(output_dir, "comparisons")
write_tsv(all_fractions[all_fractions$method == "epidish_rpc", ], file.path(epidish_dir, "fractions_long.tsv.gz"))
write_tsv(all_warnings, file.path(epidish_dir, "epidish_warnings_all.tsv"))
write_tsv(all_native, file.path(uxm_dir, "fractions_long.tsv.gz"))
write_tsv(all_fractions, file.path(comparison_dir, "fractions_long.tsv.gz"))
for (i in seq_along(sources)) {
  destination <- file.path(comparison_dir, source_settings$source_cohort[[i]])
  dir.create(destination)
  for (profile in c("full", "hematopoietic")) {
    if (!file.symlink(file.path(sources[[i]], "comparisons", profile), file.path(destination, profile))) stopf("Link failed")
  }
  dir.create(file.path(destination, "summary"))
  for (name in c("comparison_metrics_all.tsv", "threshold_stability.tsv", "threshold_stability.pdf",
                 "analysis_report.md", "marker_counts.tsv", "epidish_warnings_all.tsv")) {
    source_path <- assert_file(file.path(comparison_summary_dir(sources[[i]]), name))
    if (!file.symlink(source_path, file.path(destination, "summary", name))) stopf("Summary link failed")
    records[[length(records) + 1L]] <- artifact_record(source_path)
  }
}
plot_epidish_standalone(output_dir, epidish_dir, coverage = coverage,
  fractions_path = file.path(epidish_dir, "fractions_long.tsv.gz"), metadata_path = plot_metadata_path,
  warning_path = file.path(epidish_dir, "epidish_warnings_all.tsv"),
  subtitle = "Saved cohort fits shown together | separate source masks and iteration limits; see source_settings.tsv")
# Native UXM has per-sample marker handling and no EpiDISH warning attribution.
empty_warnings <- all_warnings[FALSE, , drop = FALSE]
write_tsv(empty_warnings, file.path(uxm_dir, "plot_warnings.tsv"))
plot_epidish_standalone(output_dir, uxm_dir,
  fractions_path = file.path(uxm_dir, "fractions_long.tsv.gz"), metadata_path = plot_metadata_path,
  warning_path = file.path(uxm_dir, "plot_warnings.tsv"), method = "uxm_native_context",
  prefix = "uxm", title = "UXM native", subtitle = "Saved native UXM estimates | all samples; no refitting")
writeLines(c("# Saved deconvolution estimates: all samples", "",
  "epidish_r/: EpiDISH RPC only. uxm/: native UXM only. comparisons/: method comparisons, separated by source cohort.",
  sprintf("The plots show all %d samples. EpiDISH plots show saved coverage >=%d; tables retain every source threshold.", nrow(metadata), coverage),
  "The cohorts retain their original marker masks and RPC iteration limits (source_settings.tsv). This is a display bundle, not a joint fit.",
  "Annotations retain saved source metadata; metadata_differences.tsv records any differences from current metadata.",
  "Source-run paths accompany every fraction. Cohort-specific method comparisons are linked in place; no combined composition overview is generated.",
  "sample_coverage.tsv checks the source union against the supplied metadata. display_manifest.json records timestamps and input checksums."),
  file.path(output_dir, "README.md"))
provenance$inputs <- records
provenance$status <- "completed"
provenance$completed_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
jsonlite::write_json(provenance, manifest_path, pretty = TRUE, auto_unbox = TRUE, null = "null")
writeLines(provenance$completed_utc, file.path(output_dir, "DISPLAY_COMPLETED"))
finish_output_bundle(output_dir, script_path)
message("Saved-estimate display completed: ", normalizePath(output_dir))
