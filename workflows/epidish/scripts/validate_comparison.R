#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L, 2L)) stop("Usage: validate_comparison.R RESULT_DIR [SAMPLESHEET_CSV]", call. = FALSE)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
source(file.path(dirname(script_path), "epidish_common.R"))

result_dir <- assert_dir(args[[1L]], "completed comparison result")
summary_dir <- comparison_summary_dir(result_dir)
# Derive cohort and marker counts from verified run records, not historical constants.
input_manifest <- jsonlite::read_json(
  file.path(result_dir, "provenance", "input_manifest.json"), simplifyVector = FALSE
)
expected_samples <- unlist(input_manifest$samples, use.names = FALSE)
if (!length(expected_samples) || anyNA(expected_samples) || any(!nzchar(expected_samples)) ||
    anyDuplicated(expected_samples)) stopf("Invalid input sample set")
if (length(args) == 2L) {
  requested <- read.csv(assert_file(args[[2L]]), stringsAsFactors = FALSE)$sample
  if (!identical(requested, expected_samples)) stopf("Samplesheet and comparison sample order differ")
}
expected_counts <- read_tsv(file.path(summary_dir, "marker_counts.tsv"))
keys <- paste(expected_counts$profile, expected_counts$coverage_threshold)
thresholds <- as.integer(unlist(input_manifest$coverage_thresholds))
if (length(thresholds) != 3L || anyNA(thresholds) ||
    any(thresholds < 1L) || anyDuplicated(thresholds)) {
  stopf("Expected three distinct positive coverage thresholds in input manifest")
}
required_keys <- as.vector(outer(c("full", "hematopoietic"), thresholds, paste))
if (nrow(expected_counts) != 6L || anyDuplicated(keys) ||
    !setequal(keys, required_keys) || anyNA(expected_counts$marker_count) ||
    any(expected_counts$marker_count <= 0) ||
    any(expected_counts$marker_count != as.integer(expected_counts$marker_count))) {
  stopf("Expected exactly six positive integral marker counts")
}
for (profile in c("full", "hematopoietic")) {
  counts <- expected_counts[expected_counts$profile == profile, ]
  if (any(diff(counts$marker_count[order(counts$coverage_threshold)]) > 0)) {
    stopf("Marker count increases with a stricter coverage threshold")
  }
}

if (!file.exists(file.path(result_dir, "EPIDISHR_COMPLETED"))) {
  stopf("Comparison result lacks its completion marker")
}
summary_required <- c(
  "fractions_long.tsv.gz", "comparison_metrics_all.tsv", "marker_counts.tsv",
  "threshold_stability.tsv", "threshold_stability.pdf",
  "analysis_report.md", "epidish_warnings_all.tsv"
)
for (name in summary_required) assert_file(file.path(summary_dir, name))
for (name in c(
  "commands.log", "comparison_jobs.tsv", "input_manifest.json", "pat_manifest.tsv",
  "run_manifest.json", "sample_metadata.tsv", "sessionInfo.txt"
)) {
  assert_file(file.path(result_dir, "provenance", name))
}
input_manifest <- jsonlite::read_json(
  file.path(result_dir, "provenance", "input_manifest.json"), simplifyVector = FALSE
)
has_unweighted_nnls <- !is.null(input_manifest$unweighted_nnls)
if (has_unweighted_nnls) {
  assert_file(file.path(result_dir, "provenance", "unweighted_nnls.json"))
}

required <- c(
  "marker_mask.tsv", "coverage_summary.tsv", "sample_regional_values.tsv.gz",
  "reference_matrix.tsv", "epidish_rpc_fractions.tsv", "epidish_warnings.tsv",
  "uxm_weighted_matched_fractions.tsv", "uxm_native_fractions.tsv",
  "comparison_metrics.tsv", "composition_comparison.pdf", "fraction_scatter.pdf",
  "manifest.json"
)
if (has_unweighted_nnls) required <- c(required, "nnls_unweighted_fractions.tsv")
for (profile in c("full", "hematopoietic")) {
  for (threshold in thresholds) {
    comparison_dir <- file.path(result_dir, "comparisons", profile, threshold_label(threshold))
    assert_dir(comparison_dir)
    for (name in required) assert_file(file.path(comparison_dir, name))
    marker_mask <- read_tsv(file.path(comparison_dir, "marker_mask.tsv"))
    expected_marker_count <- expected_counts$marker_count[
      expected_counts$profile == profile &
        expected_counts$coverage_threshold == threshold
    ]
    if (nrow(marker_mask) != expected_marker_count ||
        any(marker_mask$minimum_cohort_coverage < threshold)) {
      stopf("Comparison marker mask is invalid in %s", comparison_dir)
    }
    reference <- read_tsv(file.path(comparison_dir, "reference_matrix.tsv"))
    regional <- read_tsv(file.path(comparison_dir, "sample_regional_values.tsv.gz"))
    if (nrow(reference) != expected_marker_count || nrow(regional) != expected_marker_count ||
        !identical(reference$region_id, regional$region_id) ||
        !identical(reference$region_id, marker_mask$region_id) ||
        anyDuplicated(reference$region_id) ||
        !identical(names(regional), c("region_id", expected_samples))) {
      stopf("Reference and regional matrices are not aligned in %s", comparison_dir)
    }
    fraction_files <- c("epidish_rpc_fractions.tsv")
    if (has_unweighted_nnls) {
      fraction_files <- c(fraction_files, "nnls_unweighted_fractions.tsv")
    }
    fraction_files <- c(
      fraction_files, "uxm_weighted_matched_fractions.tsv", "uxm_native_fractions.tsv"
    )
    fraction_tables <- lapply(
      fraction_files, function(name) read_tsv(file.path(comparison_dir, name))
    )
    for (fractions in fraction_tables) {
      if (!identical(names(fractions), c("CellType", expected_samples)) ||
          !identical(fractions$CellType, fraction_tables[[1L]]$CellType)) {
        stopf("Comparison sample or cell-type order is invalid in %s", comparison_dir)
      }
      values <- as.matrix(fractions[, expected_samples, drop = FALSE])
      storage.mode(values) <- "double"
      validate_fraction_matrix(
        values, expected_samples, fractions$CellType, comparison_dir,
        cells_by_samples = TRUE
      )
    }
    metrics <- read_tsv(file.path(comparison_dir, "comparison_metrics.tsv"))
    required_metrics <- c(
      "total_variation", "mae", "rmse", "max_absolute_difference", "pearson",
      "spearman", "unweighted_reconstruction_rmse.epidish_rpc",
      "coverage_weighted_reconstruction_rmse.epidish_rpc",
      "unweighted_reconstruction_rmse.uxm_weighted_matched",
      "coverage_weighted_reconstruction_rmse.uxm_weighted_matched"
    )
    if (has_unweighted_nnls) {
      required_metrics <- c(
        required_metrics,
        "unweighted_reconstruction_rmse.nnls_unweighted",
        "coverage_weighted_reconstruction_rmse.nnls_unweighted"
      )
    }
    expected_comparisons <- c(
      "epidish_rpc_vs_uxm_weighted_matched",
      "epidish_rpc_vs_uxm_native_context"
    )
    if (has_unweighted_nnls) {
      expected_comparisons <- c(
        expected_comparisons,
        "epidish_rpc_vs_nnls_unweighted",
        "uxm_weighted_matched_vs_nnls_unweighted"
      )
    }
    if (nrow(metrics) != length(expected_samples) * length(expected_comparisons) ||
        !identical(sort(unique(metrics$sample)), sort(expected_samples)) ||
        !identical(sort(unique(metrics$comparison)), sort(expected_comparisons)) ||
        any(table(metrics$sample) != length(expected_comparisons)) ||
        !all(required_metrics %in% names(metrics)) ||
        any(!is.finite(as.matrix(metrics[, required_metrics, drop = FALSE])))) {
      stopf("Comparison metrics are incomplete in %s", comparison_dir)
    }
    comparison_manifest <- jsonlite::read_json(
      file.path(comparison_dir, "manifest.json"), simplifyVector = FALSE
    )
    for (record in comparison_manifest$files) {
      path <- file.path(comparison_dir, record$name)
      if (sha256_file(path) != record$sha256 ||
          unname(file.info(path)$size) != record$bytes) {
        stopf("Comparison manifest checksum is invalid for %s", path)
      }
    }
  }
}

verify_artifact_record(input_manifest$atlas, "canonical atlas")
verify_artifact_record(input_manifest$metadata, "sample metadata")
verify_artifact_record(input_manifest$environment, "EpiDISH Conda environment specification")
if (input_manifest$r_version != as.character(getRversion()) ||
    input_manifest$epidish_version != as.character(utils::packageVersion("EpiDISH"))) {
  stopf("Recorded R or EpiDISH version does not match the validation environment")
}
uxm_commit <- system2(
  "git", c("-C", input_manifest$uxm_checkout$path, "rev-parse", "HEAD"),
  stdout = TRUE
)
if (length(uxm_commit) != 1L || uxm_commit != input_manifest$uxm_checkout$commit ||
    uxm_commit != "8d0bb456742ece21802bc8e651e080869c67c944") {
  stopf("Pinned UXM checkout provenance is invalid")
}
for (record in input_manifest$regional_tables) {
  verify_artifact_record(record, sprintf("regional table for %s", record$sample))
}
if (length(input_manifest$pat_artifacts) != length(expected_samples) * 3L) {
  stopf("Input manifest lacks PAT/CSI/LBETA checksum records")
}
for (record in input_manifest$pat_artifacts) {
  verify_artifact_record(record, sprintf("%s artifact for %s", record$kind, record$sample))
}
source_records <- list(
  list(
    path = file.path(input_manifest$full_uxm_run$path, "pat_manifest.tsv"),
    sha256 = input_manifest$full_uxm_run$pat_manifest_sha256
  ),
  list(
    path = file.path(input_manifest$full_uxm_run$path, "results", "lineage_fractions.csv"),
    sha256 = input_manifest$full_uxm_run$native_result_sha256
  ),
  list(
    path = file.path(input_manifest$hematopoietic_uxm_run$path, "pat_manifest.tsv"),
    sha256 = input_manifest$hematopoietic_uxm_run$pat_manifest_sha256
  ),
  list(
    path = file.path(
      input_manifest$hematopoietic_uxm_run$path, "results", "lineage_fractions.csv"
    ),
    sha256 = input_manifest$hematopoietic_uxm_run$native_result_sha256
  )
)
for (i in seq_along(source_records)) {
  verify_artifact_record(source_records[[i]], sprintf("UXM source artifact %d", i))
}

run_manifest_path <- file.path(result_dir, "provenance", "run_manifest.json")
run_manifest <- jsonlite::read_json(run_manifest_path, simplifyVector = FALSE)
if (run_manifest$comparison_count != 6L ||
    !identical(unlist(run_manifest$samples), expected_samples) ||
    length(run_manifest$comparisons) != 6L ||
    run_manifest$input_manifest_sha256 != sha256_file(
      file.path(result_dir, "provenance", "input_manifest.json")
    ) ||
    run_manifest$commands_log_sha256 != sha256_file(
      file.path(result_dir, "provenance", "commands.log")
    )) {
  stopf("Top-level run manifest is invalid")
}
if (has_unweighted_nnls &&
    run_manifest$unweighted_nnls_sha256 != sha256_file(file.path(
      result_dir, "provenance", "unweighted_nnls.json"
    ))) {
  stopf("Unweighted NNLS provenance checksum is invalid")
}
if (identical(run_manifest$output_layout, "separate-methods-v2")) {
  if (file.exists(file.path(summary_dir, "composition_overview.pdf"))) {
    stopf("Separate-method layout must not contain a composition overview")
  }
  if (!length(run_manifest$epidish_files)) stopf("Missing EpiDISH-only file records")
  for (record in run_manifest$epidish_files) {
    path <- file.path(result_dir, "epidish_r", record$name)
    if (sha256_file(path) != record$sha256 || unname(file.info(path)$size) != record$bytes) {
      stopf("EpiDISH-only manifest checksum is invalid for %s", path)
    }
  }
}
for (record in run_manifest$summary_files) {
  path <- file.path(summary_dir, record$name)
  if (sha256_file(path) != record$sha256 || unname(file.info(path)$size) != record$bytes) {
    stopf("Summary manifest checksum is invalid for %s", path)
  }
}

copied_data <- list.files(
  result_dir, recursive = TRUE, full.names = TRUE,
  pattern = "[.](pat[.]gz|bam|cram)$", ignore.case = TRUE
)
if (length(copied_data)) stopf("Result unexpectedly contains copied PAT/BAM/CRAM data")

message("Comparison EpiDISH-UXM result passed acceptance validation")
