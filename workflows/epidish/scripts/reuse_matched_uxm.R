#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: reuse_matched_uxm.R WORK_DIR BASELINE_RESULT_DIR", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]]))
source(file.path(dirname(script_path), "epidish_common.R"))
work_dir <- assert_dir(args[[1L]], "prepared sensitivity work directory")
baseline_dir <- assert_dir(args[[2L]], "completed matched-UXM baseline")
assert_file(file.path(baseline_dir, "EPIDISHR_COMPLETED"), "baseline completion marker")

read_json <- function(path) jsonlite::read_json(assert_file(path), simplifyVector = FALSE)
if (!check_stage_receipt(work_dir, "preparation", file.path(work_dir, "provenance", c(
  "comparison_jobs.tsv", "pat_manifest.tsv", "sample_metadata.tsv", "input_manifest.json"
)))) stopf("Missing preparation receipt")
reuse_provenance_path <- file.path(work_dir, "provenance", "matched_uxm_reuse.json")
reuse_provenance_cached <- check_stage_receipt(work_dir, "matched_reuse_provenance", reuse_provenance_path)
baseline_manifest_path <- file.path(baseline_dir, "provenance", "run_manifest.json")
frozen_contract <- read_json(file.path(work_dir, "provenance", "recovery_contract.json"))
selection <- frozen_contract$baseline
if (is.null(selection) || !identical(selection$directory, baseline_dir) ||
    is.null(selection$run_manifest$sha256) || is.null(selection$run_manifest$bytes) ||
    selection$run_manifest$sha256 != sha256_file(baseline_manifest_path) ||
    selection$run_manifest$bytes != unname(file.info(baseline_manifest_path)$size)) {
  stopf("Baseline selection or manifest changed since the recovery contract was frozen")
}
baseline_manifest <- read_json(baseline_manifest_path)
baseline_input_path <- file.path(baseline_dir, "provenance", "input_manifest.json")
if (is.null(baseline_manifest$input_manifest_sha256) ||
    baseline_manifest$input_manifest_sha256 != sha256_file(baseline_input_path) ||
    is.null(baseline_manifest$commands_log_sha256) ||
    baseline_manifest$commands_log_sha256 != sha256_file(file.path(baseline_dir, "provenance", "commands.log")) ||
    is.null(baseline_manifest$completed_utc)) {
  stopf("Baseline run manifest or recorded input checksum is invalid")
}
baseline_input <- read_json(baseline_input_path)
current_input <- read_json(file.path(work_dir, "provenance", "input_manifest.json"))
require_equal <- function(a, b, label) {
  if (is.null(a) || is.null(b) || !identical(a, b)) stopf("Baseline compatibility evidence differs or is missing: %s", label)
}
for (field in c("samples", "profiles", "coverage_thresholds", "regional_measure")) {
  require_equal(current_input[[field]], baseline_input[[field]], field)
}
require_equal(current_input$atlas$sha256, baseline_input$atlas$sha256, "atlas identity")
require_equal(current_input$uxm_checkout$commit, baseline_input$uxm_checkout$commit, "UXM version")
require_equal(current_input$samples, baseline_manifest$samples, "manifest sample order")
# RPC iteration limit is deliberately absent: it does not enter the imported UXM fit.
artifact_identities <- function(input) {
  records <- input$pat_artifacts
  samples <- unlist(input$samples, use.names = FALSE)
  expected <- unlist(lapply(samples, function(sample) paste(sample, c("pat", "csi", "lbeta"), sep = "/")))
  observed <- vapply(records, function(x) paste(x$sample, x$kind, sep = "/"), character(1))
  if (length(records) != length(expected) || anyDuplicated(observed) || !setequal(expected, observed)) {
    stopf("Baseline or current input lacks unique PAT/CSI/LBETA identities")
  }
  lapply(expected, function(key) {
    record <- records[[match(key, observed)]]
    if (is.null(record$sha256) || !grepl("^[0-9a-f]{64}$", record$sha256) ||
        is.null(record$bytes) || !is.numeric(record$bytes) || record$bytes <= 0) {
      stopf("Missing recorded PAT/CSI/LBETA hash or size for %s", key)
    }
    list(key = key, sha256 = record$sha256, bytes = record$bytes)
  })
}
require_equal(artifact_identities(current_input), artifact_identities(baseline_input), "PAT/CSI/LBETA hashes and sizes")

verify_named_files <- function(directory, records) {
  if (!length(records)) stopf("Baseline lacks manifest file evidence: %s", directory)
  names <- vapply(records, `[[`, character(1), "name")
  if (anyDuplicated(names)) stopf("Duplicate baseline file records")
  for (record in records) {
    if (basename(record$name) != record$name || record$name %in% c(".", "..") ||
        is.null(record$bytes) || is.null(record$sha256)) stopf("Invalid baseline file record")
    path <- assert_file(file.path(directory, record$name))
    if (sha256_file(path) != record$sha256 || unname(file.info(path)$size) != record$bytes) {
      stopf("Baseline manifest checksum is invalid: %s", path)
    }
  }
  names
}
verify_named_files(comparison_summary_dir(baseline_dir), baseline_manifest$summary_files)

# Byte-for-byte equality after gzip decompression; gzip headers and relocation do not matter.
identical_table_contents <- function(left, right) {
  a <- if (grepl("[.]gz$", left)) gzfile(left, "rb") else file(left, "rb")
  b <- if (grepl("[.]gz$", right)) gzfile(right, "rb") else file(right, "rb")
  on.exit({close(a); close(b)}, add = TRUE)
  repeat {
    x <- readBin(a, "raw", n = 1024L * 1024L)
    y <- readBin(b, "raw", n = 1024L * 1024L)
    if (!identical(x, y)) return(FALSE)
    if (!length(x)) return(TRUE)
  }
}

jobs <- resolve_comparison_jobs(work_dir)
samples <- read_tsv(file.path(work_dir, "provenance", "pat_manifest.tsv"))$sample
require_equal(samples, unlist(current_input$samples, use.names = FALSE), "current sample order")
if (length(baseline_manifest$comparisons) != nrow(jobs) ||
    baseline_manifest$comparison_count != nrow(jobs)) stopf("Baseline comparison count differs")
imports <- list()
# This entire pass is read-only. One failed comparison prevents every import.
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, , drop = FALSE]
  current_dir <- job$comparison_dir
  if (!check_stage_receipt(work_dir, paste("prepared", job$profile, job$threshold_label, sep = "/"),
                           file.path(current_dir, prepared_stage_files))) stopf("Missing prepared-stage receipt")
  source_dir <- assert_dir(file.path(baseline_dir, "comparisons", job$profile, job$threshold_label))
  manifest <- read_json(file.path(source_dir, "manifest.json"))
  matching <- Filter(function(x) identical(x$profile, job$profile) &&
    identical(as.integer(x$coverage_threshold), as.integer(job$coverage_threshold)), baseline_manifest$comparisons)
  if (length(matching) != 1L || !identical(matching[[1L]], manifest)) {
    stopf("Baseline comparison manifest does not match its run manifest: %s", source_dir)
  }
  require_equal(manifest$samples, current_input$samples, "comparison samples")
  if (manifest$marker_count != job$marker_count || manifest$sample_count != length(samples)) {
    stopf("Baseline comparison marker/sample count differs")
  }
  verified <- verify_named_files(source_dir, manifest$files)
  for (name in c("marker_mask.tsv", "reference_matrix.tsv", "sample_regional_values.tsv.gz")) {
    if (!name %in% verified || !identical_table_contents(file.path(source_dir, name), file.path(current_dir, name))) {
      stopf("Baseline comparison input differs or is unverified: %s/%s", source_dir, name)
    }
  }
  source_path <- file.path(source_dir, "uxm_weighted_matched_fractions.tsv")
  if (!basename(source_path) %in% verified) stopf("Baseline fractions have no manifest checksum")
  expected_cells <- setdiff(names(read_tsv(file.path(current_dir, "reference_matrix.tsv"))), "region_id")
  fractions <- read_tsv(source_path)
  if (!identical(names(fractions), c("CellType", samples)) || !identical(fractions$CellType, expected_cells)) {
    stopf("Baseline matched-UXM ordered cells or samples differ: %s", source_path)
  }
  values <- as.matrix(fractions[, samples, drop = FALSE])
  storage.mode(values) <- "double"
  validate_fraction_matrix(values, samples, expected_cells, source_path, cells_by_samples = TRUE)
  raw_path <- file.path(current_dir, ".uxm_weighted_matched.raw.csv")
  stage <- paste("matched", job$profile, job$threshold_label, sep = "/")
  cached <- check_stage_receipt(work_dir, stage, raw_path)
  imports[[i]] <- list(profile = job$profile, coverage_threshold = as.integer(job$coverage_threshold),
                       source = source_path, source_sha256 = sha256_file(source_path),
                       raw_path = raw_path, stage = stage, cached = cached, fractions = fractions)
}

reused <- lapply(imports, function(item) {
  if (!item$cached) {
    utils::write.csv(item$fractions, item$raw_path, row.names = FALSE, quote = TRUE)
    write_stage_receipt(work_dir, item$stage, item$raw_path)
  }
  item[c("profile", "coverage_threshold", "source", "source_sha256")]
})
provenance_path <- file.path(work_dir, "provenance", "matched_uxm_reuse.json")
provenance <- list(baseline_result = baseline_dir,
                   baseline_run_manifest_sha256 = sha256_file(file.path(baseline_dir, "provenance", "run_manifest.json")),
                   reused_tables = reused)
if (reuse_provenance_cached) {
  require_equal(read_json(provenance_path), provenance, "baseline reuse provenance")
} else {
  jsonlite::write_json(provenance, provenance_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  write_stage_receipt(work_dir, "matched_reuse_provenance", provenance_path)
}
message(sprintf("Reused %d matched-UXM tables from a fully verified baseline", length(reused)))
