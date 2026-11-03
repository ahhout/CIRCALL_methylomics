#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    paste(
      "Usage: prepare_epidish.R ATLAS FULL_UXM_RUN HEME_UXM_RUN",
      "METADATA WORK_DIR UXM_CHECKOUT RPC_MAXIT"
    ),
    call. = FALSE
  )
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
source(file.path(dirname(script_path), "epidish_common.R"))

atlas_path <- assert_file(args[[1L]], "canonical hg38 UXM atlas")
full_run <- assert_dir(args[[2L]], "full-atlas UXM run")
heme_run <- assert_dir(args[[3L]], "haematopoietic UXM run")
metadata_path <- assert_file(args[[4L]], "sample metadata")
work_dir <- normalizePath(args[[5L]], mustWork = TRUE)
assert_file(file.path(work_dir, "provenance", "recovery_contract.json"), "frozen recovery contract")
preparation_paths <- file.path(work_dir, "provenance", c(
  "comparison_jobs.tsv", "pat_manifest.tsv", "sample_metadata.tsv", "input_manifest.json"
))
if (check_stage_receipt(work_dir, "preparation", preparation_paths)) {
  jobs <- resolve_comparison_jobs(work_dir)
  for (i in seq_len(nrow(jobs))) {
    job <- jobs[i, , drop = FALSE]
    if (!check_stage_receipt(work_dir, paste("prepared", job$profile, job$threshold_label, sep = "/"),
                             file.path(job$comparison_dir, prepared_stage_files))) {
      stopf("Missing prepared comparison receipt")
    }
  }
  message("Reused completed preparation with validated receipts")
  quit(status = 0L)
}
uxm_checkout <- assert_dir(args[[6L]], "pinned UXM checkout")
rpc_maxit <- suppressWarnings(as.integer(args[[7L]]))
if (is.na(rpc_maxit) || rpc_maxit < 1L || as.character(rpc_maxit) != args[[7L]]) {
  stopf("RPC_MAXIT must be a positive integer")
}
assert_dir(file.path(uxm_checkout, ".git"), "pinned UXM Git metadata")
uxm_commit <- system2("git", c("-C", uxm_checkout, "rev-parse", "HEAD"), stdout = TRUE)
if (length(uxm_commit) != 1L || !grepl("^[0-9a-f]{40}$", uxm_commit)) {
  stopf("Could not resolve pinned UXM commit")
}

if (as.character(utils::packageVersion("EpiDISH")) != "2.26.0") {
  stopf(
    "Expected EpiDISH 2.26.0, found %s",
    as.character(utils::packageVersion("EpiDISH"))
  )
}

atlas <- read_atlas(atlas_path)
pat_manifest <- read_pat_manifest(full_run)
samples <- pat_manifest$sample

heme_manifest <- read_pat_manifest(heme_run)
if (!identical(heme_manifest, pat_manifest)) {
  stopf("Full and haematopoietic UXM PAT manifests are not identical")
}

metadata <- read_tsv(metadata_path)
if (!"sample" %in% names(metadata) || anyDuplicated(metadata$sample)) {
  stopf("Metadata must contain unique sample identifiers")
}
if (!all(samples %in% metadata$sample)) {
  stopf("One or more PAT samples are absent from metadata")
}
metadata <- metadata[match(samples, metadata$sample), , drop = FALSE]

homog_paths <- setNames(vapply(samples, function(sample) {
  find_homog(full_run, sample)
}, character(1L)), samples)

coverage <- matrix(
  0, nrow = nrow(atlas), ncol = length(samples),
  dimnames = list(atlas$region_id, samples)
)
regional_values <- matrix(
  NA_real_, nrow = nrow(atlas), ncol = length(samples),
  dimnames = list(atlas$region_id, samples)
)

for (sample in samples) {
  counts <- read_homog(homog_paths[[sample]], atlas)
  regional <- regional_measure(counts, atlas$direction)
  coverage[, sample] <- regional$coverage
  regional_values[, sample] <- regional$value
}

native_paths <- list(
  full = file.path(full_run, "results", "lineage_fractions.csv"),
  hematopoietic = file.path(heme_run, "results", "lineage_fractions.csv")
)

thresholds <- c(1L, 3L, 5L)
jobs <- list()
job_index <- 0L

for (profile in c("full", "hematopoietic")) {
  definition <- profile_definition(profile, atlas)
  cells <- definition$cells
  row_scope <- definition$rows

  reference_df <- atlas[, cells, drop = FALSE]
  reference <- matrix(NA_real_, nrow = nrow(atlas), ncol = length(cells))
  colnames(reference) <- cells
  rownames(reference) <- atlas$region_id
  for (cell in cells) {
    raw <- reference_df[[cell]]
    numeric <- suppressWarnings(as.numeric(raw))
    if (any(!is.na(raw) & is.na(numeric))) {
      stopf("Atlas cell column contains a non-numeric value: %s", cell)
    }
    reference[, cell] <- numeric
  }
  native <- read_native_fractions(native_paths[[profile]], samples, cells)

  for (threshold in thresholds) {
    label <- threshold_label(threshold)
    comparison_dir <- file.path(work_dir, "comparisons", profile, label)
    stage <- paste("prepared", profile, label, sep = "/")
    stage_paths <- file.path(comparison_dir, prepared_stage_files)
    prepared <- check_stage_receipt(work_dir, stage, stage_paths)
    if (!prepared) dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

    keep <- build_marker_mask(reference, coverage, row_scope, threshold)
    n_markers <- sum(keep)
    if (n_markers <= length(cells)) {
      stopf(
        "%s/%s retains %d markers for %d cell types; the system is not overdetermined",
        profile, label, n_markers, length(cells)
      )
    }

    selected_ids <- atlas$region_id[keep]
    ref_selected <- reference[keep, cells, drop = FALSE]
    sample_selected <- regional_values[keep, samples, drop = FALSE]
    coverage_selected <- coverage[keep, samples, drop = FALSE]
    if (any(!is.finite(sample_selected)) || any(sample_selected < 0) ||
        any(sample_selected > 1)) {
      stopf("Selected sample regional values are invalid for %s/%s", profile, label)
    }

    if (!prepared) {
      fit <- run_epidish_rpc(
        sample_selected, ref_selected, samples, cells, maxit = rpc_maxit
      )
      epidish_values <- fit$fractions
      epidish_wide <- fractions_wide(epidish_values, samples, cells)
      epidish_warnings <- fit$warnings
      epidish_warnings$profile <- rep(profile, nrow(epidish_warnings))
      epidish_warnings$coverage_threshold <- rep(threshold, nrow(epidish_warnings))

      marker_mask <- atlas[keep, c(
        "chr", "start", "end", "startCpG", "endCpG", "target", "name", "direction"
      ), drop = FALSE]
      marker_mask$region_id <- selected_ids
      marker_mask$minimum_cohort_coverage <- apply(coverage_selected, 1L, min)
      marker_mask <- cbind(
        marker_mask,
        setNames(as.data.frame(coverage_selected, check.names = FALSE), paste0("coverage_", samples))
      )
      write_tsv(marker_mask, file.path(comparison_dir, "marker_mask.tsv"))

      coverage_summary <- do.call(rbind, lapply(samples, function(sample) {
        x <- coverage_selected[, sample]
        data.frame(
          sample = sample,
          marker_count = length(x),
          total_fragments = sum(x),
          minimum = min(x),
          q25 = as.numeric(stats::quantile(x, 0.25, names = FALSE)),
          median = stats::median(x),
          mean = mean(x),
          q75 = as.numeric(stats::quantile(x, 0.75, names = FALSE)),
          maximum = max(x),
          stringsAsFactors = FALSE
        )
      }))
      write_tsv(coverage_summary, file.path(comparison_dir, "coverage_summary.tsv"))

      write_tsv(
        data.frame(region_id = selected_ids, sample_selected, check.names = FALSE),
        file.path(comparison_dir, "sample_regional_values.tsv.gz")
      )
      write_tsv(
        data.frame(region_id = selected_ids, ref_selected, check.names = FALSE),
        file.path(comparison_dir, "reference_matrix.tsv")
      )
      write_tsv(epidish_wide, file.path(comparison_dir, "epidish_rpc_fractions.tsv"))
      write_tsv(epidish_warnings, file.path(comparison_dir, "epidish_warnings.tsv"))
      write_tsv(native, file.path(comparison_dir, "uxm_native_fractions.tsv"))

      matched_atlas <- atlas[keep, c(
        "chr", "start", "end", "startCpG", "endCpG", "target", "name", "direction",
        cells
      ), drop = FALSE]
      matched_atlas_path <- file.path(comparison_dir, "uxm_matched_atlas.l4.tsv")
      write_tsv(matched_atlas, matched_atlas_path)
      write_stage_receipt(work_dir, stage, stage_paths)
    }
    matched_atlas_path <- file.path(comparison_dir, "uxm_matched_atlas.l4.tsv")

    job_index <- job_index + 1L
    jobs[[job_index]] <- data.frame(
      profile = profile,
      coverage_threshold = threshold,
      threshold_label = label,
      marker_count = n_markers,
      comparison_dir = file.path("comparisons", profile, label),
      matched_atlas = file.path("comparisons", profile, label, "uxm_matched_atlas.l4.tsv"),
      native_source = normalizePath(native_paths[[profile]], mustWork = TRUE),
      stringsAsFactors = FALSE
    )
  }
}

jobs <- do.call(rbind, jobs)
dir.create(file.path(work_dir, "provenance"), recursive = TRUE, showWarnings = FALSE)
write_tsv(jobs, file.path(work_dir, "provenance", "comparison_jobs.tsv"))
write_tsv(pat_manifest, file.path(work_dir, "provenance", "pat_manifest.tsv"))
write_tsv(metadata, file.path(work_dir, "provenance", "sample_metadata.tsv"))

frozen <- jsonlite::read_json(file.path(work_dir, "provenance", "recovery_contract.json"), simplifyVector = FALSE)
input_manifest <- list(
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  started_utc = jsonlite::read_json(file.path(work_dir, "provenance", "execution_attempts.json"))$started_utc,
  recovery_contract_sha256 = sha256_file(file.path(work_dir, "provenance", "recovery_contract.json")),
  atlas = frozen$atlas,
  environment = artifact_record(file.path(dirname(script_path), "..", "environment.yml")),
  r_version = as.character(getRversion()),
  uxm_checkout = list(path = uxm_checkout, commit = uxm_commit),
  full_uxm_run = list(
    path = full_run,
    pat_manifest_sha256 = frozen$full_manifest$sha256,
    native_result_sha256 = frozen$native_results[[1L]]$sha256
  ),
  hematopoietic_uxm_run = list(
    path = heme_run,
    pat_manifest_sha256 = frozen$heme_manifest$sha256,
    native_result_sha256 = frozen$native_results[[2L]]$sha256
  ),
  metadata = frozen$metadata,
  regional_tables = frozen$regional_tables,
  pat_artifacts = lapply(frozen$pat_artifacts, function(record) {
    # Retain the existing descriptive mtime field without hashing payloads again.
    record$modified_utc <- format(file.info(record$path)$mtime, tz = "UTC", usetz = TRUE)
    record
  }),
  samples = samples,
  profiles = c("full", "hematopoietic"),
  coverage_thresholds = thresholds,
  regional_measure = "U/(U+X+M) for U markers; M/(U+X+M) for M markers",
  epidish_method = "RPC",
  epidish_maxit = rpc_maxit,
  epidish_version = as.character(utils::packageVersion("EpiDISH")),
  unweighted_nnls = paste(
    "scipy.optimize.nnls on the same unweighted marker rows;",
    "non-negative coefficients normalized to sum one"
  )
)
jsonlite::write_json(
  input_manifest,
  file.path(work_dir, "provenance", "input_manifest.json"),
  auto_unbox = TRUE, pretty = TRUE, null = "null"
)

message(sprintf("Prepared %d EpiDISH comparisons for %d samples", nrow(jobs), length(samples)))

write_stage_receipt(work_dir, "preparation", preparation_paths)
