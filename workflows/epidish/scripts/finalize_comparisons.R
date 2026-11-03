#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: finalize_comparisons.R WORK_DIR RUN_NAME", call. = FALSE)
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
source(file.path(dirname(script_path), "epidish_common.R"))

work_dir <- assert_dir(args[[1L]], "EpiDISH comparison work directory")
run_name <- args[[2L]]
assert_safe_id(run_name, "run name")

if (!check_stage_receipt(work_dir, "preparation", file.path(work_dir, "provenance", c(
  "comparison_jobs.tsv", "pat_manifest.tsv", "sample_metadata.tsv", "input_manifest.json"
)))) stopf("Missing preparation receipt")
jobs_path <- file.path(work_dir, "provenance", "comparison_jobs.tsv")
jobs <- resolve_comparison_jobs(work_dir)
required_job_columns <- c(
  "profile", "coverage_threshold", "threshold_label", "marker_count",
  "comparison_dir", "matched_atlas", "native_source"
)
if (!identical(names(jobs), required_job_columns) || nrow(jobs) != 6L) {
  stopf("Expected six valid comparison jobs in %s", jobs_path)
}

pat_manifest <- read_tsv(file.path(work_dir, "provenance", "pat_manifest.tsv"))
samples <- pat_manifest$sample

final_comparison_files <- c(
  "uxm_weighted_matched_fractions.tsv", "comparison_metrics.tsv",
  "composition_comparison.pdf", "fraction_scatter.pdf", "manifest.json"
)
final_summary_files <- c(
  "fractions_long.tsv.gz", "comparison_metrics_all.tsv", "epidish_warnings_all.tsv",
  "threshold_stability.tsv", "threshold_stability.pdf", "marker_counts.tsv", "analysis_report.md"
)
final_epidish_files <- c(
  "epidish_full.pdf", "epidish_full.png", "epidish_hematopoietic.pdf", "epidish_hematopoietic.png",
  "epidish_annotated_fractions.tsv", "epidish_plot_provenance.json", "fractions_long.tsv.gz",
  "epidish_warnings_all.tsv", "marker_counts.tsv"
)
final_relative <- c(unlist(lapply(seq_len(nrow(jobs)), function(i) file.path(
  "comparisons", jobs$profile[[i]], jobs$threshold_label[[i]], final_comparison_files
))), file.path("comparisons", "summary", final_summary_files),
  file.path("epidish_r", final_epidish_files), "provenance/sessionInfo.txt", "provenance/run_manifest.json")
if (check_stage_receipt(work_dir, "finalization", file.path(work_dir, final_relative))) {
  sealed <- jsonlite::read_json(file.path(work_dir, "provenance", "run_manifest.json"), simplifyVector = FALSE)
  linked <- c(input_manifest_sha256 = "input_manifest.json", commands_log_sha256 = "commands.log",
              unweighted_nnls_sha256 = "unweighted_nnls.json", recovery_contract_sha256 = "recovery_contract.json")
  if (!is.null(sealed$matched_uxm_reuse_sha256)) linked <- c(linked, matched_uxm_reuse_sha256 = "matched_uxm_reuse.json")
  for (field in names(linked)) {
    if (is.null(sealed[[field]]) || sealed[[field]] != sha256_file(file.path(work_dir, "provenance", linked[[field]]))) {
      stopf("Finalized provenance checksum mismatch: %s", linked[[field]])
    }
  }
  if (!file.exists(file.path(work_dir, "EPIDISHR_COMPLETED"))) {
    writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), file.path(work_dir, "EPIDISHR_COMPLETED"))
  }
  message("Reused completed finalization with validated receipt")
  quit(status = 0L)
}
attempts <- jsonlite::read_json(file.path(work_dir, "provenance", "execution_attempts.json"))$attempts
attempt_id <- attempts[[length(attempts)]]$attempt_id
assert_safe_id(attempt_id, "execution attempt ID")
finalization_root <- file.path(work_dir, ".finalization_attempts", attempt_id)
if (dir.exists(finalization_root)) stopf("Finalization staging already exists for this attempt")
if (!dir.create(finalization_root, recursive = TRUE, showWarnings = FALSE)) stopf("Cannot reserve finalization staging")

read_fraction_table <- function(path, samples, cell_types, label) {
  x <- read_tsv(path)
  if (!identical(names(x), c("CellType", samples)) || !identical(x$CellType, cell_types)) {
    stopf("Unexpected fraction table shape or labels for %s: %s", label, path)
  }
  values <- as.matrix(x[, samples, drop = FALSE])
  storage.mode(values) <- "double"
  rownames(values) <- x$CellType
  colnames(values) <- samples
  validate_fraction_matrix(values, samples, cell_types, label, cells_by_samples = TRUE)
  list(wide = x, values = values)
}

all_fractions <- list()
all_metrics <- list()
all_manifests <- list()
all_warnings <- list()
fraction_index <- metric_index <- manifest_index <- warning_index <- 0L

for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, , drop = FALSE]
  comparison_dir <- assert_dir(job$comparison_dir, "comparison directory")
  derived_dir <- file.path(finalization_root, "comparisons", job$profile, job$threshold_label)
  dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
  profile <- job$profile
  if (!check_stage_receipt(work_dir, paste("prepared", profile, job$threshold_label, sep = "/"),
                           file.path(comparison_dir, prepared_stage_files)) ||
      !check_stage_receipt(work_dir, paste("nnls", profile, job$threshold_label, sep = "/"),
                           file.path(comparison_dir, "nnls_unweighted_fractions.tsv"))) {
    stopf("Missing prepared or NNLS stage receipt")
  }
  threshold <- as.integer(job$coverage_threshold)
  label <- job$threshold_label

  reference_df <- read_tsv(file.path(comparison_dir, "reference_matrix.tsv"))
  sample_df <- read_tsv(file.path(comparison_dir, "sample_regional_values.tsv.gz"))
  marker_mask <- read_tsv(file.path(comparison_dir, "marker_mask.tsv"))
  if (!identical(reference_df$region_id, sample_df$region_id) ||
      !identical(reference_df$region_id, marker_mask$region_id)) {
    stopf("Reference, sample, and marker rows differ for %s/%s", profile, label)
  }
  if (nrow(reference_df) != as.integer(job$marker_count)) {
    stopf("Marker count changed for %s/%s", profile, label)
  }
  cell_types <- setdiff(names(reference_df), "region_id")

  reference <- as.matrix(reference_df[, cell_types, drop = FALSE])
  sample_values <- as.matrix(sample_df[, samples, drop = FALSE])
  coverage_columns <- paste0("coverage_", samples)
  if (!all(coverage_columns %in% names(marker_mask))) {
    stopf("Marker mask lacks sample coverage columns for %s/%s", profile, label)
  }
  coverage <- as.matrix(marker_mask[, coverage_columns, drop = FALSE])
  storage.mode(reference) <- storage.mode(sample_values) <- storage.mode(coverage) <- "double"
  rownames(reference) <- rownames(sample_values) <- rownames(coverage) <- reference_df$region_id
  colnames(sample_values) <- colnames(coverage) <- samples

  epidish <- read_fraction_table(
    file.path(comparison_dir, "epidish_rpc_fractions.tsv"),
    samples, cell_types, "EpiDISH RPC"
  )
  unweighted <- read_fraction_table(
    file.path(comparison_dir, "nnls_unweighted_fractions.tsv"),
    samples, cell_types, "unweighted NNLS"
  )
  epidish_warnings <- read_tsv(file.path(comparison_dir, "epidish_warnings.tsv"))
  if (!identical(names(epidish_warnings), c(
    "sample", "warning", "profile", "coverage_threshold"
  ))) {
    stopf("Invalid EpiDISH warning table for %s/%s", profile, label)
  }
  if (nrow(epidish_warnings)) {
    if (any(!epidish_warnings$sample %in% samples) ||
        any(epidish_warnings$profile != profile) ||
        any(epidish_warnings$coverage_threshold != threshold)) {
      stopf("EpiDISH warning attribution is invalid for %s/%s", profile, label)
    }
    warning_index <- warning_index + 1L
    all_warnings[[warning_index]] <- epidish_warnings
  }
  native <- read_fraction_table(
    file.path(comparison_dir, "uxm_native_fractions.tsv"),
    samples, cell_types, "native weighted UXM"
  )

  raw_matched_path <- file.path(comparison_dir, ".uxm_weighted_matched.raw.csv")
  if (!check_stage_receipt(work_dir, paste("matched", profile, label, sep = "/"), raw_matched_path)) {
    stopf("Matched UXM output has no completion receipt")
  }
  raw_matched <- assert_file(raw_matched_path, "matched weighted UXM raw result")
  matched_csv <- utils::read.csv(raw_matched, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(names(matched_csv), c("CellType", samples)) ||
      !identical(matched_csv$CellType, cell_types)) {
    stopf("Matched UXM result labels differ for %s/%s", profile, label)
  }
  write_tsv(matched_csv, file.path(derived_dir, "uxm_weighted_matched_fractions.tsv"))
  matched <- read_fraction_table(
    file.path(derived_dir, "uxm_weighted_matched_fractions.tsv"),
    samples, cell_types, "threshold-matched weighted UXM"
  )

  matched_distance <- fraction_distance(
    epidish$wide, matched$wide, samples, "epidish_rpc_vs_uxm_weighted_matched"
  )
  epidish_unweighted_distance <- fraction_distance(
    epidish$wide, unweighted$wide, samples, "epidish_rpc_vs_nnls_unweighted"
  )
  weighted_unweighted_distance <- fraction_distance(
    matched$wide, unweighted$wide, samples, "uxm_weighted_matched_vs_nnls_unweighted"
  )
  native_distance <- fraction_distance(
    epidish$wide, native$wide, samples, "epidish_rpc_vs_uxm_native_context"
  )
  distance <- rbind(
    matched_distance, epidish_unweighted_distance,
    weighted_unweighted_distance, native_distance
  )
  distance$profile <- profile
  distance$coverage_threshold <- threshold

  epidish_reconstruction <- reconstruction_metrics(
    reference, sample_values, coverage, epidish$values, "epidish_rpc"
  )
  matched_reconstruction <- reconstruction_metrics(
    reference, sample_values, coverage, matched$values, "uxm_weighted_matched"
  )
  unweighted_reconstruction <- reconstruction_metrics(
    reference, sample_values, coverage, unweighted$values, "nnls_unweighted"
  )
  reconstruction <- rbind(
    epidish_reconstruction, matched_reconstruction, unweighted_reconstruction
  )
  reconstruction$profile <- profile
  reconstruction$coverage_threshold <- threshold

  metrics <- merge(
    distance,
    reshape(
      reconstruction,
      idvar = c("sample", "profile", "coverage_threshold"),
      timevar = "method",
      direction = "wide"
    ),
    by = c("sample", "profile", "coverage_threshold"), all.x = TRUE
  )
  write_tsv(metrics, file.path(derived_dir, "comparison_metrics.tsv"))

  composition <- rbind(
    fractions_long(epidish$wide, profile, threshold, "epidish_rpc"),
    fractions_long(unweighted$wide, profile, threshold, "nnls_unweighted"),
    fractions_long(matched$wide, profile, threshold, "uxm_weighted_matched"),
    fractions_long(native$wide, profile, threshold, "uxm_native_context")
  )
  composition$sample <- factor(composition$sample, levels = samples)
  composition$method <- factor(
    composition$method,
    levels = c(
      "epidish_rpc", "nnls_unweighted", "uxm_weighted_matched", "uxm_native_context"
    )
  )
  composition_plot <- ggplot2::ggplot(
    composition,
    ggplot2::aes(x = sample, y = fraction, fill = cell_type)
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
  ggplot2::ggsave(
    file.path(derived_dir, "composition_comparison.pdf"),
    composition_plot, width = 11, height = 10
  )

  scatter_pair <- function(a, b, x_method, y_method, comparison) {
    x <- merge(
      fractions_long(a, profile, threshold, x_method),
      fractions_long(b, profile, threshold, y_method),
      by = c("sample", "profile", "coverage_threshold", "cell_type"),
      suffixes = c("_x", "_y")
    )
    data.frame(
      x[, c("sample", "profile", "coverage_threshold", "cell_type")],
      fraction_x = x$fraction_x, fraction_y = x$fraction_y,
      comparison = comparison, stringsAsFactors = FALSE
    )
  }
  scatter <- rbind(
    scatter_pair(
      matched$wide, epidish$wide, "uxm_weighted_matched", "epidish_rpc",
      "EpiDISH RPC vs weighted UXM"
    ),
    scatter_pair(
      unweighted$wide, epidish$wide, "nnls_unweighted", "epidish_rpc",
      "EpiDISH RPC vs unweighted NNLS"
    ),
    scatter_pair(
      unweighted$wide, matched$wide, "nnls_unweighted", "uxm_weighted_matched",
      "Weighted UXM vs unweighted NNLS"
    )
  )
  scatter_plot <- ggplot2::ggplot(
    scatter,
    ggplot2::aes(x = fraction_x, y = fraction_y)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60") +
    ggplot2::geom_point(alpha = 0.75, size = 1.5) +
    ggplot2::facet_grid(comparison ~ sample) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = sprintf("Matched-method comparisons: %s, coverage >= %d", profile, threshold),
      x = "First method fraction", y = "Second method fraction"
    ) +
    ggplot2::theme_bw(base_size = 10)
  ggplot2::ggsave(
    file.path(derived_dir, "fraction_scatter.pdf"),
    scatter_plot, width = 11, height = 10
  )


  required_files <- c(
    "marker_mask.tsv", "coverage_summary.tsv", "sample_regional_values.tsv.gz",
    "reference_matrix.tsv", "epidish_rpc_fractions.tsv", "epidish_warnings.tsv",
    "nnls_unweighted_fractions.tsv", ".uxm_weighted_matched.raw.csv",
    "uxm_weighted_matched_fractions.tsv", "uxm_native_fractions.tsv",
    "comparison_metrics.tsv", "composition_comparison.pdf", "fraction_scatter.pdf"
  )
  file_records <- lapply(required_files, function(name) {
    path <- file.path(if (name %in% final_comparison_files) derived_dir else comparison_dir, name)
    list(name = name, bytes = unname(file.info(assert_file(path))$size), sha256 = sha256_file(path))
  })
  comparison_manifest <- list(
    profile = profile,
    coverage_threshold = threshold,
    marker_count = nrow(reference_df),
    sample_count = length(samples),
    samples = samples,
    methods = c(
      "epidish_rpc", "nnls_unweighted", "uxm_weighted_matched", "uxm_native_context"
    ),
    epidish_warning_count = nrow(epidish_warnings),
    native_source = job$native_source,
    native_source_sha256 = sha256_file(job$native_source),
    note = paste(
      "Tool-concordance and robustness comparison only; estimates are not",
      "validated BCP-ALL blast fractions. Native UXM uses its default marker handling."
    ),
    files = file_records
  )
  jsonlite::write_json(
    comparison_manifest, file.path(derived_dir, "manifest.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null"
  )

  fraction_index <- fraction_index + 1L
  all_fractions[[fraction_index]] <- composition
  metric_index <- metric_index + 1L
  all_metrics[[metric_index]] <- metrics
  manifest_index <- manifest_index + 1L
  all_manifests[[manifest_index]] <- comparison_manifest
}

summary_dir <- file.path(finalization_root, "comparisons", "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
all_fractions <- do.call(rbind, all_fractions)
all_fractions$sample <- as.character(all_fractions$sample)
all_metrics <- do.call(rbind, all_metrics)
write_tsv(all_fractions, file.path(summary_dir, "fractions_long.tsv.gz"))
write_tsv(all_metrics, file.path(summary_dir, "comparison_metrics_all.tsv"))
warnings_all <- if (length(all_warnings)) {
  do.call(rbind, all_warnings)
} else {
  data.frame(
    sample = character(), warning = character(), profile = character(),
    coverage_threshold = integer(), stringsAsFactors = FALSE
  )
}
write_tsv(warnings_all, file.path(summary_dir, "epidish_warnings_all.tsv"))

matched_methods <- c("epidish_rpc", "nnls_unweighted", "uxm_weighted_matched")
thresholds <- sort(unique(as.integer(jobs$coverage_threshold)))
threshold_pairs <- combn(thresholds, 2L, simplify = FALSE)
stability_rows <- list()
stability_index <- 0L
for (profile in c("full", "hematopoietic")) {
  for (method in matched_methods) {
    for (sample in samples) {
      subset_data <- all_fractions[
        all_fractions$profile == profile & all_fractions$method == method &
          all_fractions$sample == sample,
        , drop = FALSE
      ]
      for (pair in threshold_pairs) {
        a <- subset_data[subset_data$coverage_threshold == pair[[1L]], c("cell_type", "fraction")]
        b <- subset_data[subset_data$coverage_threshold == pair[[2L]], c("cell_type", "fraction")]
        joined <- merge(a, b, by = "cell_type", suffixes = c("_a", "_b"))
        if (nrow(joined) == 0L) stopf("Missing threshold stability values")
        delta <- joined$fraction_a - joined$fraction_b
        stability_index <- stability_index + 1L
        stability_rows[[stability_index]] <- data.frame(
          profile = profile,
          method = method,
          sample = sample,
          threshold_a = pair[[1L]],
          threshold_b = pair[[2L]],
          total_variation = 0.5 * sum(abs(delta)),
          max_absolute_difference = max(abs(delta)),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
stability <- do.call(rbind, stability_rows)
write_tsv(stability, file.path(summary_dir, "threshold_stability.tsv"))

stability$threshold_pair <- sprintf(">=%d vs >=%d", stability$threshold_a, stability$threshold_b)
stability_plot <- ggplot2::ggplot(
  stability,
  ggplot2::aes(x = sample, y = total_variation, fill = threshold_pair)
) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::facet_grid(profile ~ method, scales = "free_y") +
  ggplot2::labs(
    title = "Sensitivity of estimated fractions to regional coverage threshold",
    x = NULL, y = "Total-variation distance", fill = "Thresholds"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggplot2::ggsave(
  file.path(summary_dir, "threshold_stability.pdf"),
  stability_plot, width = 12, height = 8
)

marker_counts <- jobs[, c("profile", "coverage_threshold", "marker_count")]
write_tsv(marker_counts, file.path(summary_dir, "marker_counts.tsv"))

report_lines <- c(
  sprintf("# EpiDISH-UXM comparison: %s", run_name),
  "",
  sprintf("Samples: %s", paste(samples, collapse = ", ")),
  "",
  "The analysis compares EpiDISH RPC, unweighted NNLS, and coverage-weighted UXM using identical",
  sprintf("cohort-wide marker masks at minimum regional coverages of %s.",
          paste(thresholds, collapse = ", ")),
  "Native UXM results are included only as contextual results because their marker",
  "handling differs from the matched comparisons.",
  sprintf("EpiDISH emitted %d captured warning(s); see epidish_warnings_all.tsv.", nrow(warnings_all)),
  "",
  "These results assess tool concordance and threshold robustness. They are not",
  "validated measurements of cell-fraction accuracy or BCP-ALL blast fraction.",
  "",
  "## Retained markers",
  "",
  paste0(
    "- ", marker_counts$profile, " at coverage >= ", marker_counts$coverage_threshold,
    ": ", marker_counts$marker_count
  )
)
writeLines(report_lines, file.path(summary_dir, "analysis_report.md"), useBytes = TRUE)

dir.create(file.path(finalization_root, "provenance"), recursive = TRUE, showWarnings = FALSE)
capture.output(sessionInfo(), file = file.path(finalization_root, "provenance", "sessionInfo.txt"))

input_manifest_path <- file.path(work_dir, "provenance", "input_manifest.json")
input_manifest <- jsonlite::read_json(assert_file(input_manifest_path), simplifyVector = FALSE)
verify_artifact_record(input_manifest$atlas, "canonical atlas")
verify_artifact_record(input_manifest$metadata, "sample metadata")
verify_artifact_record(input_manifest$environment, "EpiDISH Conda environment specification")
if (input_manifest$r_version != as.character(getRversion()) ||
    input_manifest$epidish_version != as.character(utils::packageVersion("EpiDISH"))) {
  stopf("Recorded R or EpiDISH version changed during the run")
}
observed_uxm_commit <- system2(
  "git", c("-C", input_manifest$uxm_checkout$path, "rev-parse", "HEAD"),
  stdout = TRUE
)
if (length(observed_uxm_commit) != 1L ||
    observed_uxm_commit != input_manifest$uxm_checkout$commit) {
  stopf("Pinned UXM checkout changed during the run")
}
verify_artifact_record(
  list(
    path = file.path(input_manifest$full_uxm_run$path, "pat_manifest.tsv"),
    sha256 = input_manifest$full_uxm_run$pat_manifest_sha256
  ),
  "full UXM PAT manifest"
)
verify_artifact_record(
  list(
    path = file.path(input_manifest$full_uxm_run$path, "results", "lineage_fractions.csv"),
    sha256 = input_manifest$full_uxm_run$native_result_sha256
  ),
  "full native UXM result"
)
verify_artifact_record(
  list(
    path = file.path(input_manifest$hematopoietic_uxm_run$path, "pat_manifest.tsv"),
    sha256 = input_manifest$hematopoietic_uxm_run$pat_manifest_sha256
  ),
  "haematopoietic UXM PAT manifest"
)
verify_artifact_record(
  list(
    path = file.path(
      input_manifest$hematopoietic_uxm_run$path, "results", "lineage_fractions.csv"
    ),
    sha256 = input_manifest$hematopoietic_uxm_run$native_result_sha256
  ),
  "haematopoietic native UXM result"
)
for (record in input_manifest$regional_tables) {
  verify_artifact_record(record, sprintf("regional table for %s", record$sample))
}
if (length(input_manifest$pat_artifacts) != length(samples) * 3L) {
  stopf("Input manifest does not record every PAT, CSI, and LBETA artifact")
}
for (record in input_manifest$pat_artifacts) {
  verify_artifact_record(
    record, sprintf("%s artifact for %s", record$kind, record$sample)
  )
}

relative_jobs <- jobs
relative_jobs$comparison_dir <- file.path(
  "comparisons", relative_jobs$profile, relative_jobs$threshold_label
)
relative_jobs$matched_atlas <- file.path(
  relative_jobs$comparison_dir, "uxm_matched_atlas.l4.tsv"
)
# The jobs table was frozen with relative paths during preparation; keep it unchanged.
if (!identical(relative_jobs, read_tsv(jobs_path))) stopf("Prepared jobs changed during finalization")

source(file.path(dirname(script_path), "plot_epidish_standalone.R"))
epidish_dir <- file.path(finalization_root, "epidish_r")
plot_epidish_standalone(work_dir, epidish_dir,
  fractions_path = file.path(summary_dir, "fractions_long.tsv.gz"),
  warning_path = file.path(summary_dir, "epidish_warnings_all.tsv"))
write_tsv(all_fractions[all_fractions$method == "epidish_rpc", , drop = FALSE],
          file.path(epidish_dir, "fractions_long.tsv.gz"))
write_tsv(warnings_all, file.path(epidish_dir, "epidish_warnings_all.tsv"))
write_tsv(marker_counts, file.path(epidish_dir, "marker_counts.tsv"))
plot_provenance_path <- file.path(epidish_dir, "epidish_plot_provenance.json")
plot_provenance <- jsonlite::read_json(plot_provenance_path, simplifyVector = FALSE)
plot_provenance$source_result <- "."
plot_provenance$metadata$path <- "provenance/sample_metadata.tsv"
plot_provenance$fractions$path <- "comparisons/summary/fractions_long.tsv.gz"
plot_provenance$warnings$path <- "comparisons/summary/epidish_warnings_all.tsv"
jsonlite::write_json(plot_provenance, plot_provenance_path, auto_unbox = TRUE, pretty = TRUE)

summary_files <- list.files(summary_dir, full.names = TRUE)
run_manifest <- list(
  run_name = run_name,
  output_layout = "separate-methods-v2",
  summary_directory = "comparisons/summary",
  epidish_files = lapply(list.files(epidish_dir, full.names = TRUE), function(path) list(
    name = basename(path), bytes = unname(file.info(path)$size), sha256 = sha256_file(path)
  )),
  started_utc = jsonlite::read_json(file.path(work_dir, "provenance", "execution_attempts.json"))$started_utc,
  recovery_contract_sha256 = sha256_file(file.path(work_dir, "provenance", "recovery_contract.json")),
  completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  samples = samples,
  comparison_count = nrow(jobs),
  comparisons = all_manifests,
  summary_files = lapply(summary_files, function(path) list(
    name = basename(path), bytes = unname(file.info(path)$size), sha256 = sha256_file(path)
  )),
  input_manifest_sha256 = sha256_file(input_manifest_path),
  commands_log_sha256 = sha256_file(file.path(work_dir, "provenance", "commands.log")),
  matched_uxm_reuse_sha256 = {
    reuse_path <- file.path(work_dir, "provenance", "matched_uxm_reuse.json")
    if (file.exists(reuse_path)) sha256_file(reuse_path) else NULL
  },
  unweighted_nnls_sha256 = sha256_file(file.path(
    work_dir, "provenance", "unweighted_nnls.json"
  )),
  disclaimer = paste(
    "Tool-concordance and robustness comparison only; results do not establish",
    "cell-fraction accuracy and are not validated BCP-ALL blast fractions."
  )
)
jsonlite::write_json(
  run_manifest, file.path(finalization_root, "provenance", "run_manifest.json"),
  auto_unbox = TRUE, pretty = TRUE, null = "null"
)
# Publish only complete derived artifacts. A failed generation leaves its attempt staging
# intact; any partially published section without this receipt blocks recovery.
for (relative in final_relative) {
  source <- assert_file(file.path(finalization_root, relative), "staged finalization output")
  destination <- file.path(work_dir, relative)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  status <- system2("mv", c("-nT", "--", shQuote(source), shQuote(destination)))
  if (status != 0L || file.exists(source)) stopf("Finalization publication collision; retained staging: %s", source)
}
write_stage_receipt(work_dir, "finalization", file.path(work_dir, final_relative))
writeLines(
  format(Sys.time(), tz = "UTC", usetz = TRUE),
  file.path(work_dir, "EPIDISHR_COMPLETED"),
  useBytes = TRUE
)

message(sprintf("Finalized %d comparison folders", nrow(jobs)))
