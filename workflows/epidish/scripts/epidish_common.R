stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

assert_file <- function(path, label = "file") {
  if (length(path) != 1L || is.na(path) || !nzchar(path) ||
      !file.exists(path) || file.info(path)$isdir || file.info(path)$size <= 0) {
    stopf("Missing or empty %s: %s", label, path)
  }
  normalizePath(path, mustWork = TRUE)
}

assert_dir <- function(path, label = "directory") {
  if (length(path) != 1L || is.na(path) || !nzchar(path) ||
      !dir.exists(path)) {
    stopf("Missing %s: %s", label, path)
  }
  normalizePath(path, mustWork = TRUE)
}

assert_safe_id <- function(value, label = "identifier") {
  if (length(value) != 1L || is.na(value) ||
      !grepl("^[A-Za-z0-9._-]+$", value)) {
    stopf("Unsafe %s: %s", label, value)
  }
  invisible(value)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- if (grepl("[.]gz$", path)) gzfile(path, open = "wt") else file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(
    x, file = con, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "NA"
  )
  invisible(path)
}

read_tsv <- function(path, ...) {
  path <- assert_file(path)
  if (grepl("[.]gz$", path)) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    return(utils::read.delim(
      con, header = TRUE, sep = "\t", quote = "", comment.char = "",
      na.strings = c("NA", ""), stringsAsFactors = FALSE,
      check.names = FALSE, ...
    ))
  }
  data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE,
                    na.strings = c("NA", ""), check.names = FALSE, ...)
}

sha256_file <- function(path) {
  path <- assert_file(path)
  output <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stopf("sha256sum failed for %s: %s", path, paste(output, collapse = " "))
  }
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

artifact_record <- function(path) {
  path <- assert_file(path)
  info <- file.info(path)
  list(
    path = path,
    bytes = unname(info$size),
    modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE),
    sha256 = sha256_file(path)
  )
}

verify_artifact_record <- function(record, label = "recorded artifact") {
  if (is.null(record$path) || is.null(record$sha256)) {
    stopf("Invalid checksum record for %s", label)
  }
  path <- assert_file(record$path, label)
  observed <- sha256_file(path)
  if (!identical(observed, record$sha256)) {
    stopf("Checksum changed for %s: %s", label, path)
  }
  if (!is.null(record$bytes) && unname(file.info(path)$size) != record$bytes) {
    stopf("Size changed for %s: %s", label, path)
  }
  invisible(TRUE)
}

region_key <- function(chr, start, end, start_cpg, end_cpg) {
  paste(chr, start, end, start_cpg, end_cpg, sep = ":")
}

atlas_region_id <- function(name, direction) {
  paste(name, direction, sep = "|")
}

regional_measure <- function(counts, direction) {
  if (!is.matrix(counts) || !identical(colnames(counts), c("U", "X", "M")) ||
      nrow(counts) != length(direction)) {
    stopf("Counts must be an n-region by U/X/M matrix matching direction")
  }
  if (any(!direction %in% c("U", "M"))) stopf("Direction must contain only U or M")
  total <- rowSums(counts)
  numerator <- ifelse(direction == "U", counts[, "U"], counts[, "M"])
  list(
    coverage = total,
    value = ifelse(total > 0, numerator / total, NA_real_)
  )
}

build_marker_mask <- function(reference, coverage, row_scope, threshold) {
  if (!is.matrix(reference) || !is.matrix(coverage) || nrow(reference) != nrow(coverage) ||
      length(row_scope) != nrow(reference)) {
    stopf("Reference, coverage, and profile scope have incompatible dimensions")
  }
  if (length(threshold) != 1L || is.na(threshold) || threshold < 1 ||
      threshold != floor(threshold)) {
    stopf("Coverage threshold must be a positive integer")
  }
  reference_valid <- apply(reference, 1L, function(x) {
    all(is.finite(x)) && all(x >= 0) && all(x <= 1)
  })
  coverage_valid <- apply(coverage >= threshold, 1L, all)
  as.logical(row_scope) & reference_valid & coverage_valid
}

read_atlas <- function(path) {
  atlas <- read_tsv(path)
  required <- c(
    "chr", "start", "end", "startCpG", "endCpG", "target", "name",
    "direction"
  )
  if (ncol(atlas) <= length(required) || !identical(names(atlas)[seq_along(required)], required)) {
    stopf("Atlas does not begin with the expected eight metadata columns: %s", path)
  }
  if (any(!atlas$direction %in% c("U", "M"))) {
    stopf("Atlas contains directions other than U or M")
  }
  atlas$key <- region_key(
    atlas$chr, atlas$start, atlas$end, atlas$startCpG, atlas$endCpG
  )
  atlas$region_id <- atlas_region_id(atlas$name, atlas$direction)
  if (anyDuplicated(atlas$key) || anyDuplicated(atlas$region_id)) {
    stopf("Atlas region coordinates and name/direction identifiers must be unique")
  }
  atlas
}

read_pat_manifest <- function(run_dir) {
  manifest_path <- file.path(run_dir, "pat_manifest.tsv")
  manifest <- read_tsv(manifest_path)
  expected <- c("sample", "pat", "csi", "lbeta")
  if (!identical(names(manifest), expected) || nrow(manifest) == 0L) {
    stopf("Invalid PAT manifest: %s", manifest_path)
  }
  if (anyDuplicated(manifest$sample)) {
    stopf("PAT manifest sample identifiers are not unique")
  }
  for (sample in manifest$sample) assert_safe_id(sample, "sample identifier")
  for (column in c("pat", "csi", "lbeta")) {
    for (relpath in manifest[[column]]) {
      if (grepl("^/", relpath) || grepl("(^|/)\\.\\.(/|$)", relpath)) {
        stopf("PAT manifest contains a non-relative or unsafe path: %s", relpath)
      }
      assert_file(file.path(run_dir, relpath), sprintf("PAT manifest %s artifact", column))
    }
  }
  manifest
}

find_homog <- function(run_dir, sample) {
  homog_dir <- file.path(run_dir, "intermediates", "l4")
  assert_dir(homog_dir, "UXM l4 intermediate directory")
  escaped <- gsub("([][{}()+*^$|\\?.])", "\\\\\\1", sample)
  candidates <- list.files(
    homog_dir,
    pattern = sprintf("^%s[.]mem[.][A-Za-z0-9]+[.]homog[.]gz$", escaped),
    full.names = TRUE
  )
  if (length(candidates) != 1L) {
    stopf("Expected exactly one regional UXM table for %s; found %d", sample, length(candidates))
  }
  assert_file(candidates[[1L]], sprintf("regional UXM table for %s", sample))
}

read_homog <- function(path, atlas) {
  con <- gzfile(assert_file(path), open = "rt")
  on.exit(close(con), add = TRUE)
  values <- utils::read.delim(
    con, header = FALSE, comment.char = "#", quote = "",
    stringsAsFactors = FALSE, check.names = FALSE
  )
  expected_names <- c("chr", "start", "end", "startCpG", "endCpG", "U", "X", "M")
  if (ncol(values) != length(expected_names)) {
    stopf("Regional UXM table must contain eight columns: %s", path)
  }
  names(values) <- expected_names
  values$key <- region_key(
    values$chr, values$start, values$end, values$startCpG, values$endCpG
  )
  if (anyDuplicated(values$key)) stopf("Duplicate coordinates in regional UXM table: %s", path)
  index <- match(atlas$key, values$key)
  if (anyNA(index) || nrow(values) != nrow(atlas)) {
    stopf("Regional UXM coordinates do not exactly match the canonical atlas: %s", path)
  }
  values <- values[index, , drop = FALSE]
  counts <- as.matrix(values[, c("U", "X", "M"), drop = FALSE])
  storage.mode(counts) <- "double"
  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) {
    stopf("Regional UXM counts must be finite non-negative integers: %s", path)
  }
  counts
}

read_native_fractions <- function(path, samples, cell_types) {
  path <- assert_file(path, "native UXM fraction table")
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(names(x), c("CellType", samples))) {
    stopf("Native UXM sample columns do not match the PAT manifest: %s", path)
  }
  if (!identical(x$CellType, cell_types)) {
    stopf("Native UXM cell types do not match the requested profile: %s", path)
  }
  values <- as.matrix(x[, samples, drop = FALSE])
  storage.mode(values) <- "double"
  validate_fraction_matrix(values, samples, cell_types, path, cells_by_samples = TRUE)
  x
}

validate_fraction_matrix <- function(
  values, samples, cell_types, label, cells_by_samples = TRUE, tolerance = 1e-4
) {
  expected_dim <- if (cells_by_samples) c(length(cell_types), length(samples)) else c(length(samples), length(cell_types))
  if (!identical(dim(values), expected_dim)) {
    stopf("Unexpected fraction matrix dimensions for %s", label)
  }
  if (any(!is.finite(values)) || any(values < -1e-10)) {
    stopf("Fractions must be finite and non-negative for %s", label)
  }
  sums <- if (cells_by_samples) colSums(values) else rowSums(values)
  if (any(abs(sums - 1) > tolerance)) {
    stopf("Fractions do not sum to one within %.1e for %s", tolerance, label)
  }
  invisible(TRUE)
}

orient_epidish_fractions <- function(est, samples, cell_types) {
  est <- as.matrix(est)
  storage.mode(est) <- "double"
  if (!is.null(rownames(est)) && !is.null(colnames(est)) &&
      setequal(rownames(est), samples) && setequal(colnames(est), cell_types)) {
    est <- t(est[samples, cell_types, drop = FALSE])
  } else if (!is.null(rownames(est)) && !is.null(colnames(est)) &&
             setequal(rownames(est), cell_types) && setequal(colnames(est), samples)) {
    est <- est[cell_types, samples, drop = FALSE]
  } else if (identical(dim(est), c(length(samples), length(cell_types)))) {
    est <- t(est)
    rownames(est) <- cell_types
    colnames(est) <- samples
  } else if (identical(dim(est), c(length(cell_types), length(samples)))) {
    rownames(est) <- cell_types
    colnames(est) <- samples
  } else {
    stopf("Could not orient EpiDISH estF output")
  }
  est <- est[cell_types, samples, drop = FALSE]
  validate_fraction_matrix(est, samples, cell_types, "EpiDISH RPC", cells_by_samples = TRUE)
  est
}

run_epidish_rpc <- function(sample_values, reference, samples, cell_types, maxit = 50L) {
  estimates <- matrix(
    NA_real_, nrow = length(cell_types), ncol = length(samples),
    dimnames = list(cell_types, samples)
  )
  warning_rows <- list()
  warning_index <- 0L
  for (sample in samples) {
    messages <- character()
    fit <- withCallingHandlers(
      EpiDISH::epidish(
        beta.m = sample_values[, sample, drop = FALSE],
        ref.m = reference,
        method = "RPC",
        maxit = maxit
      ),
      warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    estimates[, sample] <- orient_epidish_fractions(
      fit$estF, sample, cell_types
    )[, sample]
    if (length(messages)) {
      for (message in unique(messages)) {
        warning_index <- warning_index + 1L
        warning_rows[[warning_index]] <- data.frame(
          sample = sample, warning = message, stringsAsFactors = FALSE
        )
      }
    }
  }
  warnings <- if (length(warning_rows)) {
    do.call(rbind, warning_rows)
  } else {
    data.frame(sample = character(), warning = character(), stringsAsFactors = FALSE)
  }
  list(fractions = estimates, warnings = warnings)
}

fractions_wide <- function(values, samples, cell_types) {
  values <- values[cell_types, samples, drop = FALSE]
  data.frame(CellType = cell_types, values, check.names = FALSE)
}

fractions_long <- function(wide, profile, threshold, method) {
  samples <- setdiff(names(wide), "CellType")
  out <- do.call(rbind, lapply(samples, function(sample) {
    data.frame(
      sample = sample,
      profile = profile,
      coverage_threshold = threshold,
      method = method,
      cell_type = wide$CellType,
      fraction = as.numeric(wide[[sample]]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

safe_cor <- function(x, y, method) {
  if (stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x, y, method = method))
}

fraction_distance <- function(a, b, samples, comparison) {
  stopifnot(identical(a$CellType, b$CellType))
  do.call(rbind, lapply(samples, function(sample) {
    x <- as.numeric(a[[sample]])
    y <- as.numeric(b[[sample]])
    delta <- x - y
    data.frame(
      sample = sample,
      comparison = comparison,
      total_variation = 0.5 * sum(abs(delta)),
      mae = mean(abs(delta)),
      rmse = sqrt(mean(delta^2)),
      max_absolute_difference = max(abs(delta)),
      pearson = safe_cor(x, y, "pearson"),
      spearman = safe_cor(x, y, "spearman"),
      stringsAsFactors = FALSE
    )
  }))
}

reconstruction_metrics <- function(reference, sample_values, coverage, fractions, method) {
  predicted <- reference %*% fractions
  residual <- predicted - sample_values
  do.call(rbind, lapply(seq_len(ncol(sample_values)), function(i) {
    weights <- coverage[, i]
    data.frame(
      sample = colnames(sample_values)[[i]],
      method = method,
      unweighted_reconstruction_rmse = sqrt(mean(residual[, i]^2)),
      coverage_weighted_reconstruction_rmse = sqrt(
        sum((weights * residual[, i])^2) / sum(weights^2)
      ),
      stringsAsFactors = FALSE
    )
  }))
}

threshold_label <- function(value) sprintf("coverage-ge%d", as.integer(value))

profile_definition <- function(profile, atlas) {
  meta_columns <- c(
    "chr", "start", "end", "startCpG", "endCpG", "target", "name", "direction"
  )
  cell_columns <- setdiff(names(atlas), c(meta_columns, "key", "region_id"))
  if (profile == "full") {
    list(cells = cell_columns, rows = rep(TRUE, nrow(atlas)))
  } else if (profile == "hematopoietic") {
    cells <- c(
      "Blood-B", "Blood-T", "Blood-NK", "Blood-Mono+Macro",
      "Blood-Granul", "Eryth-prog", "Megakaryocytes"
    )
    missing <- setdiff(cells, cell_columns)
    if (length(missing)) stopf("Atlas lacks haematopoietic cell types: %s", paste(missing, collapse = ", "))
    list(cells = cells, rows = atlas$target %in% cells)
  } else {
    stopf("Unknown analysis profile: %s", profile)
  }
}

# New runs isolate cross-method summaries; retain read access to historical layouts.
comparison_summary_dir <- function(result_dir) {
  current <- file.path(result_dir, "comparisons", "summary")
  if (dir.exists(current)) current else file.path(result_dir, "summary")
}

# Receipts use run-root-relative paths so final publication never invalidates them.
stage_receipt_path <- function(root, stage) {
  file.path(root, "provenance", "receipts", paste0(gsub("/", "__", stage), ".json"))
}

stage_input_fingerprint <- function(root, stage) {
  fingerprint <- sha256_file(file.path(root, "provenance", "recovery_contract.json"))
  if (grepl("^(matched|nnls)/", stage)) {
    prepared <- sub("^(matched|nnls)/", "prepared/", stage)
    fingerprint <- paste(fingerprint, sha256_file(stage_receipt_path(root, prepared)), sep = ":")
  }
  fingerprint
}

stage_relative_paths <- function(root, paths) {
  prefix <- paste0(normalizePath(root, mustWork = TRUE), "/")
  vapply(paths, function(path) {
    absolute <- normalizePath(path, mustWork = FALSE)
    if (!startsWith(absolute, prefix)) stopf("Stage output is outside run root: %s", path)
    substring(absolute, nchar(prefix) + 1L)
  }, character(1))
}

check_stage_receipt <- function(root, stage, paths) {
  receipt <- stage_receipt_path(root, stage)
  if (!file.exists(receipt)) {
    links <- Sys.readlink(paths)
    if (any(file.exists(paths) | (!is.na(links) & nzchar(links)))) {
      stopf("Unreceipted partial stage %s; use a fresh run directory", stage)
    }
    return(FALSE)
  }
  data <- jsonlite::read_json(receipt, simplifyVector = FALSE)
  relative <- stage_relative_paths(root, paths)
  if (!identical(data$stage, stage) ||
      !identical(data$input_fingerprint, stage_input_fingerprint(root, stage)) ||
      !identical(vapply(data$outputs, `[[`, character(1), "path"), unname(relative))) {
    stopf("Receipt input or output list mismatch for %s", stage)
  }
  for (i in seq_along(paths)) {
    record <- data$outputs[[i]]
    if (sha256_file(paths[[i]]) != record$sha256 ||
        unname(file.info(paths[[i]])$size) != record$bytes) {
      stopf("Receipt output checksum mismatch: %s", paths[[i]])
    }
  }
  TRUE
}

write_stage_receipt <- function(root, stage, paths) {
  receipt <- stage_receipt_path(root, stage)
  if (file.exists(receipt)) stopf("Receipt already exists: %s", receipt)
  relative <- stage_relative_paths(root, paths)
  dir.create(dirname(receipt), recursive = TRUE, showWarnings = FALSE)
  records <- lapply(seq_along(paths), function(i) list(
    path = unname(relative[[i]]), bytes = unname(file.info(assert_file(paths[[i]]))$size),
    sha256 = sha256_file(paths[[i]])
  ))
  jsonlite::write_json(list(stage = stage, input_fingerprint = stage_input_fingerprint(root, stage),
                           completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                           outputs = records), paste0(receipt, ".tmp"),
                      auto_unbox = TRUE, pretty = TRUE, null = "null")
  if (!file.rename(paste0(receipt, ".tmp"), receipt)) stopf("Could not publish receipt: %s", receipt)
}

prepared_stage_files <- c(
  "marker_mask.tsv", "coverage_summary.tsv", "sample_regional_values.tsv.gz",
  "reference_matrix.tsv", "epidish_rpc_fractions.tsv", "epidish_warnings.tsv",
  "uxm_native_fractions.tsv", "uxm_matched_atlas.l4.tsv"
)

resolve_comparison_jobs <- function(root) {
  jobs <- read_tsv(file.path(root, "provenance", "comparison_jobs.tsv"))
  expected <- file.path("comparisons", jobs$profile, jobs$threshold_label)
  if (!identical(jobs$comparison_dir, expected) ||
      !identical(jobs$matched_atlas, file.path(expected, "uxm_matched_atlas.l4.tsv"))) {
    stopf("Comparison jobs must use consistent paths relative to the run root")
  }
  jobs$comparison_dir <- file.path(root, jobs$comparison_dir)
  jobs$matched_atlas <- file.path(root, jobs$matched_atlas)
  jobs
}

# A run bundle is reserved once; internal plotting sections are not separate runs.
reserve_output_bundle <- function(output_dir, script_path, inputs = character()) {
  helper <- normalizePath(file.path(dirname(script_path), "..", "..", "output_contract.py"), mustWork = TRUE)
  args <- c(helper, "reserve", output_dir, "--workflow", script_path)
  for (input in inputs) args <- c(args, "--input", input)
  status <- system2("python3", shQuote(args))
  if (status != 0L) stopf("Output reservation failed; existing outputs were retained: %s", output_dir)
  invisible(output_dir)
}
finish_output_bundle <- function(output_dir, script_path, success = TRUE) {
  helper <- normalizePath(file.path(dirname(script_path), "..", "..", "output_contract.py"), mustWork = TRUE)
  status <- system2("python3", shQuote(c(helper, if (success) "complete" else "fail", output_dir)))
  if (status != 0L) stopf("Could not record bundle state: %s", output_dir)
}
