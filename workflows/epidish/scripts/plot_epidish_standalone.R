#!/usr/bin/env Rscript
# Sourceable plotting function; CLI adds versioned figures without changing analysis outputs.
plot_epidish_standalone <- function(result_dir, output_dir, coverage = NULL,
    fractions_path = file.path(comparison_summary_dir(result_dir), "fractions_long.tsv.gz"),
    metadata_path = file.path(result_dir, "provenance", "sample_metadata.tsv"),
    warning_path = file.path(comparison_summary_dir(result_dir), "epidish_warnings_all.tsv"),
    method = "epidish_rpc", prefix = "epidish", title = "EpiDISH RPC",
    subtitle = basename(normalizePath(result_dir))) {
  started_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  fractions <- read_tsv(fractions_path)
  x <- fractions[fractions$method == method, , drop = FALSE]
  if (!is.null(coverage)) {
    if (length(coverage) != 1L || is.na(coverage) || !coverage %in% x$coverage_threshold) {
      stopf("Requested coverage is not present in saved estimates")
    }
    x <- x[x$coverage_threshold == coverage, , drop = FALSE]
  }
  metadata <- read_tsv(metadata_path)
  required <- c("sample", "patient", "compartment", "blast_pct")
  if (!all(required %in% names(metadata)) || anyDuplicated(metadata$sample)) {
    stopf("Sample metadata must contain unique samples, patient, compartment and blast_pct")
  }
  samples <- unique(as.character(x$sample))
  if (!length(samples) || !all(samples %in% metadata$sample)) stopf("Missing EpiDISH sample metadata")
  metadata <- metadata[match(samples, metadata$sample), required, drop = FALSE]
  patient_order <- unique(metadata$patient)
  compartment_order <- c("BM cellular DNA", "PB cellular DNA", "Plasma cfDNA")
  metadata <- metadata[order(match(metadata$patient, patient_order),
                             match(metadata$compartment, compartment_order)), , drop = FALSE]
  samples <- metadata$sample
  burden <- suppressWarnings(as.numeric(metadata$blast_pct))
  if (any(!is.na(metadata$blast_pct) & is.na(burden)) ||
      any(burden < 0 | burden > 100, na.rm = TRUE)) stopf("Invalid recorded blast percentage")
  cellular <- metadata$compartment %in% c("PB cellular DNA", "BM cellular DNA")
  annotation <- rep("Recorded blasts:\nNot recorded", nrow(metadata))
  annotation[cellular & !is.na(burden)] <- paste0("Recorded blasts:\n", burden[cellular & !is.na(burden)], "%")
  short_compartment <- c("BM cellular DNA" = "BM", "PB cellular DNA" = "PB", "Plasma cfDNA" = "cfDNA")
  labels <- setNames(paste(short_compartment[metadata$compartment], metadata$sample, annotation, sep = "\n"), samples)
  x$patient <- factor(metadata$patient[match(x$sample, metadata$sample)], levels = unique(metadata$patient))
  x$sample <- factor(x$sample, levels = samples)
  x$compartment <- metadata$compartment[match(x$sample, metadata$sample)]
  x$recorded_blast_pct <- burden[match(x$sample, metadata$sample)]
  thresholds <- sort(unique(x$coverage_threshold))
  x$coverage_panel <- factor(x$coverage_threshold, levels = thresholds,
                            labels = if (method == "uxm_native_context") "Native marker handling" else paste0("Coverage >=", thresholds))
  warnings <- read_tsv(warning_path)
  keys <- function(d) paste(d$profile, d$coverage_threshold, d$sample, sep = ":")
  x$fit_warning <- keys(x) %in% keys(warnings)
  cells <- sort(unique(x$cell_type))
  # Fixed atlas-wide ramps keep each cell's shade stable across profiles.
  compartments <- list(
    lymphoid = c("Blood-B", "Blood-T", "Blood-NK"),
    myeloid = c("Blood-Mono+Macro", "Blood-Granul", "Eryth-prog", "Megakaryocytes"),
    epithelial = c("Gastric-Ep", "Small-Int-Ep", "Colon-Ep", "Gallbladder", "Liver-Hep",
                   "Pancreas-Acinar", "Pancreas-Duct", "Pancreas-Alpha", "Pancreas-Beta",
                   "Pancreas-Delta", "Lung-Ep-Alveo", "Lung-Ep-Bron", "Kidney-Ep",
                   "Bladder-Ep", "Prostate-Ep", "Ovary-Ep", "Fallopian-Ep",
                   "Breast-Basal-Ep", "Breast-Luminal-Ep", "Head-Neck-Ep", "Thyroid-Ep"),
    stromal = c("Colon-Fibro", "Dermal-Fibro", "Heart-Fibro", "Adipocytes", "Bone-Osteob"),
    muscular = c("Heart-Cardio", "Skeletal-Musc", "Smooth-Musc"),
    endothelial = "Endothel",
    neural = c("Neuron", "Oligodend"),
    epidermal = "Epid-Kerat"
  )
  endpoints <- list(c("#2166AC", "#B3D5EE"), c("#B8E0B0", "#238443"),
                    c("#DEDEDE", "#555555"), c("#FFE082", "#E88724"),
                    c("#D3BCE8", "#7541A0"), "#EF9999",
                    c("#F5B6D2", "#D65B9A"), "#CBAF98")
  colors <- unlist(Map(function(types, ends) {
    setNames(if (length(types) == 1L) ends else
      grDevices::colorRampPalette(ends)(length(types)), types)
  }, unname(compartments), endpoints))
  if (length(setdiff(cells, names(colors)))) {
    stopf("Cell types missing from compartment palette: %s",
          paste(setdiff(cells, names(colors)), collapse = ", "))
  }
  cell_order <- intersect(names(colors), cells)
  legend_order <- c(setdiff(cell_order, compartments$epithelial),
                    intersect(compartments$epithelial, cell_order))
  x$cell_type <- factor(x$cell_type, levels = cell_order)
  x <- x[order(x$profile, x$coverage_threshold, x$patient, x$sample, x$cell_type), ]
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stopf("Could not create figure output directory: %s", output_dir)
  }
  paths <- file.path(output_dir, paste0(prefix, c("_full.pdf", "_full.png",
                                 "_hematopoietic.pdf", "_hematopoietic.png",
                                 "_annotated_fractions.tsv", "_plot_provenance.json")))
  if (any(file.exists(paths))) stopf("Standalone output already exists; choose a new output directory")
  progress <- list(started_utc = started_utc, status = "running", command = commandArgs())
  jsonlite::write_json(progress, paths[[6L]], auto_unbox = TRUE, pretty = TRUE)
  successful <- FALSE
  on.exit({
    if (!successful) {
      progress$status <- "failed"
      progress$ended_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
      jsonlite::write_json(progress, paths[[6L]], auto_unbox = TRUE, pretty = TRUE)
    }
  }, add = TRUE)
  for (profile in c("full", "hematopoietic")) {
    d <- x[x$profile == profile, , drop = FALSE]
    if (!nrow(d)) stopf("Missing EpiDISH profile: %s", profile)
    d <- d[order(d$sample, d$coverage_threshold, -d$fraction, d$cell_type), ]
    bar_key <- interaction(d$sample, d$coverage_threshold, drop = TRUE)
    d$stack_rank <- ave(seq_len(nrow(d)), bar_key, FUN = seq_along)
    flagged <- unique(d[d$fit_warning, c("sample", "patient", "coverage_panel")])
    plot <- ggplot2::ggplot(d, ggplot2::aes(sample, fraction, fill = cell_type, group = stack_rank)) +
      ggplot2::geom_col(width = 0.72, position = ggplot2::position_stack(reverse = TRUE)) +
      ggplot2::facet_grid(coverage_panel ~ patient, scales = "free_x", space = "free_x") +
      ggplot2::scale_x_discrete(labels = labels) +
      ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25), labels = function(v) paste0(v * 100, "%"),
                                 limits = c(0, 1.07), expand = c(0, 0)) +
      ggplot2::scale_fill_manual(values = colors,
                                 breaks = intersect(legend_order, as.character(d$cell_type))) +
      ggplot2::geom_text(data = flagged, ggplot2::aes(x = sample, y = 1.035, label = "*"),
                         inherit.aes = FALSE, size = 5) +
      ggplot2::labs(title = paste0(title, " | ", if (profile == "full") "Full atlas" else "Hematopoietic atlas"),
                    subtitle = paste0(subtitle,
                                      if (nrow(flagged)) "\n* Fit warning" else ""),
                    x = NULL, y = "Estimated reference fraction", fill = "Reference cell type") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(size = 8, lineheight = 1.15),
                     strip.text = ggplot2::element_text(face = "bold"),
                     legend.text = ggplot2::element_text(size = 8))
    for (extension in c("pdf", "png")) {
      ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_", profile, ".", extension)),
                      plot, width = 17, height = if (length(thresholds) == 1L) 7 else 10, dpi = 180, bg = "white")
    }
  }
  write_tsv(x[, c(intersect(c("source_run", "source_cohort", "rpc_maxit"), names(x)), "profile", "coverage_threshold", "patient", "sample", "compartment",
                  "recorded_blast_pct", "cell_type", "fraction", "fit_warning")], paths[[5L]])
  jsonlite::write_json(list(source_result = normalizePath(result_dir),
                           metadata = artifact_record(metadata_path),
                           fractions = artifact_record(fractions_path),
                           started_utc = started_utc,
                           status = "completed",
                           completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                           method = method,
                           warnings = artifact_record(warning_path),
                           coverage_thresholds = thresholds,
                           cell_colors = as.list(colors),
                           legend_order = legend_order,
                           color_compartments = compartments,
                           stack_order = "Descending fraction within each bar, largest at baseline",
                           note = "Display only; no fitting, aggregation or filtering of reference cell types."),
                      paths[[6L]], pretty = TRUE, auto_unbox = TRUE)
  successful <- TRUE
  invisible(paths)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  coverage <- NULL
  flag <- which(args == "--coverage")
  if (length(flag)) {
    if (length(flag) != 1L || flag == length(args) ||
        !grepl("^[1-9][0-9]*$", args[[flag + 1L]])) {
      stop("--coverage requires one positive integer", call. = FALSE)
    }
    coverage <- as.integer(args[[flag + 1L]])
    args <- args[-c(flag, flag + 1L)]
  }
  if (!length(args) %in% c(1L, 2L)) stop("Usage: plot_epidish_standalone.R COMPLETED_RESULT_DIR [NEW_OUTPUT_DIR] [--coverage N]", call. = FALSE)
  script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]])
  source(file.path(dirname(normalizePath(script_path)), "epidish_common.R"))
  assert_file(file.path(args[[1L]], "EPIDISHR_COMPLETED"), "completion marker")
  result_dir <- assert_dir(args[[1L]], "completed result directory")
  output_dir <- if (length(args) == 2L) args[[2L]] else file.path(
    dirname(result_dir), paste0(format(Sys.time(), "%Y-%m-%d_%H-%M", tz = "UTC"), "-epidish-figures-v1")
  )
  output_name <- basename(sub("/+$", "", output_dir))
  stamp <- substr(output_name, 1L, 16L)
  parsed_stamp <- as.POSIXct(strptime(stamp, "%Y-%m-%d_%H-%M", tz = "UTC"))
  if (!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}_([01][0-9]|2[0-3])-[0-5][0-9]-[A-Za-z0-9][A-Za-z0-9._-]*$", output_name) ||
      is.na(parsed_stamp) || format(parsed_stamp, "%Y-%m-%d_%H-%M", tz = "UTC") != stamp) {
    stopf("New figure directory must begin with a valid UTC YYYY-MM-DD_HH-MM timestamp: %s", output_name)
  }
  if (file.exists(output_dir) || (!is.na(Sys.readlink(output_dir)) && nzchar(Sys.readlink(output_dir)))) {
    stopf("Timestamped figure directory already exists: %s", output_dir)
  }
  reserve_output_bundle(output_dir, script_path, file.path(result_dir, "provenance", "run_manifest.json"))
  tryCatch({
    plot_epidish_standalone(result_dir, output_dir, coverage = coverage)
    finish_output_bundle(output_dir, script_path)
  }, error = function(error) {
    finish_output_bundle(output_dir, script_path, FALSE)
    stop(error)
  })
  message("Standalone figures: ", normalizePath(output_dir))
}
