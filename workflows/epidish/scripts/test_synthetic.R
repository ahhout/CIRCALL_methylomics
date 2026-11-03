#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
source(file.path(dirname(script_path), "epidish_common.R"))

expect_error <- function(expr) {
  failed <- FALSE
  tryCatch(force(expr), error = function(e) failed <<- TRUE)
  if (!failed) stop("Expected an error but expression succeeded", call. = FALSE)
}

if (as.character(utils::packageVersion("EpiDISH")) != "2.26.0") {
  stop("Synthetic tests require EpiDISH 2.26.0", call. = FALSE)
}

counts <- matrix(
  c(8, 1, 1, 1, 1, 8, 0, 0, 0),
  nrow = 3L, byrow = TRUE,
  dimnames = list(NULL, c("U", "X", "M"))
)
measured <- regional_measure(counts, c("U", "M", "U"))
stopifnot(identical(measured$coverage, c(10, 10, 0)))
stopifnot(isTRUE(all.equal(measured$value[1:2], c(0.8, 0.8))))
stopifnot(is.na(measured$value[[3L]]))
expect_error(regional_measure(counts, c("U", "X", "M")))

reference <- matrix(
  c(
    0.1, 0.9,
    0.2, 0.8,
    NA, 0.7,
    0.4, 0.6,
    0.5, 0.5
  ),
  nrow = 5L, byrow = TRUE
)
coverage <- matrix(
  c(
    10, 10,
    5, 9,
    10, 10,
    9, 10,
    10, 10
  ),
  nrow = 5L, byrow = TRUE
)
scope <- c(TRUE, TRUE, TRUE, TRUE, FALSE)
stopifnot(identical(build_marker_mask(reference, coverage, scope, 5L), c(TRUE, TRUE, FALSE, TRUE, FALSE)))
stopifnot(identical(build_marker_mask(reference, coverage, scope, 10L), c(TRUE, FALSE, FALSE, FALSE, FALSE)))
expect_error(build_marker_mask(reference, coverage, scope, 0L))

temp_dir <- tempfile("epidishr-synthetic-")
dir.create(temp_dir)
on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

compressed_tsv <- file.path(temp_dir, "roundtrip.tsv.gz")
roundtrip <- data.frame(id = c("a", "b"), value = c(1, 2), stringsAsFactors = FALSE)
write_tsv(roundtrip, compressed_tsv)
stopifnot(isTRUE(all.equal(read_tsv(compressed_tsv), roundtrip, check.attributes = FALSE)))

cell_types <- c(
  "Blood-B", "Blood-T", "Blood-NK", "Blood-Mono+Macro",
  "Blood-Granul", "Eryth-prog", "Megakaryocytes"
)
atlas_frame <- data.frame(
  chr = rep("chr1", 8L),
  start = seq(100L, 800L, 100L),
  end = seq(150L, 850L, 100L),
  startCpG = seq_len(8L),
  endCpG = seq_len(8L) + 4L,
  target = c(cell_types, "Other"),
  name = paste0("chr1:", seq(100L, 800L, 100L), "-", seq(150L, 850L, 100L)),
  direction = rep("U", 8L),
  check.names = FALSE
)
for (i in seq_along(cell_types)) atlas_frame[[cell_types[[i]]]] <- seq(0.1, 0.8, 0.1)
atlas_frame$Other <- rev(seq(0.1, 0.8, 0.1))
atlas_path <- file.path(temp_dir, "atlas.tsv")
write_tsv(atlas_frame, atlas_path)
atlas <- read_atlas(atlas_path)
heme <- profile_definition("hematopoietic", atlas)
stopifnot(identical(heme$cells, cell_types), sum(heme$rows) == 7L)

duplicate_atlas <- rbind(atlas_frame, atlas_frame[1L, , drop = FALSE])
duplicate_path <- file.path(temp_dir, "duplicate.tsv")
write_tsv(duplicate_atlas, duplicate_path)
expect_error(read_atlas(duplicate_path))

homog_path <- file.path(temp_dir, "sample.homog.gz")
con <- gzfile(homog_path, "wt")
writeLines("# synthetic", con)
homog <- atlas_frame[, c("chr", "start", "end", "startCpG", "endCpG")]
homog$U <- 10L
homog$X <- 0L
homog$M <- 5L
utils::write.table(homog, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
close(con)
stopifnot(nrow(read_homog(homog_path, atlas)) == nrow(atlas))

bad_homog_path <- file.path(temp_dir, "bad.homog.gz")
homog$start[[1L]] <- homog$start[[1L]] + 1L
con <- gzfile(bad_homog_path, "wt")
writeLines("# synthetic", con)
utils::write.table(homog, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
close(con)
expect_error(read_homog(bad_homog_path, atlas))

native <- data.frame(
  CellType = c("A", "B"), sample_2 = c(0.4, 0.6), sample_1 = c(0.5, 0.5),
  check.names = FALSE
)
native_path <- file.path(temp_dir, "native.csv")
utils::write.csv(native, native_path, row.names = FALSE, quote = FALSE)
expect_error(read_native_fractions(native_path, c("sample_1", "sample_2"), c("A", "B")))

set.seed(20260904L)
n_features <- 500L
synthetic_ref <- matrix(
  stats::runif(n_features * 3L, min = 0.02, max = 0.98),
  nrow = n_features,
  dimnames = list(paste0("region_", seq_len(n_features)), c("A", "B", "C"))
)
true_fractions <- matrix(
  c(0.7, 0.2, 0.1, 0.1, 0.35, 0.55),
  nrow = 3L,
  dimnames = list(c("A", "B", "C"), c("sample_1", "sample_2"))
)
synthetic_samples <- synthetic_ref %*% true_fractions
fit <- run_epidish_rpc(
  synthetic_samples, synthetic_ref,
  colnames(synthetic_samples), colnames(synthetic_ref), maxit = 50L
)
estimated <- fit$fractions
if (max(abs(estimated - true_fractions)) > 0.02) {
  stop("EpiDISH did not recover noiseless synthetic mixtures within 0.02", call. = FALSE)
}
stopifnot(identical(names(fit$warnings), c("sample", "warning")))

message("All synthetic EpiDISH comparison tests passed")
