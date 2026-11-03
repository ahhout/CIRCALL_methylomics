#!/usr/bin/env Rscript
# Targeted runtime regression checks; no cohort inference.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L)
dest <- args[3]
stamp <- function() format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC")
writeLines(paste0("started_utc=",stamp()),file.path(dest,"execution.tsv"))
Sys.setenv(DMRCATE_TEST_MODE="1")
source(file.path(args[1],"analysis.R"))
runtime_preflight(args[2],dest)
good <- file.path(dest,"tiny.tsv.gz")
con <- gzfile(good,"wb"); writeLines(c("chr22\t101\t1","chr22\t201\t2"),con); close(con)
observed <- read_gzip_table(good,header=FALSE,col.names=c("chr","pos","index"))
stopifnot(nrow(observed)==2L,identical(observed$pos,c(101L,201L)))
expect_error <- function(expr,pattern) {
  message <- tryCatch({force(expr); NA_character_},error=function(e) conditionMessage(e))
  stopifnot(!is.na(message),grepl(pattern,message,fixed=TRUE))
}
bad <- file.path(dest,"truncated.tsv.gz")
raw <- readBin(good,"raw",n=file.info(good)$size)
writeBin(head(raw,-8L),bad)
expect_error(read_gzip_table(bad),"Decompression failed")
expect_error(read_gzip_table(file.path(dest,"missing.gz")),"Decompression failed")
# A missing container reference must be diagnosed before scanning sample data.
cfg <- fromJSON(args[2]); cfg$reference_manifest <- file.path(dest,"missing-reference.tsv")
bad_config <- file.path(dest,"missing-reference.json")
write_json(cfg,bad_config,auto_unbox=TRUE)
expect_error(runtime_preflight(bad_config,dest),"Unreadable container references")
cat(paste0("completed_utc=",stamp(),"\n"),file=file.path(dest,"execution.tsv"),append=TRUE)
writeLines("PASS: scratch TMPDIR, reference visibility, valid gzip, truncated gzip, missing gzip, missing reference",file.path(dest,"PASSED"))
cat(readLines(file.path(dest,"PASSED")),"\n")
