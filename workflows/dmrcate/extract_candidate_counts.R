#!/usr/bin/env Rscript
started <- format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC")
args <- commandArgs(TRUE)
stopifnot(length(args)==3L)
run <- args[1]; targets_path <- args[2]; out <- args[3]
writeLines(paste0("started_utc=",started),file.path(out,"patient-extraction-started.tsv"))
suppressPackageStartupMessages({library(bsseq); library(HDF5Array); library(jsonlite)})
targets <- read.delim(targets_path,stringsAsFactors=FALSE,check.names=FALSE)
samples <- read.delim(file.path(run,"prepared/samples.tsv"),stringsAsFactors=FALSE)
b <- loadHDF5SummarizedExperiment(readRDS(file.path(run,"prepared/count_directory.rds")))
stopifnot(identical(colnames(b),samples$sample))
rows <- integer(nrow(targets))
for(ch in unique(targets$chr)) {
  at <- which(targets$chr==ch)
  chromosome <- which(seqnames(b)==ch)
  local <- match(targets$position_1based[at],start(b)[chromosome])
  stopifnot(!anyNA(local))
  rows[at] <- chromosome[local]
}
for(key in unique(targets$source_contrast)) {
  eligible <- readRDS(file.path(run,"comparisons",key,"eligible_rows.rds"))
  stopifnot(all(rows[targets$source_contrast==key] %in% eligible))
  rm(eligible); gc()
}
m <- as.matrix(getCoverage(b,type="M")[rows,,drop=FALSE])
cv <- as.matrix(getCoverage(b,type="Cov")[rows,,drop=FALSE])
stopifnot(all(is.finite(m)),all(is.finite(cv)),all(m>=0),all(cv>=m),
  all(m==floor(m)),all(cv==floor(cv)))
result <- do.call(rbind,lapply(seq_len(nrow(samples)),function(j) {
  data.frame(targets,sample=samples$sample[j],methylated=m[,j],coverage=cv[,j],
    stringsAsFactors=FALSE)
}))
write.table(result,file.path(out,"patient_cpg_counts.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"patient-sessionInfo.txt"))
write_json(list(started_utc=started,
  completed_utc=format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC"),
  source_run=run,targets=nrow(targets),samples=nrow(samples),
  source_md5=as.list(tools::md5sum(c(targets_path,
    file.path(run,"prepared/count_directory.rds"),file.path(run,"prepared/samples.tsv"))))),
  file.path(out,"patient-extraction.json"),auto_unbox=TRUE,pretty=TRUE)
