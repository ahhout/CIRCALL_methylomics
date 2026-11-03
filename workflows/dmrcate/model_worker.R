#!/usr/bin/env Rscript
started <- format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC")
args <- commandArgs(TRUE)
stopifnot(length(args)==4L)
config <- args[1]; out <- args[2]; key <- args[3]; logdir <- args[4]
writeLines(paste0("started_utc=",started),file.path(logdir,"model-execution.tsv"))
Sys.setenv(DMRCATE_TEST_MODE="1")
source(file.path(out,"code/analysis.R"))
main(config,out,model_only=TRUE,selected_contrast=key)
writeLines(capture.output(sessionInfo()),file.path(logdir,"model-sessionInfo.txt"))
write_json(list(started_utc=started,
  completed_utc=format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC"),
  stage="model",contrast=key),file.path(logdir,"model-completed.json"),
  auto_unbox=TRUE,pretty=TRUE)
