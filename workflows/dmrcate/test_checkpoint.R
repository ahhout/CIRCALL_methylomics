#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
stopifnot(length(args)==3L)
workflow <- args[1]; runtime <- args[2]; dest <- args[3]
Sys.setenv(DMRCATE_TEST_MODE="1")
source(file.path(workflow,"analysis.R"))
stamp <- function() format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC")
started <- stamp()
dir.create(dest,recursive=TRUE)
writeLines(started,file.path(dest,"started_utc"))
dir.create(file.path(dest,"code")); dir.create(file.path(dest,"prepared"))
file.copy(file.path(workflow,"analysis.R"),file.path(dest,"code/analysis.R"))
cfg <- fromJSON(file.path(workflow,"config/cohort.json")); cfg$runtime <- runtime
all_contrasts <- cfg$contrasts
cfg$contrasts <- matrix(c("cfDNA","BM"),nrow=1)
write_json(cfg,file.path(dest,"config.json"),auto_unbox=TRUE)
s <- as.data.frame(fread(file.path(workflow,"config/samples.tsv"))); s$analysis_report <- s$report
rownames(s) <- s$sample
fwrite(s,file.path(dest,"prepared/samples.tsv"),sep="\t")
set.seed(20261026)
n <- 300L
cv <- matrix(40L,n,12)
p <- rep(c(.15,.5,.8),length.out=12)
# A strong paired compartment effect produces actual regions for publication.
p <- ifelse(s$compartment=="cfDNA",.8,.2)
m <- matrix(rbinom(n*12,40,rep(p,each=n)),n,12)
b <- BSseq(chr=rep("chr22",n),pos=seq(20001001L,by=20L,length.out=n),M=m,Cov=cv,sampleNames=s$sample)
directory <- file.path(dest,"prepared/counts")
saveHDF5SummarizedExperiment(b,dir=directory)
saveRDS(directory,file.path(dest,"prepared/count_directory.rds"))
fwrite(data.frame(sample=s$sample,M=colSums(m),Cov=colSums(cv)),
  file.path(dest,"prepared/count_conservation.tsv"),sep="\t")
key <- "cfDNA_minus_BM"; d <- file.path(dest,"comparisons",key)
dir.create(d,recursive=TRUE)
comps <- c("cfDNA","BM"); jj <- which(s$compartment %in% comps)
jj <- jj[order(s$patient[jj],match(s$compartment[jj],rev(comps)))]
des <- design_for(s[jj,,drop=FALSE],comps)
obj <- b[,jj]; colData(obj) <- DataFrame(des$meta)
cp <- sequencing.annotate(obj,methdesign=des$methylation,all.cov=TRUE,fdr=cfg$fdr,coef=ncol(des$methylation))
expected <- call_regions(cp,cfg); stopifnot(!is.null(expected),length(expected@coord)>0L)
model_log <- file.path(dest,"model"); dir.create(model_log)
status <- system2(file.path(R.home("bin"),"Rscript"),shQuote(c(
  file.path(workflow,"model_worker.R"),file.path(dest,"config.json"),dest,key,model_log)))
stopifnot(status==0L,isTRUE(all.equal(cp,readRDS(file.path(d,"cpgannotated.rds")),tolerance=0)),
  !file.exists(file.path(d,"INFERENCE_READY.rds")))
before <- tools::md5sum(file.path(d,c("cpgannotated.rds","eligible_rows.rds")))
worker <- function(stage,path,success=TRUE) {
  status <- system2(file.path(R.home("bin"),"Rscript"),
    shQuote(c(file.path(workflow,"recover_checkpoint.R"),stage,dest,key,path)))
  stopifnot(if(success) status==0L else status!=0L)
}
recovery <- file.path(dest,"attempt"); dir.create(recovery)
worker("smooth",recovery)
actual <- readRDS(file.path(recovery,"smoothed.rds"))$called
stopifnot(isTRUE(all.equal(expected,actual,tolerance=0)))
worker("finish",recovery)
validate_outputs(file.path(dest,"config.json"),dest)
stopifnot(file.exists(file.path(dest,"validated_outputs.json")))
regions <- fread(file.path(d,"regions.tsv"))
values <- fread(file.path(d,"regional_sample_methylation.tsv"))
stopifnot(nrow(regions)==length(expected@coord),nrow(values)==8L*nrow(regions),
  file.exists(file.path(d,"INFERENCE_READY.rds")),
  identical(before,tools::md5sum(names(before))))
published <- tools::md5sum(list.files(d,full.names=TRUE))
worker("finish",recovery,FALSE)
stopifnot(identical(published,tools::md5sum(names(published))))
# Exercise empty smoothing without touching published inference.
cfg$contrasts <- all_contrasts
write_json(cfg,file.path(dest,"config.json"),auto_unbox=TRUE)
empty <- file.path(dest,"comparisons/BM_minus_PB"); dir.create(empty)
cp@ranges$is.sig <- rep(FALSE,n); cp@ranges$ind.fdr <- rep(1,n)
saveRDS(cp,file.path(empty,"cpgannotated.rds")); saveRDS(seq_len(n),file.path(empty,"eligible_rows.rds"))
empty_attempt <- file.path(dest,"empty"); dir.create(empty_attempt)
status <- system2(file.path(R.home("bin"),"Rscript"),shQuote(c(file.path(workflow,"recover_checkpoint.R"),"smooth",dest,"BM_minus_PB",empty_attempt)))
stopifnot(status==0L,is.null(readRDS(file.path(empty_attempt,"smoothed.rds"))$called))
write_json(list(started_utc=started,completed_utc=stamp(),checks=c(
  "stock DMRcate exact equality","fresh-process smoothing and publication",
  "fresh-process model exact equality","complete output contract",
  "checkpoint bytes preserved","completed contrast overwrite rejected","empty smoothing")),
  file.path(dest,"checks.json"),auto_unbox=TRUE,pretty=TRUE)
