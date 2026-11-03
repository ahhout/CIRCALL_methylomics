#!/usr/bin/env Rscript
# Recovery uses the run's frozen helper functions and unmodified DMRcate.
Sys.setenv(DMRCATE_TEST_MODE="1")
args <- commandArgs(TRUE)
stopifnot(length(args)==4L, args[1] %in% c("smooth","finish"))
stage <- args[1]; out <- normalizePath(args[2],mustWork=TRUE)
key <- args[3]; recovery <- normalizePath(args[4],mustWork=TRUE)
source(file.path(out,"code/analysis.R"))
stopifnot(as.character(packageVersion("DMRcate"))=="3.6.0")
cfg <- fromJSON(file.path(out,"config.json"))
keys <- apply(cfg$contrasts,1,paste,collapse="_minus_")
stopifnot(key %in% keys)
d <- file.path(out,"comparisons",key)
stamp <- function() format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC")
started <- stamp()
writeLines(capture.output(sessionInfo()),file.path(recovery,paste0(stage,"-sessionInfo.txt")))
write_json(list(started_utc=started,stage=stage,contrast=key,
  original_started_utc=readLines(file.path(out,"started_utc"))),
  file.path(recovery,paste0(stage,"-started.json")),auto_unbox=TRUE,pretty=TRUE)
options(timeout=1,repos=c(CRAN="file:///disabled"))
Sys.setenv(http_proxy="http://127.0.0.1:9",https_proxy="http://127.0.0.1:9")
cp_path <- file.path(d,"cpgannotated.rds")
result_path <- file.path(recovery,"smoothed.rds")
if(file.exists(file.path(d,"INFERENCE_READY.rds"))) stop("Contrast is already complete")
cp <- readRDS(cp_path)
idx <- readRDS(file.path(d,"eligible_rows.rds"))
stopifnot(is(cp,"CpGannotated"),length(cp@ranges)==length(idx),length(idx)>=2L,
  !anyDuplicated(idx),all(diff(idx)>0),
  all(as.character(seqnames(cp@ranges)) %in% cfg$chromosomes),
  all(is.finite(cp@ranges$stat)),all(is.finite(cp@ranges$rawpval)),
  all(is.finite(cp@ranges$ind.fdr)),all(cp@ranges$ind.fdr>=0 & cp@ranges$ind.fdr<=1),
  identical(cp@ranges$is.sig,cp@ranges$ind.fdr<cfg$fdr))
checkpoint_md5 <- unname(tools::md5sum(cp_path))
if(stage=="smooth") {
  if(file.exists(result_path)) stop("Refusing to overwrite smoothing output")
  reason <- "DMRcate completed"
  if(!any(cp@ranges$is.sig)) {
    called <- NULL
    reason <- "No CpGs pass BH FDR 0.05; no regions called. Thresholds unchanged."
  } else {
    called <- call_regions(cp,cfg)
    if(is.null(called)) reason <- "Fewer than two smoothed significant CpGs; no regions possible at min.cpgs=2."
  }
  saveRDS(list(called=called,reason=reason,checkpoint_md5=checkpoint_md5,
    started_utc=started,completed_utc=stamp()),paste0(result_path,".tmp"))
  stopifnot(file.rename(paste0(result_path,".tmp"),result_path))
} else {
  result <- readRDS(result_path)
  stopifnot(identical(result$checkpoint_md5,checkpoint_md5))
  # Refuse partial publication rather than overwriting historical outputs.
  publish <- c("dmrcate.rds","regions.tsv","regions.zero_based_half_open.bed",
    "regional_sample_methylation.tsv","patient_effects.tsv","regional_mean_effects.tsv",
    "gene_overlaps.tsv","nearest_tss.tsv","inference_status.txt","INFERENCE_READY.rds")
  if(any(file.exists(file.path(d,publish)))) stop("Existing result files require review before publication")
  samples <- as.data.frame(fread(file.path(out,"prepared/samples.tsv"),na.strings="NA"))
  comps <- cfg$contrasts[match(key,keys),]
  jj <- which(samples$compartment %in% comps)
  jj <- jj[order(samples$patient[jj],match(samples$compartment[jj],rev(comps)))]
  meta <- design_for(samples[jj,,drop=FALSE],comps)$meta
  b <- loadHDF5SummarizedExperiment(readRDS(file.path(out,"prepared/count_directory.rds")))
  stopifnot(identical(colnames(b),samples$sample))
  gr <- rowRanges(b)
  stopifnot(identical(paste(seqnames(gr)[idx],start(gr)[idx],sep=":"),names(cp@ranges)))
  Cov <- getCoverage(b,type="Cov"); beta <- bsseq::getMeth(b,type="raw")
  rm(cp); gc()
  tab <- function(x,p) fwrite(as.data.table(x),p,sep="\t",na="NA",quote=FALSE)
  checkpoint <- function(x,p) {
    if(file.exists(p)) stop("Refusing to overwrite ",p)
    saveRDS(x,paste0(p,".tmp")); stopifnot(file.rename(paste0(p,".tmp"),p))
  }
  regions <- data.frame(coord=character(),no.cpgs=integer(),min_smoothed_fdr=numeric(),
    Stouffer=numeric(),HMFDR=numeric(),Fisher=numeric(),maxdiff=numeric(),meandiff=numeric())
  called <- result$called; reason <- result$reason
  if(!is.null(called)) {
    regions <- as.data.frame(setNames(lapply(names(regions),function(n) slot(called,n)),names(regions)))
    checkpoint(called,file.path(d,"dmrcate.rds"))
  }
  # Annotation and publication below match the frozen analysis.R implementation.
a <- import(file.path(cfg$runtime,"annotations/gencode.v50.basic.annotation.gtf.gz"))
a <- a[as.character(seqnames(a)) %in% cfg$chromosomes]
genes <- a[a$type=="gene"]; tss <- a[a$type=="transcript"]
ranges(tss) <- IRanges(ifelse(as.character(strand(tss))=="-",end(tss),start(tss)),width=1)
rm(a); gc()
annotate_regions <- function(reg,path) {
  h <- findOverlaps(reg,genes,ignore.strand=TRUE); f <- genes[subjectHits(h)]
  tab(data.frame(region=as.character(names(reg)[queryHits(h)]),gene_id=f$gene_id,gene_name=f$gene_name,
                 chr=as.character(seqnames(f)),start=start(f),end=end(f)),file.path(path,"gene_overlaps.tsv"))
  h <- distanceToNearest(reg,tss,ignore.strand=TRUE,select="all"); f <- tss[subjectHits(h)]
  tab(data.frame(region=as.character(names(reg)[queryHits(h)]),gene_id=f$gene_id,gene_name=f$gene_name,
                 transcript_id=f$transcript_id,tss=start(f),distance=mcols(h)$distance),file.path(path,"nearest_tss.tsv"))
}
  if(nrow(regions)) {
    pieces <- tstrsplit(regions$coord,"[:-]"); reg <- GRanges(pieces[[1]],IRanges(as.integer(pieces[[2]]),as.integer(pieces[[3]])))
    o <- order(regions$min_smoothed_fdr,match(as.character(seqnames(reg)),cfg$chromosomes),start(reg),end(reg))
    regions <- regions[o,,drop=FALSE]; reg <- reg[o]
  } else reg <- GRanges()
  names(reg) <- regions$coord
  tab(regions,file.path(d,"regions.tsv"))
  tab(data.frame(chr=as.character(seqnames(reg)),start0=start(reg)-1L,end=end(reg),region=as.character(names(reg))),file.path(d,"regions.zero_based_half_open.bed"))
  values <- data.frame(region=character(),sample=character(),tested_cpgs=integer(),mean_fraction=numeric(),total_coverage=numeric(),minimum_coverage=numeric())
  effects <- data.frame(region=character(),patient=character(),difference_pp=numeric())
  for(i in seq_along(reg)) {
    hits <- idx[queryHits(findOverlaps(gr[idx],reg[i],ignore.strand=TRUE))]
    cv <- as.matrix(Cov[hits,jj,drop=FALSE]); be <- as.matrix(beta[hits,jj,drop=FALSE]); means <- colMeans(be)
    values <- rbind(values,data.frame(region=names(reg)[i],sample=meta$sample,tested_cpgs=length(hits),mean_fraction=means,total_coverage=colSums(cv),minimum_coverage=apply(cv,2,min)))
    for(p in levels(meta$patient)) {
      diff <- means[meta$patient==p & meta$compartment==comps[1]]-means[meta$patient==p & meta$compartment==comps[2]]
      effects <- rbind(effects,data.frame(region=names(reg)[i],patient=p,difference_pp=100*diff))
    }
  }
  tab(values,file.path(d,"regional_sample_methylation.tsv")); tab(effects,file.path(d,"patient_effects.tsv"))
  stopifnot(nrow(values)==8L*length(reg),nrow(effects)==4L*length(reg),
            all(is.finite(values$mean_fraction)),all(values$mean_fraction>=0 & values$mean_fraction<=1),
            all(values$minimum_coverage>=cfg$min_coverage),all(values$tested_cpgs>=cfg$min_cpgs))
  means <- if(nrow(effects)) as.data.table(effects)[,.(mean_difference_pp=mean(difference_pp)),by=region] else data.frame(region=character(),mean_difference_pp=numeric())
  tab(means,file.path(d,"regional_mean_effects.tsv")); annotate_regions(reg,d)
  writeLines(reason,file.path(d,"inference_status.txt"))
  checkpoint(list(completed_utc=stamp(),eligible=length(idx),regions=nrow(regions)),file.path(d,"INFERENCE_READY.rds"))
}
write_json(list(started_utc=started,completed_utc=stamp(),stage=stage,contrast=key,
  checkpoint_md5=checkpoint_md5),file.path(recovery,paste0(stage,"-completed.json")),
  auto_unbox=TRUE,pretty=TRUE)
