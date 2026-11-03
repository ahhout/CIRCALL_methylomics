#!/usr/bin/env Rscript
# All biological settings come from the frozen JSON; no downloads/installations.
suppressPackageStartupMessages({
  library(bsseq); library(DMRcate); library(edgeR)
  library(HDF5Array); library(GenomicRanges); library(rtracklayer)
  library(data.table); library(jsonlite)
})
read_gzip_table <- function(path, ...) {
  # fread(cmd=...) can parse a partial file after gzip fails. Check the
  # producer explicitly before giving any decompressed bytes to the parser.
  expanded <- tempfile("dmrcate-table-", tmpdir=tempdir())
  errors <- tempfile("dmrcate-gzip-", tmpdir=tempdir())
  on.exit(unlink(c(expanded, errors)), add=TRUE)
  status <- system2("gzip", c("-dc", "--", shQuote(path)),
                    stdout=expanded, stderr=errors)
  if (status != 0L) stop("Decompression failed for ", path,
    " in ", tempdir(), " (exit ", status, "): ",
    paste(readLines(errors, warn=FALSE), collapse="; "), call.=FALSE)
  withCallingHandlers(fread(file=expanded, ...),
    warning=function(w) stop("Invalid decompressed table ", path, ": ",
                             conditionMessage(w), call.=FALSE))
}

runtime_preflight <- function(config_path, output_path) {
  started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC")
  cfg <- fromJSON(config_path)
  expected <- normalizePath(file.path(output_path,"tmp"), mustWork=TRUE)
  actual <- normalizePath(tempdir(), mustWork=TRUE)
  if (!startsWith(actual, paste0(expected,"/")))
    stop("R temporary directory must be inside ", expected, "; got ", actual)
  probe <- tempfile("write-probe-",tmpdir=actual)
  on.exit(unlink(probe),add=TRUE)
  writeLines("temporary storage check",probe)
  stopifnot(identical(readLines(probe),"temporary storage check"))
  paths <- unlist(cfg[c("namespace","chrom_sizes","reference_manifest","reference_fasta")])
  unreadable <- paths[file.access(paths,4L)!=0L]
  if(length(unreadable)) stop("Unreadable container references: ",paste(unreadable,collapse=", "))
  manifest <- fread(cfg$reference_manifest)
  controls <- manifest$value[manifest$field=="source_controls"]
  if(length(controls)!=1L || file.access(controls,4L)!=0L)
    stop("Unreadable control reference from reference manifest")
  disk <- system2("df",c("-Pk",shQuote(actual)),stdout=TRUE)
  if(!is.null(attr(disk,"status"))) stop("Cannot inspect temporary filesystem")
  write_json(list(started_utc=started,
    completed_utc=format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC"),
    tempdir=actual, references=c(paths,controls), disk=disk),
    file.path(output_path,paste0("runtime_preflight-",Sys.getpid(),".json")),
    pretty=TRUE,auto_unbox=TRUE)
}
design_for <- function(meta, comps) {
  meta$patient <- factor(meta$patient)
  meta$compartment <- factor(meta$compartment, levels=rev(comps))
  d <- model.matrix(~patient+compartment,meta)
  md <- edgeR::modelMatrixMeth(d)
  rownames(md) <- paste(rep(meta$sample,each=2),rep(c("C","T"),nrow(meta)),sep=".")
  stopifnot(qr(md)$rank==ncol(md), nrow(md)-ncol(md)==3L,
            all(md[seq(1,nrow(md),2),ncol(md)]==as.integer(meta$compartment==comps[1])),
            all(md[seq(2,nrow(md),2),ncol(md)]==0))
  list(meta=meta,sample=d,methylation=md)
}
validate_loci <- function(gr, namespace, primary_names, fasta) {
  loci <- data.table(chr=as.character(seqnames(gr)),pos=start(gr))
  if(anyDuplicated(loci)) stop("Duplicate collapsed CpG coordinates")
  index <- Biostrings::fasta.index(fasta,seqtype="DNA")
  ids <- sub("[[:space:]].*$","",index$desc)
  if(anyDuplicated(ids)) stop("Duplicate reference FASTA contig IDs")
  unknown <- setdiff(unique(loci$chr),ids)
  if(length(unknown)) stop("Contigs absent from approved FASTA: ",paste(unknown,collapse=", "))
  lengths <- setNames(index$seqlength,ids)
  if(any(loci$pos<1L | loci$pos>=lengths[loci$chr])) stop("CpG coordinates outside reference bounds")
  human <- loci[loci$chr %in% primary_names]
  if(nrow(namespace[human,on=.(chr,pos),nomatch=0L])!=nrow(human))
    stop("Primary-chromosome CpGs absent from the approved namespace")
  # Read the non-primary subset in one native call. Opening a new file for
  # each indexed record can repeatedly scan preceding sequence on this runtime.
  other <- loci[!loci$chr %in% primary_names]
  contigs <- unique(other$chr)
  if(length(contigs)) {
    dna <- Biostrings::readDNAStringSet(index[match(contigs,ids),,drop=FALSE])
    for(i in seq_along(contigs)) {
      valid <- start(Biostrings::matchPattern("CG",dna[[i]]))
      if(!all(other$pos[other$chr==contigs[i]] %in% valid)) stop("Non-CpG reference coordinates on ",contigs[i])
    }
  }
  data.frame(chr=unique(loci$chr),cpgs=as.integer(table(factor(loci$chr,levels=unique(loci$chr)))),
             validation=ifelse(unique(loci$chr) %in% primary_names,"CpG namespace","combined FASTA"))
}
# Match the native bsseq HDF5 row chunks to avoid repeatedly decompressing them.
chunks <- function(n,size=1000000L) split(seq_len(n),ceiling(seq_len(n)/size))
call_regions <- function(cp,cfg) {
  nchr <- table(seqnames(cp@ranges))
  smoothable <- as.character(seqnames(cp@ranges)) %in% names(nchr)[nchr>=2L]
  if(!any(cp@ranges$is.sig[smoothable])) return(NULL)
  tryCatch(dmrcate(cp,lambda=cfg$lambda,C=cfg$C,pcutoff="fdr",consec=FALSE,
                  betacutoff=NULL,min.cpgs=cfg$min_cpgs),error=function(e) {
    if(identical(conditionMessage(e),"K >= 2 is not TRUE")) return(NULL)
    stop(e)
  })
}
main <- function(config_path=commandArgs(trailingOnly=TRUE)[1],
                 output_path=commandArgs(trailingOnly=TRUE)[2],
                 model_only=FALSE, selected_contrast=NULL) {
cfg <- fromJSON(config_path); out <- normalizePath(output_path)
if (Sys.getenv("DMRCATE_TEST_MODE") != "1") {
  ownership <- fromJSON(file.path(out, ".run-ownership.json"))
  if (!identical(ownership$output, out) || !ownership$status %in% c("running", "failed") || file.exists(file.path(out, "COMPLETE")))
    stop("Analysis requires an owned incomplete run from run_dmrcate.sh")
}
keys <- apply(cfg$contrasts,1,paste,collapse="_minus_")
if(!is.null(selected_contrast) && !selected_contrast %in% keys) stop("Unknown contrast")
stopifnot(as.character(packageVersion("DMRcate"))=="3.6.0")
options(timeout=1, repos=c(CRAN="file:///disabled"))
Sys.setenv(http_proxy="http://127.0.0.1:9", https_proxy="http://127.0.0.1:9")
writeLines(capture.output(sessionInfo()), file.path(out,"sessionInfo.txt"))
stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC")
tab <- function(x,p) fwrite(as.data.table(x),p,sep="\t",na="NA",quote=FALSE)
# Keep disk-backed counts in the native bsseq/HDF5Array representation.
checkpoint <- function(x,p) {base::saveRDS(x,paste0(p,".tmp")); stopifnot(file.rename(paste0(p,".tmp"),p))}
samples <- as.data.frame(fread(file.path(out,"prepared/samples.tsv"),na.strings="NA"))
rownames(samples) <- samples$sample
stopifnot(nrow(samples)==12L, !anyDuplicated(samples$sample),
          all(table(samples$patient,samples$compartment)==1L))
count_pointer <- file.path(out,"prepared/count_directory.rds")
if (file.exists(count_pointer)) {
  b <- loadHDF5SummarizedExperiment(readRDS(count_pointer))
} else {
  countdir <- file.path(out,"prepared",paste0("counts-",format(Sys.time(),"%Y%m%dT%H%M%S",tz="UTC")))
  b <- read.bismark(samples$analysis_report,colData=DataFrame(samples),
                    strandCollapse=TRUE,rmZeroCov=FALSE,BACKEND="HDF5Array",
                    dir=countdir,BPPARAM=BiocParallel::SerialParam(),nThread=8L,
                    verbose=TRUE)
  stopifnot(identical(colnames(b),samples$sample))
  # CpG.bed.gz is the wgbstools one-based CpG dictionary, NOT a BED interval file.
  ns <- read_gzip_table(cfg$namespace,header=FALSE,
              col.names=c("chr","pos","index"))
  gr <- rowRanges(b)
  known_human <- fread(cfg$chrom_sizes,header=FALSE)[[1]]
  tab(validate_loci(gr,ns,known_human,cfg$reference_fasta),file.path(out,"prepared/reference_loci_validation.tsv"))
  M <- getCoverage(b,type="M"); Cov <- getCoverage(b,type="Cov")
  sums <- matrix(0,2,ncol(b))
  for (ii in chunks(nrow(b))) {
    m <- as.matrix(M[ii,,drop=FALSE]); cv <- as.matrix(Cov[ii,,drop=FALSE])
    stopifnot(all(is.finite(m)),all(is.finite(cv)),all(m>=0),all(cv>=m),
              all(m==floor(m)),all(cv==floor(cv)))
    sums <- sums+rbind(colSums(m),colSums(cv))
  }
  raw <- fread(file.path(out,"prepared/source_count_totals.tsv"))
  if (is.null(cfg$smoke_interval)) {
    for (j in seq_len(ncol(b))) {
      rr <- raw[raw$sample==samples$sample[j]]
      stopifnot(sums[1,j]==sum(rr$M),sums[2,j]==sum(rr$Cov))
    }
  } else {
    expected <- fread(file.path(out,"prepared/smoke_expected_counts.tsv"))
    keys <- paste(seqnames(gr),start(gr))
    for(j in seq_len(ncol(b))) {
      e <- expected[expected$sample==samples$sample[j]]
      idx <- match(paste(e$chr,e$pos),keys)
      stopifnot(!anyNA(idx), all(as.vector(M[idx,j])==e$M),all(as.vector(Cov[idx,j])==e$Cov),
                sums[1,j]==sum(e$M),sums[2,j]==sum(e$Cov))
    }
  }
  tab(data.frame(sample=samples$sample,M=sums[1,],Cov=sums[2,]),file.path(out,"prepared/count_conservation.tsv"))
  checkpoint(countdir,count_pointer)
  rm(ns,M,Cov); gc()
}
gr <- rowRanges(b); primary <- as.character(seqnames(gr)) %in% cfg$chromosomes
M <- getCoverage(b,type="M"); Cov <- getCoverage(b,type="Cov")
beta <- bsseq::getMeth(b,type="raw") # Zero coverage remains missing.
qcdir <- file.path(out,"summary"); dir.create(qcdir,showWarnings=FALSE)
if (!file.exists(file.path(qcdir,"QC_READY.rds"))) {
  thresholds <- c(1L,3L,5L,10L); retention <- list()
  for(ii in chunks(nrow(b))) {
    cv <- as.matrix(Cov[ii,,drop=FALSE]); chr <- as.character(seqnames(gr))[ii]
    for(ch in unique(chr)) {
      at <- which(chr==ch)
      for(t in thresholds) {
        present <- cv[at,,drop=FALSE]>=t
        z <- data.frame(chr=ch,threshold=t,type="sample",id=samples$sample,n=colSums(present))
        retention[[length(retention)+1L]] <- z
        for(k in seq_len(nrow(cfg$contrasts))) {
          comp <- cfg$contrasts[k,]; jj <- which(samples$compartment %in% comp)
          key <- paste(comp,collapse="_minus_")
          retention[[length(retention)+1L]] <- data.frame(chr=ch,threshold=t,type="contrast",id=key,n=sum(rowSums(present[,jj,drop=FALSE])==8L))
          for(p in sort(unique(samples$patient))) {
            j <- jj[samples$patient[jj]==p]
            retention[[length(retention)+1L]] <- data.frame(chr=ch,threshold=t,type="pair",id=paste(p,key,sep=":"),n=sum(rowSums(present[,j,drop=FALSE])==2L))
          }
        }
      }
    }
  }
  rr <- rbindlist(retention)[,.(n=sum(n)),by=.(chr,threshold,type,id)]
  tab(rr,file.path(qcdir,"coverage_retention_by_chromosome.tsv"))
  tab(rr[,.(n=sum(n)),by=.(threshold,type,id)],file.path(qcdir,"coverage_retention_all_contigs.tsv"))
  checkpoint(list(completed_utc=stamp()),file.path(qcdir,"QC_READY.rds"))
}

# Retain the approved, frozen GENCODE release; avoid ExperimentHub downloads
# and annotation changes through DMRcate::extractRanges().
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
empty_regions <- function() data.frame(coord=character(),no.cpgs=integer(),min_smoothed_fdr=numeric(),Stouffer=numeric(),HMFDR=numeric(),Fisher=numeric(),maxdiff=numeric(),meandiff=numeric())
for(k in seq_len(nrow(cfg$contrasts))) {
  comps <- cfg$contrasts[k,]; key <- paste(comps,collapse="_minus_")
  if(!is.null(selected_contrast) && key!=selected_contrast) next
  d <- file.path(out,"comparisons",key); dir.create(d,recursive=TRUE,showWarnings=FALSE)
  if(file.exists(file.path(d,"INFERENCE_READY.rds"))) next
  jj <- which(samples$compartment %in% comps)
  jj <- jj[order(samples$patient[jj],match(samples$compartment[jj],rev(comps)))]
  meta <- samples[jj,,drop=FALSE]
  designs <- design_for(meta,comps); meta <- designs$meta
  design <- designs$sample; md <- designs$methylation
  coef <- ncol(md)
  stopifnot(qr(md)$rank==ncol(md),nrow(md)-ncol(md)==3L,
            all(md[seq(1,16,2),coef]==as.integer(meta$compartment==comps[1])),
            all(md[seq(2,16,2),coef]==0))
  tab(data.frame(sample=rownames(design),design,check.names=FALSE),file.path(d,"sample_design.tsv"))
  tab(data.frame(channel=rownames(md),md,check.names=FALSE),file.path(d,"methylation_design.tsv"))
  write_json(list(numerator=comps[1],denominator=comps[2],coefficient=colnames(md)[coef],
                  model="~ patient + compartment",df=3L),file.path(d,"contrast.json"),auto_unbox=TRUE,pretty=TRUE)
  eligible <- logical(nrow(b))
  for(ii in chunks(nrow(b))) eligible[ii] <- primary[ii] & rowSums(as.matrix(Cov[ii,jj,drop=FALSE])>=cfg$min_coverage)==8L
  checkpoint(which(eligible),file.path(d,"eligible_rows.rds"))
  tab(data.frame(chr=cfg$chromosomes,available=as.integer(table(factor(as.character(seqnames(gr))[primary],levels=cfg$chromosomes))),
      eligible=as.integer(table(factor(as.character(seqnames(gr))[eligible],levels=cfg$chromosomes)))),file.path(d,"eligible_by_chromosome.tsv"))
  idx <- which(eligible)
  regions <- empty_regions(); reason <- "DMRcate completed"
  if(length(idx)>=2L) {
    # One genome-wide call. Do not split this fit or its BH correction by chromosome.
    obj <- b[idx,jj]; colData(obj) <- DataFrame(meta)
    cp <- sequencing.annotate(obj,methdesign=md,all.cov=TRUE,fdr=cfg$fdr,coef=coef)
    stopifnot(length(cp@ranges)==length(idx),
              all(is.finite(cp@ranges$ind.fdr)),all(cp@ranges$ind.fdr>=0 & cp@ranges$ind.fdr<=1))
    tab(as.data.frame(cp@ranges),file.path(d,"tested_cpg_statistics.tsv.gz"))
    checkpoint(cp,file.path(d,"cpgannotated.rds"))
    # Exit this process before smoothing; the launcher starts a clean R worker.
    if(model_only) next
    if(!any(cp@ranges$is.sig)) reason <- "No CpGs pass BH FDR 0.05; no regions called. Thresholds unchanged."
    else {
      called <- call_regions(cp,cfg)
      if(is.null(called)) reason <- "Fewer than two smoothed significant CpGs; no regions possible at min.cpgs=2."
      else {regions <- as.data.frame(setNames(lapply(names(regions),function(n) slot(called,n)),names(regions)))
            checkpoint(called,file.path(d,"dmrcate.rds"))}
    }
  } else {
    reason <- "Fewer than two eligible CpGs; inference not estimable. Thresholds unchanged."
    tab(data.frame(seqnames=character(),start=integer(),end=integer(),width=integer(),strand=character(),stat=numeric(),rawpval=numeric(),diff=numeric(),ind.fdr=numeric(),is.sig=logical()),file.path(d,"tested_cpg_statistics.tsv.gz"))
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
  rm(regions,reg); gc()
}
if(model_only) return(invisible(NULL))
validate_outputs(config_path,out)
}

validate_outputs <- function(config_path,out) {
cfg <- fromJSON(config_path)
stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC")
# Publish completion only after all essential tables and checkpoints exist.
required <- c("summary/coverage_retention_by_chromosome.tsv","summary/coverage_retention_all_contigs.tsv",
              "prepared/count_conservation.tsv","sessionInfo.txt")
for(k in seq_len(nrow(cfg$contrasts))) required <- c(required,file.path("comparisons",paste(cfg$contrasts[k,],collapse="_minus_"),c("INFERENCE_READY.rds","regions.tsv","regions.zero_based_half_open.bed","tested_cpg_statistics.tsv.gz","sample_design.tsv","methylation_design.tsv","contrast.json","eligible_rows.rds","eligible_by_chromosome.tsv","gene_overlaps.tsv","nearest_tss.tsv","regional_sample_methylation.tsv","patient_effects.tsv","regional_mean_effects.tsv","inference_status.txt")))
stopifnot(all(file.exists(file.path(out,required))), all(file.info(file.path(out,required))$size>0))
write_json(list(started_utc=readLines(file.path(out,"started_utc")),completed_utc=stamp(),required_files=required),file.path(out,"validated_outputs.json"),pretty=TRUE,auto_unbox=TRUE)
}
if(Sys.getenv("DMRCATE_TEST_MODE")!="1") main()
