#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
Sys.setenv(DMRCATE_TEST_MODE="1")
source(file.path(args[1],"analysis.R"))
dest <- args[2]; dir.create(dest,recursive=TRUE,showWarnings=FALSE)
started <- format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC")
writeLines(paste("started_utc",started,sep="\t"),file.path(dest,"execution.tsv"))
validated_files <- c("analysis.R","prepare.py","test_fixtures.R","test_prepare.py",
                     "model_worker.R","recover_checkpoint.R","run_dmrcate.sh")
write.table(data.frame(file=validated_files,
           md5=unname(tools::md5sum(file.path(args[1],validated_files)))),
           file.path(dest,"validated_code.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
f <- file.path(dest,"tiny.CpG_report.txt")
writeLines(c("chr22\t101\t+\t0\t2\tCG\tCGA","chr22\t102\t-\t0\t1\tCG\tCGT",
             "chr22\t201\t+\t2\t1\tCG\tCGA","chr22\t202\t-\t3\t4\tCG\tCGT",
             "chr22\t301\t+\t0\t0\tCG\tCGA","chr22\t302\t-\t0\t0\tCG\tCGT"),f)
b <- read.bismark(f,strandCollapse=TRUE,rmZeroCov=FALSE,BPPARAM=BiocParallel::SerialParam(),
                  BACKEND="HDF5Array",dir=file.path(dest,"fixture_counts"))
m <- as.vector(getCoverage(b,type="M")); cv <- as.vector(getCoverage(b,type="Cov"))
stopifnot(identical(start(b),c(101L,201L,301L)),all(m==c(0,5,0)),all(cv==c(3,10,0)),
          sum(m)==5,sum(cv)==13,(m/cv)[1]==0,is.nan((m/cv)[3]),identical(which(cv>=3),1:2))
stopifnot(isTRUE(all.equal(as.vector(getMeth(b,type="raw")),m/cv)))
# Whole-genome reports include human alternate/unplaced contigs and chrM,
# alongside the two controls. Validate these without expanding inference.
ref <- file.path(dest,"contigs.fa")
primary_dna <- rep("A",304L)
primary_dna[c(101L,201L,301L)] <- "C"
primary_dna[c(102L,202L,302L)] <- "G"
reference <- Biostrings::DNAStringSet(c(chr22=paste(primary_dna,collapse=""),
  chr22_alt="ACGACG",chrUn_random="ACGACG",chrM="ACGACG",
  phage_lambda="ACGACG",plasmid_puc19c="ACGACG"))
Biostrings::writeXStringSet(reference,ref)
namespace <- data.table(chr="chr22",pos=c(101L,201L,301L))
contigs <- GRanges(names(reference),IRanges(c(101L,rep(2L,5)),width=1))
checked <- validate_loci(contigs,namespace,"chr22",ref)
primary_only <- validate_loci(contigs[1],namespace,"chr22",ref)
stopifnot(nrow(checked)==6L,sum(checked$cpgs)==6L,
          nrow(primary_only)==1L,primary_only$validation=="CpG namespace",
          identical(checked$validation,c("CpG namespace",rep("combined FASTA",5))),
          sum(as.character(seqnames(contigs)) %in% "chr22")==1L)
reject <- function(x,reason,ns=namespace) {
  error <- tryCatch({validate_loci(x,ns,"chr22",ref); NULL},error=identity)
  stopifnot(inherits(error,"error"),grepl(reason,conditionMessage(error),fixed=TRUE))
}
reject(GRanges("unknown",IRanges(2L,width=1)),"Contigs absent")
reject(GRanges("chr22_alt",IRanges(3L,width=1)),"Non-CpG reference")
reject(GRanges("chrM",IRanges(6L,width=1)),"outside reference bounds")
reject(contigs,"Primary-chromosome CpGs absent",namespace[0])
reject(c(contigs,contigs[1]),"Duplicate collapsed")
meta <- data.frame(sample=paste0("s",1:8),patient=rep(paste0("p",1:4),each=2),compartment=rep(c("PB","BM"),4))
rownames(meta) <- meta$sample
set.seed(20260929)
n <- 300L
cov <- matrix(sample(25:50,n*8,replace=TRUE),n,8)
base <- rep(c(.1,.25,.55,.7),each=2)
probs <- matrix(rep(base,n),n,8,byrow=TRUE)+matrix(rep(c(0,.15),4*n),n,8,byrow=TRUE)
counts <- matrix(rbinom(n*8,as.vector(cov),as.vector(probs)),n,8)
obj <- BSseq(chr=rep("chr22",n),pos=seq(10001L,by=20L,length.out=n),M=counts,Cov=cov,sampleNames=meta$sample)
fit <- function(order=1:8,comps=c("BM","PB")) {
  des <- design_for(meta[order,,drop=FALSE],comps)
  z <- obj[,order]; colData(z) <- DataFrame(des$meta)
  sequencing.annotate(z,methdesign=des$methylation,all.cov=TRUE,fdr=.05,coef=ncol(des$methylation))
}
cp <- fit(); shuffled <- fit(c(8,3,5,2,7,1,4,6)); reverse <- fit(comps=c("PB","BM"))
stopifnot(isTRUE(all.equal(cp@ranges$diff,shuffled@ranges$diff,tolerance=1e-8)),
          isTRUE(all.equal(cp@ranges$ind.fdr,shuffled@ranges$ind.fdr,tolerance=1e-8)),
          isTRUE(all.equal(cp@ranges$diff,-reverse@ranges$diff,tolerance=1e-8)),
          mean(cp@ranges$diff)>0)
des <- design_for(meta,c("BM","PB"))
# Four patient baselines must be representable without a compartment coefficient.
y <- rep(c(-3,-1,1,3),each=2)
co <- qr.coef(qr(des$sample),y)
stopifnot(abs(tail(co,1))<1e-12,ncol(des$sample)==5L,ncol(des$methylation)==13L)
cp@ranges$is.sig <- FALSE
stopifnot(is.null(call_regions(cp,list(lambda=1000,C=2,min_cpgs=2))))
bad <- meta; bad$patient <- factor(rep("one",8))
stopifnot(inherits(try(design_for(bad,c("BM","PB")),silent=TRUE),"try-error"))
# Exercise the actual three-contrast analysis, annotation and empty-result tables
# on synthetic counts. The real namespace/count import is tested independently.
fixture <- file.path(dest,"synthetic_empty"); dir.create(file.path(fixture,"prepared"),recursive=TRUE)
config <- fromJSON(file.path(args[1],"config/smoke.json"))
config$samples <- file.path(args[1],"config/samples.tsv")
config$chrom_sizes <- file.path(args[1],"config/chrom_sizes.tsv")
s <- as.data.frame(fread(config$samples)); s$analysis_report <- s$report; rownames(s) <- s$sample
n <- 60L; cm <- matrix(0L,n,12); cc <- matrix(40L,n,12)
for(j in 1:12) {
  shift <- c(-2,2,2,-2)[match(s$patient[j],sort(unique(s$patient)))]
  direction <- switch(s$compartment[j],PB=0,BM=1,cfDNA=-1)
  cm[,j] <- rep(10:29,length.out=n)+direction*shift
}
cm[4,] <- 0L
for(j in 1:3) {cc[j,j] <- 2L; cm[j,j] <- 0L}
synthetic <- BSseq(chr=rep("chr22",n),pos=seq(20001001L,by=20L,length.out=n),M=cm,Cov=cc,sampleNames=s$sample)
colData(synthetic) <- DataFrame(s)
directory <- file.path(fixture,"prepared/counts")
saveHDF5SummarizedExperiment(synthetic,dir=directory)
saveRDS(directory,file.path(fixture,"prepared/count_directory.rds"))
fwrite(s,file.path(fixture,"prepared/samples.tsv"),sep="\t")
fwrite(data.frame(sample=s$sample,M=colSums(cm),Cov=colSums(cc)),file.path(fixture,"prepared/count_conservation.tsv"),sep="\t")
writeLines(started,file.path(fixture,"started_utc"))
dir.create(file.path(fixture,"provenance"))
for(p in c("upstream_control_qc_1.tsv","upstream_control_qc_2.tsv"))
  writeLines("synthetic_fixture_only",file.path(fixture,"provenance",p))
config$status <- "SYNTHETIC FIXTURE — not cohort results"
write_json(config,file.path(fixture,"config.json"),auto_unbox=TRUE)
main(file.path(fixture,"config.json"),fixture)
for(k in 1:3) {
  contrast <- paste(config$contrasts[k,],collapse="_minus_")
  d <- file.path(fixture,"comparisons",contrast)
  stopifnot(nrow(fread(file.path(d,"regions.tsv")))==0L,
            file.exists(file.path(d,"INFERENCE_READY.rds")))
}
universes <- lapply(1:3,function(k) readRDS(file.path(fixture,"comparisons",paste(config$contrasts[k,],collapse="_minus_"),"eligible_rows.rds")))
stopifnot(!identical(universes[[1]],universes[[2]]),!identical(universes[[2]],universes[[3]]),!identical(universes[[1]],universes[[3]]))
stopifnot(file.exists(file.path(fixture,"validated_outputs.json")),
          !dir.exists(file.path(fixture,"report")),
          !file.exists(file.path(fixture,"prepared/methrix_qc.rds")))
# Recovery must preserve already completed inference files byte-for-byte.
before <- tools::md5sum(list.files(file.path(fixture,"comparisons"),pattern="INFERENCE_READY.rds",recursive=TRUE,full.names=TRUE))
main(file.path(fixture,"config.json"),fixture)
after <- tools::md5sum(names(before)); stopifnot(identical(before,after))
writeLines(c("all reference contigs accepted; inference scope retained: PASS",
             "unknown, non-CpG, out-of-bounds, missing namespace and duplicate loci rejected: PASS",
             "strand collapsing and integer conservation: PASS","observed zero and missing coverage: PASS",
             "coverage threshold and coordinate boundary: PASS","sample-order invariance: PASS",
             "reverse contrast signs: PASS","patient blocking: PASS","empty-region handling: PASS",
             "invalid patient design rejection: PASS","three-contrast synthetic empty tables: PASS",
             "completed inference recovery: PASS"),file.path(dest,"checks.txt"))
writeLines(paste("completed_utc",format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC"),sep="\t"),file.path(dest,"COMPLETED"))
