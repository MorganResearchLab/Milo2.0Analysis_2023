#! /usr/bin/env Rscript

## Perform DA testing with Milo & edgeR

library(SingleCellExperiment)
library(miloR)
library(edgeR)
library(BiocParallel)
library(genio)
library(irlba)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-m", "--Milo"), type="character",
       	             help="The path to the Milo object")

parser <- add_option(parser, c("-r", "--random"), type="character",
                     help="Variables in meta-data to use as random effects")

parser <- add_option(parser, c("-x", "--fixed"), type="character",
                     help="Variables in meta-data to use as fixed effects")

parser <- add_option(parser, c("-s", "--solver"), type="character", default="Fisher",
                     help="GLMM Solver to use: Fisher, HE or HE-NNLS")

parser <- add_option(parser, c("-d", "--reddim"), type="character", default="PCA",
                     help="Reduced dimension in which graph was built")

parser <- add_option(parser, c("-l", "--sampleID"), type="character",
       	  	     help="Sample column used to identify individual observations")

parser <- add_option(parser, c("-g", "--glm"), action="store_true", default=FALSE,
                     help="Run with a GLM instead of GLMM - ignore random effects")

parser <- add_option(parser, c("-c", "--chunk"), type="numeric", default=1000,
                     help="How many SNPs per chunk")

parser <- add_option(parser, c("-i", "--indexChunk"), type="numeric", default=NULL,
                     help="Which SNP chunk to select")

parser <- add_option(parser, c("-p", "--plink"), type="character",
                     help="Plink file set prefix")

parser <- add_option(parser, c("-n", "--GRM"), type="character",
                     help="GCTA GRM file set prefix")

parser <- add_option(parser, c("-f", "--maf"), type="numeric",
                     help="Minimum minor allele frequency to consider SNP for analysis")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output file prefix for DA testing results")

opt <- parse_args(parser)

# GRM reader functions
sum_i=function(i) {
    return(sum(1:i))
}


GRMreader=function(filenm, flag=1) {
    xDt.bin <- paste(filenm, ".grm.bin", sep="")
    xDt.nfl <- paste(filenm, ".grm.N.bin", sep="")
    xDt.gid <- paste(filenm, ".grm.id", sep="")
    
    xDt.id <- read.table(xDt.gid)
    xDt.n <- dim(xDt.id)[1]
    xDt.grm <- readBin(file(xDt.bin, "rb"), n=xDt.n*(xDt.n+1)/2, what=numeric(0), size=4)
    
    sn <- sapply(1:xDt.n, sum_i)
    off <- xDt.grm[-sn]
    diag <- xDt.grm[sn]

    close(file(xDt.bin, "rb"))
    if(flag==1) {
        return(list(diag=diag, off=off, n=xDt.n))
    } else {
        xDt.mat <- matrix(data=NA, nrow=xDt.n, ncol=xDt.n)
        xDt.mat[upper.tri(xDt.mat)] <- off
        xDt.mat <- t(xDt.mat)
        xDt.mat[upper.tri(xDt.mat)] <- off
        diag(xDt.mat) <- diag
        dimnames(xDt.mat) <- list(xDt.id$V1, xDt.id$V2)
    close(file(xDt.bin, "rb"))
        
    return(list(mat=xDt.mat, n=xDt.n))
    }
}

# set BPPARAM
bpparam <- MulticoreParam(workers=4)
#pparam <- SerialParam()
register(bpparam)
set.seed(42)

message("Reading GRM: ", opt$GRM)
grm <- GRMreader(opt$GRM, flag=2)$mat
grm <- grm[rowSums(grm) != 0, colSums(grm) != 0]

# get rid of the duplicate donor with different IDs
drop.dup <- unique(rownames(grm)[rowSums(grm == grm[1680]) > 0])[1] # only remove one of the duplicates

grm <- grm[!rownames(grm) %in% c(drop.dup, drop.dup) ,
           !colnames(grm) %in% c(drop.dup, drop.dup) ]

message(paste0("Reading in Milo object: ", opt$Milo))
big.mnn <- readRDS(opt$Milo)

## make the meta data from the colData
meta.df <- colData(big.mnn)

message("Making meta-data")
meta.df <- meta.df[!duplicated(meta.df[, opt$sampleID]), ]
rownames(meta.df) <- meta.df[, opt$sampleID]

## read in the plink files and pull out the correct SNP chunks
in.plink <- read_plink(opt$plink)

# genotypes are stored in 'X'
n.snps <- seq_len(nrow(in.plink$bim))
snp.chunks <- split(n.snps, ceiling(seq_along(n.snps)/opt$chunk))

if(!is.null(opt$indexChunk)){
    snpchunk <- opt$indexChunk
} else{
    jobindex <- Sys.getenv('LSB_JOBINDEX') # the snp chunk is the LSB job index
    if(jobindex != '' && as.integer(jobindex) > 0) {
          snpchunk <- as.integer(jobindex)
    } else{
      snpchunk <- 1
    }
}

# select the Nth chunk
message("Subsetting to SNP chunk ", snpchunk)
x.chunk <- snp.chunks[[snpchunk]]
in.geno <- in.plink$X[x.chunk, ] # does this have dimension names?
in.snps <- in.plink$bim[x.chunk, "id", drop=TRUE]

# check for maf > threshold
maf.snps <- unlist(lapply(in.snps, FUN=function(SX, genos) {
                            sum(genos[SX, ], na.rm=TRUE)/sum(!is.na(genos[SX, ]))
			    }, genos=in.geno))
# flip the MAFs > 0.5
maf.snps[maf.snps > 0.5] <- 1 - maf.snps[maf.snps > 0.5]
names(maf.snps) <- in.snps
print(maf.snps)

# check if maf is a percentage or a proportion
if(opt$maf > 1){
  maf.t <- opt$maf/100
} else{
  maf.t <- opt$maf
}

keep.snps <- names(maf.snps)[maf.snps >= maf.t] # this will need to be user-set
message("Keeping ", length(keep.snps), " SNPs with MAF >= ", round(maf.t*100,1), "%")

# remove SNPs that segregate between covariates
all.fixed <- unlist(strsplit(opt$fixed, split=",", fixed=TRUE))
drop.snps <- c()

message("Checking for complete separation between SNPs and covariates")
drop.snps <- unique(unlist(lapply(seq_along(keep.snps), FUN=function(SNP, fixed, meta, geno, sampid){
	     		   i.geno <- data.frame("SNP"=geno[SNP, ])
			   i.geno[, sampid] <- colnames(geno)
			   colnames(i.geno) <- c(paste0("chr", SNP), sampid)
			   meta <- merge(meta, i.geno, by=sampid)

                           for(q in seq_along(fixed)){			   
			       q.tab <- table(meta[, paste0("chr", SNP)], meta[, fixed[q]])
			       print(q.tab)
                               if(any(colSums(q.tab) == 0)){
                                   return(SNP)
                               } else {
          		           return(NA)
			       }
		           }
               		   }, meta=meta.df, fixed=all.fixed, geno=in.geno, sampid=opt$sampleID)))
drop.snps <- drop.snps[!is.na(drop.snps)]

message("Dropping ", length(drop.snps), " perfectly segregated SNPs")
keep.snps <- setdiff(keep.snps, drop.snps)

if(isFALSE(opt$glm)){
    glmm.snp.list <- list()

    # write this for bplapply?
    glmm.snp.list <- lapply(seq_along(keep.snps), #BPPARAM=bpparam, BPOPTIONS=bpoptions(stop.on.error = FALSE),
                              FUN=function(i, fixed, in.snps, solver, in.pca, meta.df, in.geno,
	    		                   big.mnn, reddim, grm, BPparam){

        message("Setting up model matrix for ", in.snps[i])
    	fixed.formula <- paste(c(unlist(strsplit(opt$fixed, split=",", fixed=TRUE)), paste0("chr", in.snps[i])), collapse=" + ")

        glmm.formula <- paste0("~ ", fixed.formula)
    	message("Model formula: ", glmm.formula)

        message("Running GLMM with ", opt$solver, " solver")
	i.geno <- data.frame("SNP"=in.geno[in.snps[i], ])
	i.geno[, opt$sampleID] <- colnames(in.geno)
	colnames(i.geno) <- c(paste0("chr", in.snps[i]), opt$sampleID)
    
        i.meta.df <- merge(meta.df, i.geno, by=opt$sampleID)
	rownames(i.meta.df) <- i.meta.df[, opt$sampleID]
	# drop NA genotypes
	i.meta.df <- i.meta.df[!is.na(i.meta.df[, paste0("chr", in.snps[i])]), ]

        keep.samps <- intersect(rownames(i.meta.df), rownames(grm))
    
        glmm.res <- testNhoods(x=big.mnn, design=as.formula(glmm.formula),
    	     	               design.df=i.meta.df[keep.samps, ], reduced.dim=opt$reddim,
                               kinship=grm[keep.samps, keep.samps], BPPARAM=BPparam,
                               fdr.weighting="graph-overlap", glmm.solver=opt$solver)
        print(warnings())
	print(dim(glmm.res))
	glmm.res$SNP <- in.snps[i]
	glmm.res$N <- length(keep.samps)
	return(glmm.res)
                                  }, fixed=opt$fixed, in.snps=keep.snps, solver=opt$solver, BPparam=bpparam,
			             in.geno=in.geno, meta.df=meta.df, big.mnn=big.mnn, reddim=opt$reddim, grm=grm)
} else{
    glmm.snp.list <- list()
    # PCA first by SVD of GRM
    set.seed(42)
    grm.svd <- svdr(x=grm, k=10)
    grm.d <- matrix(0L, ncol=10, nrow=10)
    diag(grm.d) <- grm.svd$d
    in.pca <- grm.svd$u %*% grm.d
    colnames(in.pca) <- paste0("PC", c(1:10))
    rownames(in.pca) <- colnames(grm)

    set.seed(NULL)
    # write this for bplapply?
    glmm.snp.list <- lapply(seq_along(keep.snps), #BPPARAM=bpparam, BPOPTIONS=bpoptions(stop.on.error = FALSE),
                              FUN=function(i, fixed, in.snps, pca, n.pcs, meta.df,
	    		                   big.mnn, reddim, BPparam){

        message("Setting up model matrix for ", in.snps[i])
    	fixed.formula <- paste(c(unlist(strsplit(opt$fixed, split=",", fixed=TRUE)),
	                         paste0("PC", c(1:10)), paste0("chr", in.snps[i])), collapse=" + ")

        glmm.formula <- paste0("~ ", fixed.formula)
    	message("Model formula: ", glmm.formula)

        message("Running GLM")
	i.geno <- data.frame("SNP"=in.geno[in.snps[i], ])
	i.geno[, opt$sampleID] <- colnames(in.geno)
	colnames(i.geno) <- c(paste0("chr", in.snps[i]), opt$sampleID)
    
        i.meta.df <- merge(meta.df, i.geno, by=opt$sampleID)
	rownames(i.meta.df) <- i.meta.df[, opt$sampleID]
	# drop NA genotypes
	i.meta.df <- i.meta.df[!is.na(i.meta.df[, paste0("chr", in.snps[i])]), ]
	keep.samps <- intersect(rownames(i.meta.df), rownames(in.pca))
	i.meta.df <- cbind.data.frame(i.meta.df[keep.samps, ], in.pca[keep.samps, ])
    
        glmm.res <- testNhoods(x=big.mnn, design=as.formula(glmm.formula),
    	     	               design.df=i.meta.df[keep.samps, ], reduced.dim=opt$reddim,
                               BPPARAM=BPparam,
                               fdr.weighting="graph-overlap", glmm.solver=opt$solver)
        print(warnings())
	print(dim(glmm.res))
	glmm.res$SNP <- in.snps[i]
	glmm.res$N <- length(keep.samps)
	return(glmm.res)
                                  }, fixed=opt$fixed, in.snps=keep.snps, BPparam=bpparam,
			             pca=in.pca, n.pcs=10, meta.df=meta.df, big.mnn=big.mnn, reddim=opt$reddim)
}

if(opt$glm){
    out.file <- paste0(opt$output, "_", snpchunk, "_GLM_results.tsv")
} else{
    out.file <- paste0(opt$output, "_", snpchunk, "_GLMM_results.tsv")
}
warnings()

full.res <- do.call(rbind.data.frame, glmm.snp.list)

write.table(full.res,
            file=out.file, sep="\t", quote=FALSE, row.names=FALSE)
message("All done")
