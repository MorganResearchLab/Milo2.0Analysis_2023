#! /usr/bin/env Rscript

## Perform DA testing with Milo & edgeR with shuffled labels

library(SingleCellExperiment)
library(miloR)
library(edgeR)
library(BiocParallel)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-m", "--Milo"), type="character",
       	             help="The path to the Milo object")

parser <- add_option(parser, c("-a", "--Annotation"), type="character",
                     help="Column in colData to use for annotating DA testing results")

parser <- add_option(parser, c("-f", "--fdr"), type="numeric", default=0.1,
                     help="FDR threshold used for defining DA nhoods")

parser <- add_option(parser, c("-r", "--random"), type="character",
                     help="Variables in meta-data to use as random effects")

parser <- add_option(parser, c("-x", "--fixed"), type="character",
                     help="Variables in meta-data to use as fixed effects")

parser <- add_option(parser, c("-s", "--solver"), type="character", default="Fisher",
                     help="GLMM Solver to use: Fisher, HE or HE-NNLS")

parser <- add_option(parser, c("-b", "--nboots"), type="numeric", default=100,
                     help="The number of random permutations to use for label shuffling")

parser <- add_option(parser, c("-d", "--reddim"), type="character", default="PCA",
                     help="Reduced dimension in which graph was built")

parser <- add_option(parser, c("-l", "--sampleID"), type="character",
       	  	     help="Sample column used to identify individual observations")

parser <- add_option(parser, c("-g", "--glm"), action="store_true", default=FALSE,
                     help="Run with a GLM instead of GLMM - ignore random effects")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output file prefix for DA testing results")

opt <- parse_args(parser)

# set BPPARAM
bpparam <- MulticoreParam(workers=6)
#pparam <- SerialParam()
register(bpparam)
set.seed(42)

message(paste0("Reading in Milo object: ", opt$Milo))
big.mnn <- readRDS(opt$Milo)

## make the meta data from the colData
meta.df <- colData(big.mnn)

message("Making meta-data")
meta.df <- meta.df[!duplicated(meta.df[, opt$sampleID]), ]
rownames(meta.df) <- meta.df[, opt$sampleID]

message("Setting up model matrix")
all.fe <- unlist(strsplit(opt$fixed, split=",", fixed=TRUE))
test.fe <- all.fe[length(all.fe)]
fixed.formula <- paste(unlist(strsplit(opt$fixed, split=",", fixed=TRUE)), collapse=" + ")
rand.formula <- paste(unlist(lapply(strsplit(opt$random, split=",", fixed=TRUE), FUN=function(RX) paste0("(1|", RX, ")"))),
                      collapse=" + ")

if(opt$glm){
    glmm.formula <- paste0("~ ", fixed.formula)
    message("GLM formula: ", glmm.formula)
} else{
    glmm.formula <- paste0("~ ", paste(c(fixed.formula, rand.formula), collapse=" + "))
    message("GLMM formula: ", glmm.formula)
}

message("Randomly permuting fixed effect labels ", opt$nboots, " times")
rando.list <- list()
for(i in seq_len(opt$nboots)){
    message("Shuffling labels of ",test.fe)
    boot.df <- meta.df
    boot.df[, test.fe] <- sample(meta.df[, test.fe])
    print(table(meta.df[, test.fe]))
    print(table(boot.df[, test.fe]))

    message("Running GLMM with ", opt$solver, " solver")
    glmm.res <- testNhoods(x=big.mnn, design=as.formula(glmm.formula), design.df=boot.df, reduced.dim=opt$reddim,
                           fdr.weighting="graph-overlap", glmm.solver=opt$solver, BPPARAM=bpparam)
    glmm.res$NBoot <- i
    print(warnings())
    print(dim(glmm.res))
    rando.list[[i]] <- glmm.res
}

all.glmm.res <- do.call(rbind.data.frame, rando.list)

if(opt$glm){
    out.file <- paste0(opt$output, "shuffled_GLM_results.tsv")
} else{
    out.file <- paste0(opt$output, "shuffled_GLMM_results.tsv")
}
write.table(all.glmm.res,
	    file=out.file, sep="\t", quote=FALSE, row.names=FALSE)
	    