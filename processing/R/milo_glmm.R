#! /usr/bin/env Rscript

## Perform DA testing with Milo & edgeR

library(SingleCellExperiment)
library(miloR)
library(edgeR)

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

parser <- add_option(parser, c("-d", "--reddim"), type="character", default="PCA",
                     help="Reduced dimension in which graph was built")

parser <- add_option(parser, c("-l", "--sampleID"), type="character",
       	  	     help="Sample column used to identify individual observations")

parser <- add_option(parser, c("-g", "--glm"), action="store_true", default=FALSE,
                     help="Run with a GLM instead of GLMM - ignore random effects")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output file prefix for DA testing results")

opt <- parse_args(parser)

message(paste0("Reading in Milo object: ", opt$Milo))
big.mnn <- readRDS(opt$Milo)

## make the meta data from the colData
meta.df <- colData(big.mnn)

message("Making meta-data")
meta.df <- meta.df[!duplicated(meta.df[, opt$sampleID]), ]
rownames(meta.df) <- meta.df[, opt$sampleID]

message("Setting up model matrix")
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

message("Running GLMM with ", opt$solver, " solver")
glmm.res <- testNhoods(x=big.mnn, design=as.formula(glmm.formula), design.df=meta.df, reduced.dim=opt$reddim,
                       fdr.weighting="graph-overlap", glmm.solver=opt$solver)
print(warnings())
print(dim(glmm.res))
message("Annotating linear results")
print(dim(nhoodCounts(big.mnn)))
test.res <- annotateNhoods(big.mnn, glmm.res, opt$Annotation)

if(opt$glm){
    out.file <- paste0(opt$output, "GLM_results.tsv")
} else{
    out.file <- paste0(opt$output, "GLMM_results.tsv")
}
write.table(glmm.res,
	    file=out.file, sep="\t", quote=FALSE, row.names=FALSE)
	    