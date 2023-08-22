#! /usr/bin/env Rscript

## Given an input dataset of 2 batches - downsample the batches to the same median sequencing depth

library(DropletUtils)
library(SingleCellExperiment)
library(scran)
library(scater)
library(irlba)
library(optparse)
library(Matrix)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCE"), type="character",
       	             help="The path to the combined SCE object")

parser <- add_option(parser, c("-d", "--dimensions"), type="numeric",
                     help="The number of dimensions to use for batch correction")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for output files denoised data and combined SCE object")

parser <- add_option(parser, c("-c", "--center"), action="store_true",
                     help="Flag to center data before PCA")

parser <- add_option(parser, c("-l", "--scale"), action="store_true",
                     help="Flag to scale data before PCA")

parser <- add_option(parser, c("-a", "--assay"), type="character",
                     help="Gene expression measurements to use for PCA")

parser <- add_option(parser, c("-u", "--subset"), type="character", default="FALSE",
                     help="Value on which to subset cells in the Milo object, e.g. selecting only one condition")

parser <- add_option(parser, c("-e", "--value"), type="character",
                     help="Column name of colData() used to subset with `subset`")

opt <- parse_args(parser)

message(paste0("Reading in ", opt$SCE))
big.sce <- readRDS(opt$SCE)

if(!opt$subset %in% c("FALSE")){
    message("Subsetting cells using ", opt$subset, " and ", opt$value)
    big.sce <- big.sce[, colData(big.sce)[, opt$value] %in% opt$subset]
    message("Keeping ", ncol(big.sce), " subset cells")
}

message("Computing highly variable genes")
hvg.stats <- modelGeneVar(big.sce)
hvg.stats$FDR[is.na(hvg.stats$FDR)] <- 1
hvg.stats$GeneID <- rownames(hvg.stats)
write.table(hvg.stats, paste0(opt$output, "_hvgs.txt"), sep="\t", quote=FALSE, row.names=FALSE)  

common.hvg <- hvg.stats$FDR < 0.1

if(sum(common.hvg) < 300){
  message("Only found ", sum(common.hvg), "HVGs - selecting top 1500 HVGs instead")
  top.hvg <- rownames((hvg.stats[order(hvg.stats$FDR, decreasing=FALSE), ])[1:1500, ])
  common.hvg <- rownames(hvg.stats) %in% top.hvg
}

## define based on HVGs
message(paste0("Performing joint PCA with ", sum(common.hvg), " HVGs"))
all.pca <- prcomp_irlba(t(assay(big.sce[common.hvg ,], opt$assay)), n=opt$dimensions + 1, .scale=opt$scale, center=opt$center)

out.pca <- as.data.frame(all.pca$x)
out.pca$CellID <- colnames(big.sce)
write.table(out.pca, file=paste0(opt$output, "_", opt$subset, "-PCA.tsv"), sep="\t", quote=FALSE)

reducedDim(big.sce, "PCA") <- all.pca$x

ofile.sce <- paste0(opt$output, "_", opt$subset, "-PCA.RDS")

message(paste0("Saving SCE to: ", ofile.sce))
saveRDS(big.sce, file=ofile.sce)

message("All done")
