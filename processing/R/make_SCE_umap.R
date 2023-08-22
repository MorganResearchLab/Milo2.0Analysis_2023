#! /usr/bin/env Rscript

## Given an input dataset make a UMAP with certain numbers of dimensions

library(SingleCellExperiment)
library(scran)
library(optparse)
library(Matrix)
library(umap)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCE"), type="character",
       	             help="The path to the combined SCE object")

parser <- add_option(parser, c("-d", "--dimensions"), type="numeric",
                     help="The number of dimensions to use for UMAP")

parser <- add_option(parser, c("-r", "--reduced"), type="character",
                     help="The reducedDim slot to use for UMAP")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for output files")

parser <- add_option(parser, c("-b", "--barcodes"), type="character",
                     default=FALSE, help="List of barcodes to use for UMAP")

opt <- parse_args(parser)

if(opt$SCE %in% c("package")){
  library(scuttle)
  library(MouseThymusAgeing)
  big.mnn <- MouseDropletData()
} else{
  message(paste0("Reading in SCE object: ", opt$SCE))
  big.mnn <- readRDS(opt$SCE)
}

if(!opt$barcodes %in% c("FALSE")){
   keep.bcs <- read.table(opt$barcodes, sep="\t", header=FALSE, stringsAsFactors=FALSE)[, 1]
   message(paste0("Subsetting to ", length(keep.bcs), " barcodes"))
   big.mnn <- big.mnn[, intersect(colnames(big.mnn), keep.bcs)]
}

message("Computing a new UMAP in the corrected space")
set.seed(42)
mnn.umap.map <- as.data.frame(umap(as.matrix(reducedDim(big.mnn, opt$reduced)[, c(1:opt$dimensions)]),
                                   n_neighbors=30,
                                   init='random',
                                   n_components=2,
                                   min_dist=0.2,
                                   method='naive')$layout)
colnames(mnn.umap.map) <- c("UMAP1", "UMAP2")
rownames(mnn.umap.map) <- colnames(big.mnn)
mnn.umap.map$CellID <- colnames(big.mnn)

message("Saving UMAP")
umap.file <- paste0(opt$output, "_UMAPs.tsv")
write.table(mnn.umap.map, file=umap.file, sep="\t", quote=FALSE, row.names=FALSE)

## add to SCE object
reducedDim(big.mnn, "UMAP") <- as.matrix(mnn.umap.map[, c("UMAP1", "UMAP2")])
sce.out <- paste0(opt$output, "_SCE-UMAP.RDS")
saveRDS(big.mnn, file=sce.out)