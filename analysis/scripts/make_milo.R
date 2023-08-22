#! /usr/bin/env Rscript

### Build a Milo object ready for DA testing
library(SingleCellExperiment)
library(miloR)
library(BiocParallel)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCE"), type="character",
       	             help="The path to the combined SCE object")

parser <- add_option(parser, c("-k", "--keep"), type="character", default="FALSE",
                      help="File of cell barcodes to keep")

parser <- add_option(parser, c("-x", "--drop"), type="character", default="FALSE",
                     help="File of cell barcodes to drop")

parser <- add_option(parser, c("-i", "--sample"), type="character",
                     help="Column of SCE colData corresponding to experimental sample")

parser <- add_option(parser, c("-r", "--remove"), action="store_true",
                     help="Flag to remove doublets if stored in a 'Class' colData slot")

parser <- add_option(parser, c("-n", "--knn"), type="numeric",
                     help="Value to use for k-NN graph building")

parser <- add_option(parser, c("-p", "--props"), type="numeric",
                     help="Proportion of indices to sample for nhood computation")

parser <- add_option(parser, c("-e", "--reducedDim"), type="character",
                     help="Reduced dimension slot to used. Default=PCA", default="PCA")

parser <- add_option(parser, c("-d", "--dimensions"), type="numeric",
                     help="Number of dimensions to use for distance calculations")

parser <- add_option(parser, c("-l", "--overlap"), type="numeric", default=0,
                     help="Minimum cell overlap between nhoods for nhood adjacency")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output file for Milo object")

opt <- parse_args(parser)

bpparam <- MulticoreParam(workers=6)
register(bpparam)
set.seed(42)

message(paste0("Reading in SCE object: ", opt$SCE))
big.mnn <- readRDS(opt$SCE)

if(any(names(opt) %in% c("remove"))){
    message(paste0("Removing ", sum(big.mnn$Class %in% c("Doublet")), " doublets"))
    big.mnn <- big.mnn[!big.mnn$Class %in% c("Doublet")]
}

if(!opt$keep %in% c("FALSE")){
  keep.bcs <- read.table(opt$keep, sep="\t", header=FALSE, stringsAsFactors=FALSE)[, 1]
  keep.bcs <- intersect(keep.bcs, colnames(big.mnn))

  message(paste("Keeping ", length(keep.bcs), " cell barcodes"))
  big.mnn <- big.mnn[, colnames(big.mnn) %in% keep.bcs]
}

if(!opt$drop %in% c("FALSE")){
  drop.bcs <- read.table(opt$drop, sep="\t", header=FALSE, stringsAsFactors=FALSE)[, 1]
  drop.bcs <- intersect(drop.bcs, colnames(big.mnn))

  message(paste("Dropping ", length(drop.bcs), " cell barcodes"))
  big.mnn <- big.mnn[, !colnames(big.mnn) %in% drop.bcs]
}

message("Making Milo object")
milo <- Milo(big.mnn)

message("Building kNN graph")
milo <- buildGraph(milo, k=opt$knn, d=opt$dimensions, reduced.dim=opt$reducedDim, BPPARAM=bpparam)

message("Compute representative neighbourhoods")
# why is the memory exploding here?
milo <- makeNhoods(milo, prop=opt$props, k=opt$knn, d=opt$dimensions, refined=TRUE,
                   reduced_dims=opt$reducedDim, refinement_scheme="graph")

message("Count cells over neighbourhoods and samples")
milo <- countCells(milo, meta.data=as.data.frame(colData(milo)), samples=opt$sample)

#message("Calculating nhood distances")
#milo <- calcNhoodDistance(milo, d=opt$dimensions, reduced.dim=opt$reducedDim)

message("Building nhood graph")
milo <- buildNhoodGraph(milo, overlap=opt$overlap)

message("Calculating per-neighbourhood gene expression")
milo <- calcNhoodExpression(milo)

message(paste0("Saving Milo object to: ", opt$output))
saveRDS(milo, file=opt$output)