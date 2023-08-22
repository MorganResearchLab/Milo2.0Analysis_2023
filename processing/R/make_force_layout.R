#! /usr/bin/env Rscript

## Given an input dataset make a force-directed layout from a KNN graph

library(SingleCellExperiment)
library(scran)
library(scuttle)
library(scater)
library(optparse)
library(Matrix)
library(igraph)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCE"), type="character",
       	             help="The path to the combined SCE object")

parser <- add_option(parser, c("-d", "--dimensions"), type="numeric",
                     help="The number of dimensions to use for graph building")

parser <- add_option(parser, c("-r", "--reduced"), type="character",
                     help="The reducedDim slot to use for graph building")

parser <- add_option(parser, c("-k", "--knn"), type="numeric", default=21,
       	             help="k value to use for graph building")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for output files")

parser <- add_option(parser, c("-b", "--barcodes"), type="character",
                     default=FALSE, help="List of barcodes to use for graph building")

parser <- add_option(parser, c("-l", "--layout"), type="character", default="FR",
                     help="Force-direct graph layout to use. Options are FR or OPT")

opt <- parse_args(parser)

message(paste0("Reading in SCE object: ", opt$SCE))
big.mnn <- readRDS(opt$SCE)

if(!opt$barcodes %in% c("FALSE")){
   keep.bcs <- read.table(opt$barcodes, sep="\t", header=FALSE, stringsAsFactors=FALSE)[, 1]
   message(paste0("Subsetting to ", length(keep.bcs), " barcodes"))
   big.mnn <- big.mnn[, intersect(colnames(big.mnn), keep.bcs)]
}

message("Computing a KNN graph with ", opt$knn, " neighbours using ", opt$dimensions, " dimensions")
nn.graph <- buildKNNGraph(reducedDim(big.mnn, opt$reduced)[, c(1:opt$dimensions)], transposed=TRUE, k=opt$knn)

message("Computing a ", opt$layout, " force-directed layout")
set.seed(42)
if(opt$layout %in% c("FR")){
    graph.layout <- as.data.frame(layout_with_fr(nn.graph))
    colnames(graph.layout) <- c("FR1", "FR2")
    graph.layout$CellID <- colnames(big.mnn)    
} else if(opt$layout %in% c("OPT")){
    graph.layout <- layout_with_graphopt(nn.graph)
    colnames(graph.layout) <- c("OPT1", "OPT2")
    graph.layout$CellID <- colnames(big.mnn) 
} else{
    stop(opt$layout, " not recognised. Should be either FR or OPT")
}

message("Saving graph layout")
umap.file <- paste0(opt$output, "_FDL.tsv")
write.table(graph.layout, file=umap.file, sep="\t", quote=FALSE, row.names=FALSE)

## add to SCE object
reducedDim(big.mnn, opt$layout) <- as.matrix(graph.layout[, c(1:2)])
sce.out <- paste0(opt$output, "_SCE-FDL.RDS")
saveRDS(big.mnn, file=sce.out)