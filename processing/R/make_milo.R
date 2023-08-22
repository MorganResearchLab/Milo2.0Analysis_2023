#! /usr/bin/env Rscript

### Build a Milo object ready for DA testing
library(SingleCellExperiment)
library(miloR)
library(optparse)
library(MouseThymusAgeing)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCE"), type="character",
       	             help="The path to the combined SCE object")

parser <- add_option(parser, c("-i", "--sample"), type="character",
                     help="Column of SCE colData corresponding to experimental sample")

parser <- add_option(parser, c("-r", "--remove"), action="store_true", default=FALSE,
                     help="Flag to remove doublets if stored in a 'Class' colData slot")

parser <- add_option(parser, c("-x", "--drop"), action="store_true", default=FALSE,
                     help="Flag to drop replicate 3 or not")

parser <- add_option(parser, c("-b", "--batch"), type="character", default=FALSE,
                     help="Column of colData() to combine with `sample` to create batchXsample nhood counts")

parser <- add_option(parser, c("-k", "--knn"), type="numeric",
                     help="Value to use for k-NN graph building")

parser <- add_option(parser, c("-p", "--props"), type="numeric",
                     help="Proportion of indices to sample for nhood computation")

parser <- add_option(parser, c("-d", "--dimensions"), type="numeric",
                     help="Number of dimensions to use for distance calculations")

parser <- add_option(parser, c("-n", "--reduced"), type="character",
                     help="Name of reduced dimension to use for graph building, etc")

parser <- add_option(parser, c("-l", "--overlap"), type="numeric", default=0,
                     help="Minimum cell overlap between nhoods for nhood adjacency")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output file for Milo object")

opt <- parse_args(parser)

if(opt$SCE %in% c("package")){
  library(MouseThymusAgeing)
  library(scran)
  library(scuttle)
  big.mnn <- MouseDropletData()
  big.mnn <- logNormCounts(big.mnn)
  
  keep.genes <- !is.na(rowData(big.mnn)$ensembl_gene_id)
  gene.names <- rowData(big.mnn)$ensembl_gene_id
  gene.names <- gene.names[keep.genes]
  rowData(big.mnn) <- NULL
  big.mnn <- big.mnn[keep.genes, ]
  rownames(big.mnn) <- gene.names
  
  colData(big.mnn)[, opt$sample] <- paste(big.mnn$HTO, gsub(big.mnn$CellOrigin, pattern="([0-9])([A-Za-z]{2})(Run[0-9])", replacement="R\\1"), sep="_")
} else{
  message(paste0("Reading in SCE object: ", opt$SCE))
  big.mnn <- readRDS(opt$SCE)
}

if(opt$remove){
  big.mnn <- big.mnn[, big.mnn$Doublet.Class %in% c("Singlet")]
  message(paste0("Keeping ", ncol(big.mnn), " singlets"))
}

if(opt$drop){
  message("Dropping replicate 3")
  big.mnn <- big.mnn[, !grepl(big.mnn$HTO, pattern="R3")]
}

if(!opt$batch %in% c("FALSE")){
  message("Adding batch info to sample columns")
  colData(big.mnn)$Batch.Sample <- paste(colData(big.mnn)[, opt$batch], colData(big.mnn)[, opt$sample], sep="_")
  samp.col <- "Batch.Sample"
} else{
  samp.col <- opt$sample
}


message("Making Milo object")
milo <- Milo(big.mnn)

message("Building kNN graph")
milo <- buildGraph(milo, k=opt$knn, d=opt$dimensions, reduced.dim=opt$reduced)

message("Compute representative neighbourhoods")
milo <- makeNhoods(milo, prop=opt$props, k=opt$knn, d=opt$dimensions,
     	           refined=TRUE, reduced_dims=opt$reduced, refinement_scheme="graph")

message("Count cells over neighbourhoods and samples")
print(samp.col)
print(head(colData(milo)))
milo <- countCells(milo, meta.data=as.data.frame(colData(milo)), samples=samp.col)

message("Building nhood graph")
milo <- buildNhoodGraph(milo, overlap=opt$overlap)

message("Calculating per-neighbourhood gene expression")
milo <- calcNhoodExpression(milo)

message(paste0("Saving Milo object to: ", opt$output))
saveRDS(milo, file=opt$output)