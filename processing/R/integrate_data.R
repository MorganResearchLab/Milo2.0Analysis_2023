#! /usr/bin/env Rscript

########################################################
## Combine 2 data sets with the purpose of performing a
## label transfer from a reference to a query
########################################################

library(SingleCellExperiment)
library(optparse)
library(batchelor)
library(BiocNeighbors)
library(scran)
library(irlba)
library(FNN)

parser <- OptionParser()
parser <- add_option(parser, c("-x", "--SCEList"), type="character",
       	             help="A comma separated list of SCE file paths")

parser <- add_option(parser, c("-d", "--dimensions"), type="numeric",
                     help="Number of dimensions to use for integration")

parser <- add_option(parser, c("-k", "--knn"), type="numeric",
                     help="K for MNN searching")

parser <- add_option(parser, c("-l", "--l2norm"), action="store_true",
                     help="Perform cosine normalisation before PCA and integration")

parser <- add_option(parser, c("-b", "--batch"), type="character",
                     help="Column of colData to use for batch correction")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for the combined SCE object")

opt <- parse_args(parser)

sce.files <- unlist(strsplit(opt$SCEList, split=",", fixed=TRUE))

if(length(sce.files) > 1){
    message(paste0("Reading SCE files: ", sce.files))
    sce.list <- lapply(sce.files, readRDS)
    sce.names <- unlist(lapply(sce.files, FUN=function(NAME){
                              x.name <- unlist(strsplit(NAME, split="/", fixed=TRUE))
     		              x.name <- gsub(x.name[length(x.name)], pattern="_SCE\\.RDS", replacement="")
			      return(x.name)
                              }))
    names(sce.list) <- sce.names
} else{
    # split the SCE by the batch columns
    orig.sce <- readRDS(sce.files)
    sce.names <- unique(colData(orig.sce)[, opt$batch])
    sce.list <- list()
    for(x in seq_along(sce.names)){
        x.sce <- orig.sce[, colData(orig.sce)[, opt$batch] %in% sce.names[x]]
	sce.list[[sce.names[x]]] <- x.sce
    }
    names(sce.list) <- sce.names
}

message("Finding highly variable genes")
hvg.list <- sapply(sce.names, FUN=function(XCE){
    x.hvg <- modelGeneVar(sce.list[[XCE]], assay.type="logcounts")
    x.hvg$FDR[is.na(x.hvg$FDR)] <- 1
    x.hvg$GeneID <- rownames(x.hvg)
    message("Found ", sum(x.hvg$FDR < 0.1), " HVGs")
    x.hvg$GeneID[x.hvg$FDR < 0.1]
    })

if(length(hvg.list) > 1){
    common.hvgs <- Reduce(x=hvg.list, f=function(X, Y) union(X, Y))
    all.genes <- Reduce(x=lapply(sce.list, rownames), f=function(x, y) intersect(x, y))    
} else{
    common.hvgs <- hvg.list[[1]]
    all.genes <- rownames(sce.list[[1]])
}

common.hvgs <- intersect(common.hvgs, all.genes)
message(paste0("Found ", length(common.hvgs), " common HVGs"))

red.dim <- "PCA"
if(opt$l2norm){
  message("Performing cosine normalisation")
  for(x in seq_along(sce.names)){
      x.sce <- sce.list[[sce.names[x]]][all.genes, ]
      assay(x.sce, "Cosine") <- cosineNorm(logcounts(x.sce))
      rowData(x.sce) <- NULL
      sce.list[[sce.names[x]]] <- x.sce
    }
    
  print(lapply(sce.list, dim))
  print(lapply(sce.list, FUN=function(SX) dim(SX[common.hvgs, ])))
  message("Performing joint PCA on " , opt$dimensions, " dimensions")
  all.pca <- multiBatchPCA(sce.list, d=opt$dimensions, subset.row=common.hvgs,
                           assay.type="Cosine")
  names(all.pca) <- sce.names
  red.dim <- "Cosine.PCA"
} else{
  message("Performing joint PCA on " , opt$dimensions, " dimensions")
  all.pca <- multiBatchPCA(sce.list, d=opt$dimensions, subset.rows=common.hvgs,
                           assay.type="logcounts")
  names(all.pca) <- sce.names
}

for(i in seq_along(sce.names)){
    i.sce <- sce.list[[sce.names[i]]]
    reducedDim(i.sce, red.dim) <- all.pca[[sce.names[i]]][, c(1:opt$dimensions)]
    colData(i.sce)$Batch <- sce.names[i]
    sce.list[[sce.names[i]]] <- i.sce
}

message(paste0("Performing integration using ", opt$dimensions, " dimensions"))
mnn.out <- reducedMNN(all.pca, 
		      k=opt$knn)

# find the common column names
common.meta <- Reduce(x=lapply(sce.list, FUN=function(XC) colnames(colData(XC))),
                      f=function(x, y) intersect(x, y))

for(j in seq_along(sce.names)){
    j.sce <- sce.list[[sce.names[j]]]
    j.df <- colData(j.sce)[, common.meta]
    colData(j.sce) <- NULL
    for(q in seq_along(common.meta)){
      colData(j.sce)[, common.meta[q]] <- j.df[, common.meta[q]]
    }
    print(head(j.sce))
    sce.list[[sce.names[j]]] <- j.sce
}

message("Making combined SCE object")
big.logcounts <- do.call(cbind, lapply(sce.list, logcounts))
big.counts <- do.call(cbind, lapply(sce.list, counts))
big.meta <- do.call(rbind, lapply(sce.list, colData))
big.sce <- SingleCellExperiment(assay=list(logcounts=big.logcounts, counts=big.counts), colData=big.meta)

reducedDim(big.sce, red.dim) <- do.call(rbind, all.pca)
reducedDim(big.sce, gsub(red.dim, pattern="PCA", replacement="MNN")) <- mnn.out$corrected

ofile <- paste0(opt$output, "_SCE.RDS")
message(paste0("Saving combined SCE to ", ofile))
saveRDS(big.sce, file=ofile)

message("All done")