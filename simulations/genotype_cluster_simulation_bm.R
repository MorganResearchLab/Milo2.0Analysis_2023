suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(miloR)
  library(ggplot2)
  library(uwot)
  library(optparse)
})

parser <- OptionParser()
parser <- add_option(parser, "--data_RDS", type="character",
                     help = "path to RDS storing SingleCellExperiment object")
parser <- add_option(parser, "--simulation", type="character",
                     help="path to file containing simulated conditions and ground truth")
parser <- add_option(parser, "--alpha", type="numeric",
                     help="SpatialFDR threshold to define DA nhoods")
parser <- add_option(parser, "--seed", type="numeric",
                     help="Random seed to use for UMAP")
parser <- add_option(parser, "--k", type="numeric",
                     help="Number of nearest-neighbours for graph and nhoods")
parser <- add_option(parser, "--d", type="numeric",
                     help="Number of reduced dimensions to use for graph building")
parser <- add_option(parser, "--prop", type="numeric",
                     help="Propotion of initial nhoods to sample before refinement")
parser <- add_option(parser, "--output", type="character",
                     help="Output file path")
opt <- parse_args(parser)


alpha <- opt$alpha
message("Reading in ", opt$data_RDS)
milo.sim <- readRDS(opt$data_RDS)
sim.coldata <- read.csv(opt$simulation,
                        sep=",", header=TRUE, stringsAsFactors = FALSE)
rownames(sim.coldata) <- sim.coldata$rowname
milo.obj <- milo.sim$mylo
meta.sim <- milo.sim$meta

colData(milo.obj)$Block <- meta.sim[colnames(milo.obj), ]$Block
colData(milo.obj)$Condition <- meta.sim[colnames(milo.obj), ]$Condition
colData(milo.obj)$Replicate <- meta.sim[colnames(milo.obj), ]$Replicate
colData(milo.obj)$Vertex <- meta.sim[colnames(milo.obj), ]$Vertex
colData(milo.obj)$Genotype <- sim.coldata[colnames(milo.obj), ]$synth_labels
colData(milo.obj)$Sample <- sim.coldata[colnames(milo.obj), ]$synth_samples

set.seed(opt$seed)
message("Building UMAP for visualisation")
# add umap
set.seed(opt$seed)
ump.df <- as.matrix(umap(reducedDim(milo.obj, "PCA"), min_dist=0.5,
                         n_neighbors=opt$k, n_components = 2, init="spectral"))
colnames(ump.df) <- c("UMAP1", "UMAP2")
reducedDim(milo.obj, "UMAP") <- ump.df

message("Starting Milo workflow")
milo.obj <- buildGraph(milo.obj, k=opt$k, d=opt$d, reduced.dim="PCA")
milo.obj <- makeNhoods(milo.obj, k=opt$k, d=opt$d, prop=opt$prop, refined=TRUE, refinement_scheme="graph")
milo.obj <- countCells(milo.obj, samples="Sample", meta.data=colData(milo.obj))
milo.obj <- buildNhoodGraph(milo.obj)
n.samples <- length(unique(colData(milo.obj)$Sample)) * 3
x.min <- min(c(n.samples, min(rowSums(nhoodCounts(milo.obj)))))

# plotNhoodSizeHist(milo.obj, bins=20) + 
#   expand_limits(x=c(x.min)) +
#   geom_vline(xintercept=n.samples, lty=2)

message("Performing DA testing")
test.meta <- dplyr::distinct(as.data.frame(colData(milo.obj)), Sample, .keep_all=TRUE)
rownames(test.meta) <- test.meta$Sample
test.meta$GenotypeOrd <- ifelse(test.meta$Genotype == "Genotype1", 2,
                                ifelse(test.meta$Genotype == "Genotype2", 1, 0))
test.model <- model.matrix(~GenotypeOrd, data=test.meta[colnames(nhoodCounts(milo.obj)), ])
#test.model <- test.model[, c(1:2)]

sim.res <- testNhoods(milo.obj, design=test.model, design.df=test.meta[colnames(nhoodCounts(milo.obj)), ],
                      fdr.weighting="graph-overlap")
message("Identified ", sum(table(sim.res$SpatialFDR < alpha)), " DA nhoods at SpatialFDR<", alpha)

# plotNhoodGraphDA(milo.obj, milo_res = sim.res, alpha=alpha)

sim.res <- annotateNhoods(milo.obj, sim.res, coldata_col="Block")
# plotDAbeeswarm(sim.res, group.by="Block", alpha=alpha)
write.table(sim.res,
            file=paste0(opt$output, "_DAres.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

message("Comparing results to ground truth")
# assign cells to DA label
cell.label <- as.data.frame(sapply(seq_len(nrow(nhoods(milo.obj))),
                                   FUN=function(PX, alpha.fdr) {
                                     px.nh <- which(nhoods(milo.obj)[PX, ] == 1)
                                     px.fdr <- sim.res[px.nh, ]$SpatialFDR < alpha.fdr
                                     px.lfc <- sim.res[px.nh, ]$logFC
                                     px.label <- unique(ifelse(px.fdr & px.lfc > 0, "PosLFC",
                                                               ifelse(px.fdr & px.lfc < 0, "NegLFC", "NotDA")))
                                     
                                     if(any(px.label != "NotDA")){
                                       if(any(px.label == "PosLFC") & any(px.label != "NegLFC")){
                                         px.label <- "PosLFC"
                                       } else if(any(px.label == "NegLFC") & any(px.label != "PosLFC")){
                                         px.label <- "NegLFC"
                                       } else if(any(px.label == "NegLFC") & any(px.label == "PosLFC")){
                                         px.label <- "NotDA"
                                       }
                                     } else{
                                       if(any(px.label == "PosLFC") & any(px.label != "NegLFC")){
                                         px.label <- "PosLFC"
                                       } else if(any(px.label == "NegLFC") & any(px.label != "PosLFC")){
                                         px.label <- "NegLFC"
                                       } else if(any(px.label == "NegLFC") & any(px.label == "PosLFC")){
                                         px.label <- "NotDA"
                                       }
                                     }
                                     
                                     if(length(px.label)){
                                       return(px.label)
                                     } else{
                                       return("None")
                                     }
                                     
                                   }, alpha.fdr=alpha))
colnames(cell.label) <- c("test_labels")
cell.label$rowname <- colnames(milo.obj)
sim.bench.df <- merge(sim.coldata, cell.label, by="rowname")
sim.bench.df$test_labels <- factor(sim.bench.df$test_labels,
                                   levels=c("PosLFC", "NotDA", "None", "NegLFC"))
sim.bench.df$true_labels <- factor(sim.bench.df$true_labels,
                                   levels=c("PosLFC", "NotDA", "None", "NegLFC"))

message("Computing confusion matrix")
true.vs.test <-  table(sim.bench.df$true_labels, sim.bench.df$test_labels)
tp <- true.vs.test["PosLFC", "PosLFC"] # + true.vs.test["NegLFC", "NegLFC"]
tn <- sum(true.vs.test[c("NotDA", "None"), c("NotDA", "None")])
fn <- sum(true.vs.test["PosLFC", c("NotDA", "None")])
fp <- sum(true.vs.test[c("NotDA", "None"), c("PosLFC", "NegLFC")])

confuse.mat <- matrix(c(tp, fn, fp, tn), ncol=2)
print(confuse.mat)
dimnames(confuse.mat) <- list(c("True", "False"), c("True", "False"))

write.table(confuse.mat,
            file=paste0(opt$output, "_confusion.tsv"),
            row.names=TRUE, sep="\t")

fdr <- fp/(fp+tp)
tpr <- tp/(fn+tp)
bench.mat <- as.data.frame(matrix(c(fdr, tpr), ncol=1))
bench.mat <- cbind(c("FDR", "TPR"), bench.mat)
colnames(bench.mat) <- c("Metric", "Value")

write.table(bench.mat,
            file=paste0(opt$output, "_FDR_TPR.tsv"),
            sep="\t", row.names=FALSE)



