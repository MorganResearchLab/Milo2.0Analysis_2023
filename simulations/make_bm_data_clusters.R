## Make simple simulated DA on clusters ##

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(SingleCellExperiment)
  library(scran)
})

source('./simulations/benchmark_utils.R')

parser <- OptionParser()
parser <- add_option(parser, "--data_RDS", type="character",
                    help = "path to RDS storing SingleCellExperiment object")
parser <- add_option(parser, "--batch_seed", type="integer",
                    help = "Seed 4 batch effect")
parser <- add_option(parser, "--population", type="character",
                    help = "Which cell type is DA?")
parser <- add_option(parser, "--pop_enrichment", type="double", default=0.7,
                    help = "Max condition probability in DA population")
# parser <- add_option(parser, "--max_size", type="integer", default=1000,
#                     help = "Min number of cells in population to select")
parser <- add_option(parser, "--k", type="integer", default=50,
                    help = "K parameter")
parser <- add_option(parser, "--data_id", type="character", default="embryo",
                    help = "ID for the dataset used")
parser <- add_option(parser, "--make_batch_effect", type="character", default="yes",
                    help = "should synthetic batch effects be added? (yes/no)")
parser <- add_option(parser, "--output", type="character",
       	             help="Output file path")
opt <- parse_args(parser)

seed <- opt$batch_seed
data_path <- opt$data_RDS
k <- opt$k
pop <- opt$population
pop_enr <- opt$pop_enrichment
data_id <- opt$data_id
make_batch_effects <- opt$make_batch_effect

## Load data
print("Loading dataset...")
sce <- readRDS(data_path)

if(class(sce) == "list"){
  meta.df <- sce$meta
  sce <- sce$mylo
  colData(sce)$Block <- meta.df[colnames(sce), ]$Block
}

## Select population to simulate DA by size and save
# sized_pops = names(table(sce$celltype))[table(sce$celltype) < max_size]
# pop = sample(sized_pops, 1)

if (str_detect(pop, "_")) {
  pop <- str_replace(pop, "_", " ")
}

het.p <- 0.5 # additive model
sce <- add_synthetic_labels_by_cluster_genotype(sce, pop=paste0("B", pop), pop_column = "Block", seed=seed, pop_enr=pop_enr,
                                                het=het.p)
print(table(colData(sce)$Genotype1_prob))

## Consider TRUE DA effect size +- 10%
cond.probs <- c(1+pop_enr, (1 + pop_enr)*het.p, 1/3)
cond.probs <- cond.probs/sum(cond.probs)

if (pop_enr < 0.5) {
  t.prob <- min(cond.probs)
  da_lower <- t.prob - (t.prob/100) * 10
  da_upper <- t.prob + (t.prob/100) * 10
} else {
  t.prob <- max(cond.probs)
  da_lower <- t.prob - (t.prob/100) * 10
  da_upper <- t.prob + (t.prob/100) * 10
}


true_labels <- ifelse(sce$Genotype1_prob > da_lower & pop_enr >= 0.5, "PosLFC", 
                      ifelse(sce$Genotype1_prob < da_upper & pop_enr <= 0.5, "NegLFC", "NotDA"))

colData(sce)[["true_labels"]] <- true_labels
if (str_detect(pop, " ")) {
  pop <- str_replace(pop, " ", "_")
}

## Save coldata
outprefix <- str_c("benchmark_", data_id, "_pop_", pop, '_enr', pop_enr, "_seed", seed)
coldata <- data.frame(colData(sce)) %>% rownames_to_column()
ofile <- paste(opt$output, paste0(outprefix, ".coldata.csv"), sep="/")
write_csv(coldata, ofile)

## Simulate batch effects of different magnitude
set.seed(seed)
if (make_batch_effects=="yes") {
  print("Simulating batch effects...")
  bm_sce_ls <- lapply(c(0, 0.25, 0.5, 0.75, 1), function(sd){
    sce_be <- add_batch_effect(sce, batch_col = "synth_batches", norm_sd=sd)
    
    X_pca <- reducedDim(sce_be, "pca_batch")
    
    ## Save reduced dims
    write_csv(as.data.frame(X_pca) %>% rownames_to_column(), str_c(opt$output, outprefix, "_batchEffect", sd, ".pca.csv"))
  })
} else {
  X_pca <- reducedDim(sce, "PCA")
  ## Save reduced dims
  write_csv(as.data.frame(X_pca) %>% rownames_to_column(), str_c(opt$output, outprefix, "_batchEffect0.pca.csv"))
}
