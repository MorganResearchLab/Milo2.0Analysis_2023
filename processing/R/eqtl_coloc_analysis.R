#! /usr/bin/env Rscript

## Perform colocalisation analysis using coloc
library(coloc)
library(reshape2)
library(readr)
library(miloR)
library(SingleCellExperiment)
library(Matrix)
library(cowplot)
library(dplyr)
library(rtracklayer)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCE"), type="character",
                     help="Path to Milo object")

parser <- add_option(parser, c("-d", "--DA"), type="character",
                     help="Path to popDA Milo results")

parser <- add_option(parser, c("-g", "--GWASDir"), type="character",
                     help="Directory containing Milo-GWAS results")

parser <- add_option(parser, c("-i", "--window"), type="numeric",
                     default=10000, help="Window size around lead SNP to extract")

parser <- add_option(parser, c("-e", "--eqtl"), type="character",
       	             help="Results table from Matrix eQTL")

parser <- add_option(parser, c("-t", "--sds"), type="character",
                     help="RDS file containing list of gene expression standard deviations")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output file prefix for DA testing results")

opt <- parse_args(parser)

message("Reading in Milo data: ", opt$SCE)
rand.ni.milo <- readRDS(opt$SCE)
rand.ni.meta <- as.data.frame(colData(rand.ni.milo))
rand.ni.meta$CellID <- colnames(rand.ni.milo)
rand.ni.fr <- as.data.frame(reducedDim(rand.ni.milo, "UMAP"))
rand.ni.fr$CellID <- colnames(rand.ni.milo)
rand.ni.meta <- merge(rand.ni.meta, rand.ni.fr, by="CellID")

message("Computing nhood standard deviations")
rand.ni.sd <- apply(nhoodCounts(rand.ni.milo), 1, sd) # these are numbered sequentially

message("Reading in Milo popDA results")
glmm.da.res <- read.table(opt$DA,
                          sep="\t", header=TRUE, stringsAsFactors = FALSE)

message("Reading in summary statistics from Milo-GWAS")
# gwas results files - this contains results for all nhoods for each lead SNP +/- window/2
gwas.dir <- opt$GWASDir
gwas.files <- list.files(gwas.dir, pattern="\\.gz$", full.names=TRUE, recursive=TRUE, include.dirs=TRUE)

message("Found ", length(gwas.files), " GWAS files to collate")
lead.snps <- c()
for(x in seq_along(gwas.files)){
  x.file <- gwas.files[x]
  x.lead <- gsub(x.file, pattern="([[:ascii:]]*)/([0-9]+)_([0-9]+)_([ATCG]+)_([ATCG]+)(_locus\\.tsv\\.gz)", 
                 replacement="\\2_\\3_\\4_\\5", perl=TRUE)
  x.chr <- as.numeric(gsub(x.lead, pattern="([0-9]+)_([0-9]+)_([ATCG]+)_([ATCG]+)", replacement="\\1"))
  x.bp <- as.numeric(gsub(x.lead, pattern="([0-9]+)_([0-9]+)_([ATCG]+)_([ATCG]+)", replacement="\\2"))
  
  lead.snps <- c(lead.snps, x.lead)
}

lead.snps <- unique(lead.snps)

message("Reading in gene expression standard deviations")
nh.exprs.sd <- readRDS(opt$sds)

message("Reading in Matrix eQTL results table")
matrix.qtl.df <- read.table(opt$eqtl, sep="\t", header=TRUE, stringsAsFactors=FALSE)
eqtl.snps <- intersect(lead.snps, unique(matrix.qtl.df$LeadSNP)) # this will remove the csQTL SNPs that don't have any proximal genes within 1Mb

##################################################
### Read in the association summary statistics ###
##################################################
# run a system call to tabix for a particular SNP
# split this into the summary statistics for each gene

eqtl.csqtl.coloc.list <- list()

# eqtl.snps are the ones that are _not_ in gene deserts
for(x in seq_along(eqtl.snps)){
  x.snp <- eqtl.snps[x]
  message("csQTL SNP: ", x.snp)
  x.eqtls <- matrix.qtl.df[matrix.qtl.df$LeadSNP %in% x.snp, ]
  x.eqtls$BP <- as.numeric(gsub(x.eqtls$snps, pattern="([0-9]+)_([0-9]+)_([ATCG]+)_([ATCG]+)", 
                                replacement="\\2"))
  x.egenes <- unique(x.eqtls[x.eqtls$FDR < 0.01, ]$gene)
  x.nhoods <- unique(x.eqtls[x.eqtls$FDR < 0.01, ]$Nhood)
  x.esnps <- unique(x.eqtls$snps)
  
  if(length(x.egenes) & length(x.nhoods)){
    x.chr <- gsub(x.snp, pattern="([0-9]+)_([0-9]+)_([ATCG]+)_([ATCG]+)", 
                  replacement="\\1")
    x.bp <- as.numeric(gsub(x.snp, pattern="([0-9]+)_([0-9]+)_([ATCG]+)_([ATCG]+)", 
                            replacement="\\2"))
    x.csqtl.file <- gwas.files[grepl(gwas.files, pattern=x.snp)] # these are from the first code block in this notebook
    x.tabix.query <- paste0("tabix ", x.csqtl.file, " ", x.chr, ":", x.bp-opt$window, "-", x.bp+opt$window) # interval size is also from above
    x.csqtl.df <- do.call(rbind.data.frame, lapply(system(x.tabix.query, intern=TRUE), FUN=function(TX){
      unlist(strsplit(TX, split="\t", fixed=TRUE))
      }))
    
    colnames(x.csqtl.df) <- c("SNP", "BP", "CHR", "logFC", "SE", "GeneticVariance", "SpatialFDR", "Nhood", "ALT", "REF")
    x.csqtl.df$CHR <- as.numeric(x.csqtl.df$CHR)
    x.csqtl.df$BP <- as.numeric(x.csqtl.df$BP)
    x.csqtl.df$logFC <- as.numeric(x.csqtl.df$logFC)
    x.csqtl.df$SE <- as.numeric(x.csqtl.df$SE)
    x.csqtl.df$GeneticVariance <- as.numeric(x.csqtl.df$GeneticVariance)
    x.csqtl.df$SpatialFDR <- as.numeric(x.csqtl.df$SpatialFDR)
    x.csqtl.df$Nhood <- as.numeric(x.csqtl.df$Nhood)
    x.csqtl.df <- x.csqtl.df[x.csqtl.df$SNP %in% x.esnps, ]
    
    for(q in seq_along(x.nhoods)){
      q.nh <- x.nhoods[q]
      message("csQTL Nhood: ", q.nh)
      q.csqtl <- x.csqtl.df[x.csqtl.df$Nhood %in% q.nh, ]
      
      x.coloc.list <- list("beta"=q.csqtl$logFC, "varbeta"=(q.csqtl$SE)**2,
                           "snp"=q.csqtl$SNP, "position"=q.csqtl$BP,
                           "type"="quant", "sdY"=rand.ni.sd[q.nh])
      
      for(k in seq_along(x.egenes)){
        k.gene <- x.egenes[k]
        message("eQTL gene: ", k.gene)
        k.sd <- (nh.exprs.sd[[paste0(q.nh, "_", x.snp)]])[k.gene]
        k.eqtls <- x.eqtls[x.eqtls$gene %in% k.gene & x.eqtls$snps %in% q.csqtl$SNP, ]
        k.eqtls <- k.eqtls[!duplicated(k.eqtls$snps), ]
        
        k.eqtl.coloc.list <- list("beta"=k.eqtls$beta,
                                  "varbeta"=(k.eqtls$beta/k.eqtls$statistic)**2,
                                  "snp"=k.eqtls$snps,
                                  "position"=k.eqtls$BP,
                                  "type"="quant", "sdY"=k.sd)
        
        # check coloc validity
        if(is.null(check_dataset(k.eqtl.coloc.list)) & is.null(check_dataset(x.coloc.list))){
	  message("Performing coloc analysis")
          k.coloc <- coloc.abf(x.coloc.list, k.eqtl.coloc.list)
          k.df <- do.call(cbind.data.frame, as.list(k.coloc$summary))  
          k.df$Nhood <- q.nh
          k.df$Gene <- k.gene
          k.df$LeadSNP <- x.snp
          
          k.ppsum <- sum(k.df[c(2:6)])
          if(k.ppsum > 1){
            warning("PP sum > 1 for ", k.gene, " in nhood ", q.nh, " for SNP ", x.snp)
          }
          
          eqtl.csqtl.coloc.list[[paste0(x.snp, "_", q.nh, "_", k.gene)]] <- k.df
          
        }
        
      }
    }
  }
}

csqtl.eqtl.coloc.df <- do.call(rbind.data.frame, eqtl.csqtl.coloc.list)
# merge with Randolph cell type annotations
csqtl.eqtl.coloc.df <- merge(csqtl.eqtl.coloc.df, glmm.da.res[, c("Nhood", "NhoodGroup", "ident", "ident_fraction")], by="Nhood", all.x=TRUE)


ofile <- paste0(opt$output, "_MatrixQTLcoloc.tsv")
message("Writing coloc results to: ", ofile)
write.table(csqtl.eqtl.coloc.df, file=ofile,
            sep="\t", quote=FALSE, row.names=FALSE)

message("All done")






