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
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCE"), type="character",
                     help="Path to Milo object")

parser <- add_option(parser, c("-d", "--DA"), type="character",
                     help="Path to popDA Milo results")

parser <- add_option(parser, c("-g", "--GWASDir"), type="character",
                     help="Directory containing Milo-GWAS results")

parser <- add_option(parser, c("-e", "--eqtl"), type="character",
                     help="Tabix-indexed results file from eQTL study")

parser <- add_option(parser, c("-a", "--eqtlName"), type="character",
                     help="Name of eQTL study name")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output file prefix for DA testing results")

opt <- parse_args(parser)

eqtl.sampsize.list <- list("QTD000021"=200, "QTD000031"=200, "QTD000439"=91, "QTD000444"=91,
                           "QTD000449"=91, "QTD000454"=91, "QTD000464"=91, "QTD000469"=91,
                           "QTD000479"=91, "QTD000489"=91, "QTD000499"=91,
                           "QTD000504"=91)

eqtl.n <- eqtl.sampsize.list[[opt$eqtlName]]

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
# gwas results files - this contains results for all nhoods for each lead SNP +/- 250kb
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

##############################################################
### read in the association summary statistics for eQTL    ###
##############################################################
# run a system call to tabix for a particular SNP
# split this into the summary statistics for each gene

milo.coloc.list <- list()
locus.dir <- opt$GWASDir
track <- 1
for(x in seq_along(lead.snps)){
  x.snp <- lead.snps[x]
  message("Lead SNP ", x.snp)
  x.chr <- as.numeric(gsub(x.snp, pattern="([0-9]+)_([0-9]+)_([ATCG]+)_([ATCG]+)", replacement="\\1"))
  x.bp <- as.numeric(gsub(x.snp, pattern="([0-9]+)_([0-9]+)_([ATCG]+)_([ATCG]+)", replacement="\\2"))
  
  x.gwas.file <- list.files(locus.dir, pattern=x.snp, full.names=TRUE)
  x.gwas.file <- x.gwas.file[grepl(x.gwas.file, pattern="\\.tsv.gz$")]
  
  x.tab <- system(paste("tabix", x.gwas.file, paste0(x.chr, ":", x.bp-250000, "-", x.bp+250000)), intern=TRUE)
  
  x.locus <- do.call(rbind.data.frame, lapply(x.tab, FUN=function(TX){
    unlist(strsplit(TX, split="\t", fixed=TRUE))
    }))
  
  colnames(x.locus) <- c("SNP", "BP", "CHR", "logFC", "SE", "GeneticVariance", "SpatialFDR", "Nhood", "ALT", "REF")
  x.locus$BP <- as.numeric(x.locus$BP)
  x.locus$logFC <- as.numeric(x.locus$logFC)
  x.locus$CHR <- as.numeric(x.locus$CHR)
  x.locus$SE <- as.numeric(x.locus$SE)
  x.locus$GeneticVariance <- as.numeric(x.locus$GeneticVariance)
  x.locus$SpatialFDR <- as.numeric(x.locus$SpatialFDR)
  x.locus$Nhood <- as.numeric(x.locus$Nhood)
  
  # this needs to be a list for each associated nhood.
  x.nhoods <- unique(x.locus$Nhood[x.locus$SNP %in% x.snp & x.locus$SpatialFDR <= 1e-8])
  x.locus <- x.locus[x.locus$Nhood %in% x.nhoods, ]

  message("Results extracted for ", nrow(x.locus), "SNPs and nhoods")
  
  x.coloc.list <- list()
  for(k in seq_along(x.nhoods)){
    x.df <- x.locus[x.locus$Nhood %in% x.nhoods[k], ]
    x.df <- x.df[!is.na(x.df$logFC), ]

    message("Constructing coloc list for nhood ", x.nhoods[k])
    x.nhood.list <- list("beta"=x.df$logFC, "varbeta"=(x.df$SE)**2,
                         "snp"=x.df$SNP, "position"=x.df$BP,
                         "pvalue"=x.df$SpatialFDR, "type"="quant", "sdY"=rand.ni.sd[x.nhoods[k]])

    q.data <- opt$eqtlName
    message("Querying eQTL study: ", q.data)
    q.qtl.file <- opt$eqtl
    q.query <- paste0("tabix ", q.qtl.file, " ", x.chr, ":", x.bp-250000, "-", x.bp+250000)
    message("eQTL tabix query: ", q.query)
    q.df <- do.call(rbind.data.frame,
                    sapply(system(q.query, intern=TRUE),
                           strsplit, split="\t", fixed=TRUE))

    message("Found ", nrow(q.df), " eQTL results at this locus")
    if(nrow(q.df) > 0){
      # this needs to be a colocalisation analysis
      colnames(q.df) <- colnames(readr::read_tsv(q.qtl.file, n_max = 1, show_col_types=FALSE))
      q.df$SNP <- paste(q.df$chromosome, q.df$position, q.df$ref, q.df$alt, sep="_") # note that ref and alt might not be identical
      q.df$DataSet <- q.data
      q.df$position <- as.numeric(q.df$position)
      q.df$pvalue <- as.numeric(q.df$pvalue)
      q.df$beta <- as.numeric(q.df$beta)
      q.df$chromosome <- as.numeric(q.df$chromosome)
      q.df$maf <- as.numeric(q.df$maf)
      
      ## make coloc data set for eQTL for each gene
      q.genes <- unique(q.df$gene_id)
      for(i in seq_along(q.genes)){
        i.trait <- q.genes[i]
        # remove duplicated SNP entries
        iq.df <- dplyr::distinct(q.df[q.df$gene_id %in% q.genes[i], ], SNP, gene_id, .keep_all=TRUE)
        q.eqtl_coloc.list <- list("beta"=iq.df$beta, "MAF"=iq.df$maf, 
                                  "snp"=iq.df$SNP, "position"=iq.df$position,
                                  "type"="quant", "N"=eqtl.n,
                                  "pvalues"=iq.df$pvalue)
        qk.snps <- intersect(q.eqtl_coloc.list$snp, x.nhood.list$snp)
	message(length(qk.snps), " overlapping SNPs to analyse")
	message("Coloc check output: ", check_dataset(q.eqtl_coloc.list))
        if(is.null(check_dataset(q.eqtl_coloc.list)) & length(qk.snps)){
          # do the coloc here.
          qk.log <- capture.output({
            # coloc.abf is very verbose - do this to silence the printed output
            qk.coloc <- coloc.abf(x.nhood.list, q.eqtl_coloc.list)
            })
            
          qk.df <- do.call(cbind.data.frame, as.list(qk.coloc$summary))
          qk.df$Nhood <- x.nhoods[k]
          qk.df$LeadSNP <- x.snp
          qk.df$eQTLDataset <- q.data
          qk.df$eGene <- i.trait
          milo.coloc.list[[paste0(track)]] <- qk.df
          track <- track + 1
          }
        }
      }
  }
}

eqtl.coloc.df <- do.call(rbind.data.frame, milo.coloc.list)

message("Writing coloc results to: ", opt$output)
write.table(eqtl.coloc.df, file=opt$output,
            sep="\t", quote=FALSE, row.names=FALSE)

message("All done")






