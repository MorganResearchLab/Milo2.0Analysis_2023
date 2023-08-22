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

parser <- add_option(parser, c("-t", "--sumstats"), type="character",
                    help="File containing summary statistics from COVID19 HGI")

parser <- add_option(parser, c("-i", "--window"), type="numeric",
                     default=10000, help="Window size around lead SNP to extract")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output file prefix for DA testing results")

opt <- parse_args(parser)

# define a function for mapping multi-based variants to opposite strand
mapMulti <- function(ALLELE, NTLIST){
  # map each individual NT
  x.alleles <- unlist(strsplit(ALLELE, split="", fixed=TRUE))
  x.flip <- paste0(unlist(sapply(x.alleles, FUN=function(NT, nt.list) nt.list[[NT]],
  	           		 simplify=FALSE, nt.list=NTLIST)))
  return(x.flip)
}


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

message("Running analysis for COVID19 summary statistics")

# make a list that reverse complements a sequence for strand-flipping SNP alleles
nt.list <- list("A"="T", "C"="G", "T"="A", "G"="C")

##############################################################
### Read in the association summary statistics for eQTL    ###
##############################################################
# run a system call to tabix for a particular SNP
# split this into the summary statistics for each gene

milo.coloc.list <- list()
locus.dir <- opt$GWASDir
track <- 1
l.window <- floor(opt$window/2)


for(x in seq_along(lead.snps)){
  x.snp <- lead.snps[x]
  message("Lead SNP ", x.snp)
  x.chr <- as.numeric(gsub(x.snp, pattern="([0-9]+)_([0-9]+)_([ATCG]+)_([ATCG]+)", replacement="\\1"))
  x.bp <- as.numeric(gsub(x.snp, pattern="([0-9]+)_([0-9]+)_([ATCG]+)_([ATCG]+)", replacement="\\2"))
  
  x.gwas.file <- list.files(locus.dir, pattern=x.snp, full.names=TRUE)
  x.gwas.file <- x.gwas.file[grepl(x.gwas.file, pattern="\\.tsv.gz$")]
  
  x.tab <- system(paste("tabix", x.gwas.file, paste0(x.chr, ":", x.bp-l.window, "-", x.bp+l.window)), intern=TRUE)
  
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

  message("Results extracted for ", nrow(x.locus), " SNPs and nhoods")
  
  x.coloc.list <- list()
  for(k in seq_along(x.nhoods)){
    x.df <- x.locus[x.locus$Nhood %in% x.nhoods[k], ]
    x.df <- x.df[!is.na(x.df$logFC), ]

    q.query <- paste0("tabix ", opt$sumstats, " ", x.chr, ":", x.bp-l.window, "-", x.bp+l.window)
    message("GWAS sumstats tabix query: ", q.query)
    q.df <- do.call(rbind.data.frame,
                    sapply(system(q.query, intern=TRUE),
                           strsplit, split="\t", fixed=TRUE))
      
    if(nrow(q.df) > 0){
      # this needs to be a colocalisation analysis
      # what are the column names?

      colnames(q.df) <- colnames(readr::read_tsv(opt$sumstats, n_max = 1, show_col_types=FALSE))
      colnames(q.df) <- gsub(colnames(q.df), pattern="#", replacement="")
      print(head(q.df))
      q.df$SNP <- gsub(q.df$SNP, pattern=":", replacement="_")
      q.df$all_inv_var_meta_beta <- as.numeric(q.df$all_inv_var_meta_beta)
      q.df$all_inv_var_meta_sebeta <- as.numeric(q.df$all_inv_var_meta_sebeta)
      q.df$all_inv_var_meta_cases <- as.numeric(q.df$all_inv_var_meta_cases)
      q.df$all_inv_var_meta_controls <- as.numeric(q.df$all_inv_var_meta_controls)
      q.df$BP <- as.numeric(q.df$POS)
      q.df$N <- q.df$all_inv_var_meta_case + q.df$all_inv_var_meta_controls
      
      q.df <- q.df[!is.na(q.df$all_inv_var_meta_beta), ]
      q.df <- q.df[!is.na(q.df$all_inv_var_meta_sebeta), ]
      q.df <- q.df[!is.na(q.df$N), ]
      q.df <- q.df[!duplicated(q.df$SNP), ]
      
      # sum over both cases and controls
      eqtl.n <- max(na.omit(q.df$N))
      eqtl.n <- eqtl.n[!is.infinite(eqtl.n)]
      print(eqtl.n)
	
      # drop NA values in betas
      q.df <- q.df[!is.na(q.df$beta), ]
      q.df <- q.df[!is.na(q.df$maf), ]

      is.assoc <- any(q.df$all_inv_var_meta_p < 1e-6)
      if(is.assoc){

      # create a merge between the csQTL and GWAS summ stats
      # do this as its easier to resolve the strand swapping problems
      q.x.overlap <- intersect(q.df$SNP, x.locus$SNP)
      if(length(q.x.overlap) >= 5){
	
      message("Constructing coloc list for nhood ", x.nhoods[k])
      x.nhood.list <- list("beta"=x.locus$logFC, "varbeta"=(x.locus$SE)**2,
                           "snp"=x.locus$SNP, "position"=x.locus$BP,
                           "type"="quant", "sdY"=rand.ni.sd[x.nhoods[k]])

      # make coloc data set for eQTL for each gene
      # remove duplicated SNP entries
      if(length(eqtl.n)){
	    q.eqtl_coloc.list <- list("beta"=q.df$all_inv_var_meta_beta,
                                     "varbeta"=(q.df$all_inv_var_meta_sebeta)**2,
                                     "snp"=q.df$SNP,
                                     "position"=q.df$POS,
                                     "type"="cc", "N"=eqtl.n)
      } else{
	    warning("No case numbers found for ", q.data, ". Skipping analysis")
	    break
      }

      # just check the position overlap.
      x.nhood.snps <- gsub(x.nhood.list$snp, pattern="([0-9]+)_([0-9]+)_([ATCG]+)_([ATCG]+)", replacement="\\1_\\2")
      q.qtl.snps <- gsub(q.eqtl_coloc.list$snp, pattern="([0-9]+)_([0-9]+)_([ATCG]+)_([ATCG]+)", replacement="\\1_\\2")
      qk.snps <- intersect(x.nhood.snps, q.qtl.snps)
      q.which.snps <- which(gsub(q.eqtl_coloc.list$snp, pattern="([0-9]+)_([0-9]+)_([ATCG]+)_([ATCG]+)", replacement="\\1_\\2") == qk.snps, arr.ind=TRUE)

      # if this is < 10% of SNPs then check the position overlaps
      # as the discordant could be the ref/alt definitions
      message("There are ", length(qk.snps), " SNPs that overlap for analysis")
	
      message("Coloc check output: ", check_dataset(q.eqtl_coloc.list))
      if(is.null(check_dataset(q.eqtl_coloc.list)) & length(qk.snps) > 2){
            # do the coloc here.
	    # coloc.abf is very verbose - do this to silence the printed output
	    qk.coloc <- coloc.abf(x.nhood.list, q.eqtl_coloc.list)    
	    qk.df <- do.call(cbind.data.frame, as.list(qk.coloc$summary))
	    qk.df$Nhood <- x.nhoods[k]
	    qk.df$LeadSNP <- x.snp
	    qk.df$Pheno <- "SevereCovid19"
	    
	    # check if sum of PPs is > 1
	    pp.sum <- rowSums(qk.df[, c(2:6)])
	    if(any(pp.sum > 1)){
	        warning("Sum of posteriors is greater than 1")
		print(qk.df)
	    }
 	    
	    milo.coloc.list[[paste0(track)]] <- qk.df
	    track <- track + 1
	    }
      }
    }
    }
  }
}

eqtl.coloc.df <- do.call(rbind.data.frame, milo.coloc.list)

ofile <- paste0(opt$output, "_HGICovid19coloc.tsv")
message("Writing coloc results to: ", ofile)
write.table(eqtl.coloc.df, file=ofile,
            sep="\t", quote=FALSE, row.names=FALSE)

message("All done")






