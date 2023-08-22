#! /usr/bin/env Rscript

## Perform DA testing with Milo & edgeR

library(BiocParallel)
library(ggplot2)
library(cowplot)
library(scattermore)
library(reshape2)
library(readr)
library(optparse)

parser <- OptionParser()

parser <- add_option(parser, c("-i", "--indir"), type="character",
                     help="A directory containing Milo-GWAS results")

parser <- add_option(parser, c("-r", "--regex"), type="character",
                     help="Regex to select chromsome-specific files")

parser <- add_option(parser, c("-o", "--outfile"), type="character",
                     help="Output file for Manhattan plot")

opt <- parse_args(parser)

# input is a list
print(opt$indir)
print(opt$regex)
in.list <- list.files(opt$indir, pattern=opt$regex, full.names=TRUE)
print(length(in.list))
message("Found ", length(in.list), " GWAS files")
snp.list <- list()

for(x in seq_along(in.list)){
    message("Processing file ", x)
    snp.chunks <- read_delim(in.list[[x]], delim="\t", col_names=TRUE, progress=FALSE, show_col_types=FALSE)
    snp <- snp.chunks$SNP
    chrome <- gsub(snp, pattern="([0-9]+)_([0-9]+)_([ATCG])_([ATCG])$", replacement="\\1")
    bp <- gsub(snp, pattern="([0-9]+)_([0-9]+)_([ATCG])_([ATCG])$", replacement="\\2")
    snp.chunks$SpatialFDR[snp.chunks$SpatialFDR == 0.0] <- 1e-308
    snp.chunks$SpatialFDR[is.na(snp.chunks$SpatialFDR)] <- 1.0

    snp.list[[paste0(x)]] <- data.frame("SNP"=snp, "BP"=bp, "CHR"=chrome, "LFC"=as.numeric(snp.chunks$logFC),
    			     	        "SE"=as.numeric(snp.chunks$SE),
    			     	        "Genetic.Variance"=as.numeric(snp.chunks$`Genetic variance`), 
    			                "SpatialFDR"=as.numeric(snp.chunks$SpatialFDR), "Nhood"=snp.chunks$Nhood)
    message("Found ", sum(snp.chunks$SpatialFDR < 1e-8), " SNPs with FDR <= 1e-8")
}

message("Concatenating data")
manhattan.df <- do.call(rbind.data.frame, snp.list)

man.file <- gsub(opt$outfile, pattern="manhattanplot\\.png", replacement="sumStats.tsv.gz")
message("Saving Manhattan plot data file to ", man.file)
write.table(manhattan.df, file=gzfile(man.file, "w"), sep="\t",
            row.names=FALSE, quote=FALSE)

message("Flushing the environment")
sink(file="/dev/null")
rm(list=c("snp.list"))
gc()
sink(file=NULL)

#manhattan.df
n.snps <- length(unique(manhattan.df$SNP))

message("Making Manhattan plot with ", n.snps, " SNPs")
man.plot <- ggplot(manhattan.df, aes(x=BP, y=-log10(SpatialFDR))) +
    geom_scattermore(pointsize=0.5, colour="grey80") +
    geom_hline(yintercept=8, lty=2, colour='orange') +
    theme_cowplot() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    labs(x="Position (BP)", y=expression(paste("-log"[10], " Spatial FDR"))) +
    NULL

message("Saving plot to ", opt$outfile)
ggsave(man.plot, filename=opt$outfile, bg='white', dpi=300,
       height=5, width=8)
       
message("All done")
