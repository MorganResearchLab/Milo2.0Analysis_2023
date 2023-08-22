#! /usr/bin/env Rscript

library(optparse)
library(reshape2)
library(readr)
library(imager)
library(data.table)
library(Cairo)

parser <- OptionParser()

parser <- add_option(parser, c("-i", "--indir"), type="character",
                     help="A directory containing Milo-GWAS results with tabix indices")

parser <- add_option(parser, c("-r", "--regex"), type="character",
                     help="Regex to select chromsome-specific files")

parser <- add_option(parser, c("-s", "--suffix"), type="character",
                     help="File suffix for GWAS tabix index results files")

parser <- add_option(parser, c("-n", "--nchunks"), type="numeric", default=100,
                     help="Number of chunks to break each chromosome into")

parser <- add_option(parser, c("-p", "--plink"), type="character",
                     help="Path to plink bim files with SNP info")

parser <- add_option(parser, c("-t", "--threshold"), type="character",
                     help="SpatialFDR threshold for defining index SNPs")

parser <- add_option(parser, c("-w", "--window"), type="numeric", default=250000,
                     help="Window size around each index SNP to use for clumping")

parser <- add_option(parser, c("-o", "--outdir"), type="character",
                     help="Output directory for each locus results to be written to")

parser <- add_option(parser, c("-l", "--summary"), type="character",
                     help="Directory containing tabix indexed sum stat files")

opt <- parse_args(parser)

# this makes a massive manhattan plot for all variants across all traits (mean or variability)
#test_results <- test_results[grepl(test_results, pattern="mlma$")]
chromes <- paste0("Chr", 1:22)

# first loop over one trait for all chromosomes and count the number of SNPs on each contig
bim.files <- list.files(opt$plink, pattern="bim$", full.names=TRUE)
chrome_size_list <- list()

message("Measuring chromosome lengths")
for(j in seq_along(chromes)){
    check.chr <- bim.files[grepl(bim.files, pattern=chromes[j])]
    j.bim <- read_delim(check.chr, col_names=FALSE, delim="\t", show_col_types=FALSE)
    bp <- as.vector(j.bim[, 4, drop=FALSE])
    n.hits <- max(bp) - min(bp)
    chrome_size_list[[chromes[j]]] <- n.hits
    message("Chromosome ", chromes[j], " has ", n.hits, " SNP positions")
}

total.snps <- sum(unlist(chrome_size_list))
message(total.snps, " base pairs to plot in total")
plot.size <- c()

chrom.cols <- rep(c("#ca0020", "#404040"), 11)
names(chrom.cols) <- paste0("Chr", c(1:22))

## we can use tabix to subset the results
## there is a tabix file per chromosome
## is this going to be very heavy on the I/O? Maybe might want to limit the number
## of chunks for each chromosome.
test_results.list <- list()
test_results <- list.files(opt$indir, pattern=opt$regex, full.names=TRUE)

for(i in seq_along(chromes)){
    message("Collecting GWAS results for chromsome ", chromes[i])
    i.int.chr <- gsub(chromes[i], pattern="Chr", replacement="") # need the integer version for tabix
    
    chrome.results.list <- list()
    i.size <- chrome_size_list[[chromes[i]]] # this is 1-based indexing
    i.chunksize <- as.integer(floor(i.size/opt$nchunks))
    # there are many files for each chromosome
    chrome_test_results <- test_results[grepl(test_results, pattern=paste0(chromes[i], opt$suffix))]

    message("Found ", length(chrome_test_results), " results files to parse")
    if(length(chrome_test_results) > 0){
        # break the whole chromosome into chunks - 100?
	i.min <- 1
	i.max <- as.integer(floor(i.size/opt$nchunks))

	for(x in seq_len(opt$nchunks)){
	    x.command <- paste0("tabix ", chrome_test_results, " ", i.int.chr, ":", i.min, "-", i.max)
	    message("Tabix command: ", x.command)
	    x.tabix <- system(x.command, intern=TRUE)

	    if(length(x.tabix) > 0){
  	        x.df <- do.call(rbind.data.frame, sapply(x.tabix, strsplit, split="\t", fixed=TRUE))
		message(nrow(x.df), " GWAS results extracted")
		colnames(x.df) <- c("SNP", "BP", "CHR", "LFC", "SpatialFDR", "Nhood")
		x.df$BP <- as.numeric(x.df$BP)
		x.df$SpatialFDR <- as.numeric(x.df$SpatialFDR)
		chrome.results.list[[paste0(x)]] <- x.df
	    }
    	    i.min <- i.max + 1
	    i.max <- i.max + i.chunksize

	    if(i.max > i.size){
	        i.max <- i.size
	    }
        }
    }

    chrome.df <- do.call(rbind.data.frame, chrome.results.list)

    message("Flushing environment")
    sink(file="/dev/null")
    rm(list=c("chrome.results.list"))
    gc()
    sink(file=NULL)

    if(nrow(chrome.df) > 0){
        # number of unique SNPs
	chr.snps <- unique(chrome.df$SNP)
	n.snps <- length(chr.snps)
	message("GWAS hits for ", n.snps, " SNPs found")
	message("Looping over SNPs to extract surrounding loci")
	snp.locus.list <- list()

	# there needs to be 2 types of output - the same SNP across all nhoods
	# and the same nhood across SNPs - I'll do the latter first
	# as this is mechanistically easier for locuszoom plotting
    	summ.file <- list.files(opt$summary, pattern=paste0(chromes[i], "_sumStats\\.tsv\\.bgz"), full.names=TRUE)
	print(paste0(chromes[i], "_sumStats.tsv.bgz"))
	print(opt$summary)
	print(chromes[i])
	print(summ.file)

	# do this in order of most signifiant SNPs first
	# this way I can check if they are already in an extract locus
	q.order <- unique(chrome.df[order(chrome.df$SpatialFDR, decreasing=FALSE), ]$SNP)

	seen.snps <- c()
	
	for(q in seq_along(q.order)){
	    q.snp <- q.order[q]
	    message("SNP: ", q.snp)
	    print(q.snp %in% seen.snps)
	    if(q.snp %in% seen.snps){
	        message("Skipping ", q.snp, " - already included")
	    } else{
		q.nhoods <- unique(chrome.df[chrome.df$SNP %in% q.snp, ]$Nhood)
		q.bp <- unique(chrome.df[chrome.df$SNP %in% q.snp, ]$BP)
		q.start <- q.bp - floor(opt$window/2)
		q.end <- q.bp + floor(opt$window/2)
   		if(q.start < 1){
	            q.start <- 1
	        }

	        if(q.end > i.size){
	            q.end <- i.size
	        }

	        # get tabix indexed results
		print(summ.file)
		tab.command <- paste0("tabix ", summ.file, " ", i.int.chr, ":", q.start, "-", q.end)
		sum.tabix <- system(tab.command, intern=TRUE)
		message("Tabix command: ", tab.command)

		print(head(sum.tabix))
	        if(length(sum.tabix) > 0){
	      	    message("Extracting surrounding SNPs: ", paste0(i.int.chr, ":", q.start, "-", q.end, "\n"))
	            sum.df <- do.call(rbind.data.frame, sapply(sum.tabix, strsplit, split="\t", fixed=TRUE))
	            colnames(sum.df) <- c("SNP", "BP", "CHR", "LFC", "SE", "GeneticVariance", "SpatialFDR", "Nhood")
		    # add ref and alt for locuszoom plotting
		    sum.df$REF <- gsub(sum.df$SNP, pattern="([0-9]+)_([0-9]+)_([ATGC]+)_([ATCG]+)", replacement="\\3")
		    sum.df$ALT <- gsub(sum.df$SNP, pattern="([0-9]+)_([0-9]+)_([ATGC]+)_([ATCG]+)", replacement="\\4")
		    
		    q.file <- paste0(opt$outdir, paste0(q.snp, "_locus.tsv"))
		    message("Writing results to ", q.file)
		    write.table(sum.df, file=q.file, sep="\t", quote=FALSE, row.names=FALSE)
		    seen.snps <- unique(c(seen.snps, unique(sum.df$SNP)))
	        }	        
	    }
        }
	
    } else{
        message("No results - skipping ", chromes[i])
    }
}





message("All done")
