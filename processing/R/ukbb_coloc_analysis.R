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

parser <- add_option(parser, c("-m", "--manifest"), type="character",
                     help="Manifest file of UKBB phenotypes")

parser <- add_option(parser, c("-w", "--wgetCol"), type="character",
                     help="Column of manifest file with wget path")

parser <- add_option(parser, c("-e","--wgetDir"), type="character",
                     help="Directory to download wget output to")

parser <- add_option(parser, c("-n", "--NCol"), type="character",
                     help="Column of manifest file with sample size")

parser <- add_option(parser, c("-l", "--liftover"), type="character",
                     help="Path to chain file for genome build liftover")

parser <- add_option(parser, c("-c", "--chain"), type="character",
                     help="Path to chain file that reverses genome build liftover")

parser <- add_option(parser, c("-i", "--window"), type="numeric",
                     default=10000, help="Window size around lead SNP to extract")

parser <- add_option(parser, c("-x", "--index"), type="numeric",
                     default=NULL, help="Used for debugging only")

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

message("Reading in manifest file: ", opt$manifest)
manifest.df <- read.delim(opt$manifest, sep="\t", header=TRUE)
all.pheno <- unique(manifest.df$phenocode)
# get LSB_JOBINDEX
if(is.null(opt$index)){
    q.num <- as.numeric(Sys.getenv('LSB_JOBINDEX')) # the phenocode number is the LSB job index
} else{
    q.num <- opt$index
}

q.data <- all.pheno[q.num] # LSB_JOBINDEX starts at 1
q.manifest <- dplyr::distinct(manifest.df[manifest.df$phenocode %in% q.data, ], phenocode, coding, .keep_all=TRUE)
message("Running analysis for ", q.data)

message("Reading in liftOver chain file: ", opt$liftover)
chain.over <- import.chain(opt$liftover)

message("Reading reverse liftOver chain file: ", opt$chain)
chain.reverse <- import.chain(opt$chain)
run.list <- list()

for(i in seq_len(nrow(q.manifest))){
    q.coding <- q.manifest[i, ]$coding
    message("Downloading GWAS results for: ", q.data, " ", q.coding)
    q.which.col <- which(colnames(manifest.df) == opt$wgetCol)
    # wget -q switches on quiet setting. -P redirects download to specific directory
    # check if files already exist

    # the file name the the URL are different
    wget.file <- gsub(as.character(q.manifest[q.manifest$coding %in% q.coding, q.which.col]), pattern="([[:ascii:]]*)/(\\S+\\.tsv\\.bgz)", replacement="\\2", perl=TRUE)
    message("Summ stat file: ", wget.file)
    if(file.exists(paste0(opt$wgetDir, "/", wget.file))){
        message(wget.file, " already exists - skipping download")        
    } else{
        q.wget.url <- gsub(as.character(q.manifest[q.manifest$coding %in% q.coding, q.which.col]), pattern="wget ", replacement="", perl=TRUE)
	message("Downloading sumstats from ", q.wget.url)
	q.wget.down <- system(paste0("wget -q -P ", opt$wgetDir, " ", q.wget.url), intern=TRUE)
    }

    wget.tabx <- paste0(wget.file, ".tbi")    
    message("Tabix file: ", wget.tabx)
    if(file.exists(paste0(opt$wgetDir, "/", wget.tabx))){
        message(wget.tabx, " already exists - skipping download")
    } else{
        q.wget.taburl <- gsub(as.character(q.manifest[q.manifest$coding %in% q.coding, q.which.col+1]), pattern="wget ", replacement="", perl=TRUE)
	message("Downloading tabix from ", q.wget.taburl)
        q.wgettbi.down <- system(paste0("wget -q  -P", opt$wgetDir, " ", q.wget.taburl), intern=TRUE)
    }

    # check files for columns required for colocalisation
    req.cols <- c("neglog10_pval_meta", "beta_meta", "se_meta")
    q.qtl.file <- paste0(opt$wgetDir, "/", wget.file)
    q.header <- colnames(readr::read_tsv(q.qtl.file, n_max = 1, show_col_types=FALSE))
    
    if(all(req.cols %in% q.header)){
         q.check <- TRUE
    } else{
        q.check <- FALSE
    }
    run.list[[q.qtl.file]] <- q.check
}

# fail if no sub-phenotype has proper columns
sum.missing <- sum(!unlist(run.list))
if(sum.missing == length(run.list)){
    stop("No sub-phenotypes contain sufficient summary statistics for colocalisation analysis - exiting")
}

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
  # need to liftover the lead SNP coordinates to get the correct GWAS sum states.
  x.lift <- as.data.frame(liftOver(GRanges(seqnames=paste0("chr", x.chr), ranges=IRanges(start=x.bp, end=x.bp)), chain.reverse))
  x.lift.bp <- x.lift$start
  
  x.coloc.list <- list()
  for(k in seq_along(x.nhoods)){
    x.df <- x.locus[x.locus$Nhood %in% x.nhoods[k], ]
    x.df <- x.df[!is.na(x.df$logFC), ]

    for(i in seq_len(nrow(q.manifest))){
      q.coding <- q.manifest[i, ]$coding
      q.file <- q.manifest[, ]

      q.qtl.pattern <- gsub(as.character(q.manifest[q.manifest$coding %in% q.coding, q.which.col]), pattern="([[:ascii:]]*)/(\\S+\\.tsv\\.bgz)", replacement="\\2", perl=TRUE)
      q.qtl.file <- list.files(opt$wgetDir, pattern=q.qtl.pattern, full.names=TRUE)
      q.qtl.file <- unique(q.qtl.file[grepl(q.qtl.file, pattern="bgz$")])
      q.check <- run.list[[q.qtl.file]]
      
      if(isFALSE(q.check)){
          message("Insufficient summary statistic information to perform colocalisation")
	  break
      }
    
      q.query <- paste0("tabix ", q.qtl.file, " ", x.chr, ":", x.lift.bp-l.window, "-", x.lift.bp+l.window)
      message("GWAS sumstats tabix query: ", q.query)
      q.df <- do.call(rbind.data.frame,
                      sapply(system(q.query, intern=TRUE),
                             strsplit, split="\t", fixed=TRUE))
      
      if(nrow(q.df) > 0){
        # this needs to be a colocalisation analysis
        # what are the column names?

        colnames(q.df) <- colnames(readr::read_tsv(q.qtl.file, n_max = 1, show_col_types=FALSE))
	q.df$SNP <- paste(q.df$chr, q.df$pos, q.df$ref, q.df$alt, sep="_") # note that ref and alt might not be identical
	q.df$DataSet <- paste(q.data, q.coding, sep="_")
	q.df$BP <- as.numeric(q.df$pos)
	q.df$neglog10_pval_meta <- as.numeric(q.df$neglog10_pval_meta)
	q.df$beta_meta <- as.numeric(q.df$beta_meta)
	q.df$se_meta <- as.numeric(q.df$se_meta)
	q.df$CHR <- as.numeric(q.df$chr)
	q.type <- "quant"
	
	if(any(colnames(q.df) %in% c("af_meta"))){
            q.df$maf <- as.numeric(q.df$af_meta)
        } else{
            q.df$maf <- as.numeric(q.df$af_controls_meta)
	    q.type <- "cc"
        }

	eqtl.n <- unique(na.omit(q.manifest[q.manifest$coding %in% q.coding, opt$NCol]))
	eqtl.n <- eqtl.n[!is.infinite(eqtl.n)]
	print(eqtl.n)
	
	# drop NA values in betas
	q.df <- q.df[!is.na(q.df$beta_meta), ]
	q.df <- q.df[!is.na(q.df$maf), ]


	# check if there is _any_ evidence of an trait association at the locus
	# if not then no point in doing coloc analysis
	# I'll be liberal and 1e-6 for this
	is.assoc <- any(q.df$neglog10_pval_meta >= 6)

	if(is.assoc){

 	# the panUKBB positions are GRCh37 while my GWAS results are GRCh38 - go figure!
	# rtracklayer implements liftOver to convert coordinates
	# need to keep the original SNP names for comparison
	q.granges <- GRanges(seqnames=paste0("chr", q.df$CHR), ranges=IRanges(start=q.df$BP, end=q.df$BP), SNP=q.df$SNP)
	q.lift <- as.data.frame(liftOver(q.granges, chain.over))
	q.lift$CHR <- gsub(q.lift$seqnames, pattern="chr", replacement="")

	q.diff <- length(setdiff(q.df$SNP, q.lift$SNP))
	message(q.diff, " SNP positions dropped by liftOver")
	q.df <- merge(q.df, q.lift[, c("CHR", "start", "SNP")], by=c("SNP", "CHR"))
	q.df$CHR <- as.numeric(q.df$CHR)
	q.df$BP <- as.numeric(q.df$start)
	q.df$SNP <- paste(q.df$CHR, q.df$BP, q.df$ref, q.df$alt, sep="_") # note that ref and alt might not be identical

	# create a merge between the csQTL and GWAS summ stats
	# do this as its easier to resolve the strand swapping problems
	q.x.merge <- merge(x.df, q.df, by=c("CHR", "BP"))
	q.x.merge$Match.Class <- "Stick"
	q.x.merge$Match.Class[(q.x.merge$REF == q.x.merge$alt) & (q.x.merge$ref == q.x.merge$ALT)] <- "Swap"
	q.x.merge$Match.Class[(q.x.merge$REF != q.x.merge$ref) & (q.x.merge$REF != q.x.merge$alt)] <- "Strand.Flip"
	q.x.merge$beta_meta[q.x.merge$Match.Class == "Swap"] <- 1/q.x.merge$beta_meta[q.x.merge$Match.Class == "Swap"]

        q.x.merge$Flip.Alt <- unlist(sapply(seq_len(nrow(q.x.merge)), FUN=function(RX){
	        rx.class <- q.x.merge[RX, ]$Match.Class
		if(rx.class %in% "Strand.Flip"){
	          rx.a <- mapMulti(q.x.merge[RX, ]$alt, nt.list)
	        } else if(rx.class %in% "Swap"){
	          rx.a <- q.x.merge[RX, ]$ref
	        } else {
	          rx.a <- q.x.merge[RX, ]$alt
	        }

		if(length(rx.a) > 1){
		    return(rx.a[1])
		} else{
		    return(rx.a)
		}
	    }, simplify=FALSE))

        q.x.merge$Flip.Ref <- unlist(sapply(seq_len(nrow(q.x.merge)), FUN=function(RX){
	        rx.class <- q.x.merge[RX, ]$Match.Class
	        if(rx.class %in% "Strand.Flip"){
	          rx.a <- mapMulti(q.x.merge[RX, ]$ref, nt.list)
	        } else if(rx.class %in% "Swap"){
	          rx.a <- q.x.merge[RX, ]$alt
 	        } else{
	          rx.a <- q.x.merge[RX, ]$ref
	        }
		
		if(length(rx.a) > 1){
		    return(rx.a[1])
		} else{
		    return(rx.a)
		}

	    }, simplify=FALSE))
	
        q.x.merge$SNP <- paste(q.x.merge$CHR, q.x.merge$BP, q.x.merge$Flip.Ref, q.x.merge$Flip.Alt, sep="_")
	# remove inf and NAs
	q.x.merge <- q.x.merge[!is.infinite(q.x.merge$logFC), ]
	q.x.merge <- q.x.merge[!is.na(q.x.merge$logFC), ]
	q.x.merge <- q.x.merge[!is.infinite(q.x.merge$beta_meta), ]
	q.x.merge <- q.x.merge[!is.na(q.x.merge$beta_meta), ]
	q.x.merge <- q.x.merge[!is.infinite(q.x.merge$se_meta), ]
	q.x.merge <- q.x.merge[!is.na(q.x.merge$se_meta), ]
	q.x.merge <- q.x.merge[!is.infinite(q.x.merge$SpatialFDR), ]
	q.x.merge[!is.na(q.x.merge$SpatialFDR), ]$SpatialFDR <- 1

	q.x.merge <- dplyr::distinct(q.x.merge, SNP, .keep_all=TRUE)
	
        message("Constructing coloc list for nhood ", x.nhoods[k])
        x.nhood.list <- list("beta"=q.x.merge$logFC, "varbeta"=(q.x.merge$SE)**2,
                             "snp"=q.x.merge$SNP, "position"=q.x.merge$BP,
                             "type"="quant", "sdY"=rand.ni.sd[x.nhoods[k]])

      	## make coloc data set for eQTL for each gene
	# remove duplicated SNP entries
	if(length(eqtl.n)){
	    q.eqtl_coloc.list <- list("beta"=q.x.merge$beta_meta, "MAF"=q.x.merge$maf,
        		          "varbeta"=(q.x.merge$se_meta)**2,
                        	  "snp"=q.x.merge$SNP, "position"=q.x.merge$BP,
                                  "type"=q.type, "N"=eqtl.n)
	} else{
	    warning("No case numbers found for ", q.data, " ", q.coding, ". Skipping analysis")
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
	if(is.null(check_dataset(q.eqtl_coloc.list)) & length(qk.snps) > 5){
            # do the coloc here.
	    # coloc.abf is very verbose - do this to silence the printed output
	    qk.coloc <- coloc.abf(x.nhood.list, q.eqtl_coloc.list)           
	    qk.df <- do.call(cbind.data.frame, as.list(qk.coloc$summary))
	    qk.df$Nhood <- x.nhoods[k]
	    qk.df$LeadSNP <- x.snp
	    qk.df$UKBBPheno <- paste(q.data, q.coding, sep="-")
	    qk.df$PhenoDescription <- as.character(manifest.df[manifest.df$phenocode %in% q.data, "description"])
	    # check if sum of PPs is > 1
	    pp.sum <- rowSums(qk.df[, c(2:6)])
	    if(any(pp.sum > 1)){
	        warning("Sum of posteriors is greater than 1")
		print(qk.df)
	    }
 	    
	    milo.coloc.list[[paste0(track)]] <- qk.df
	    track <- track + 1
	}
     } else{
       warning("No association found with ", q.data, " at locus with ", x.snp)
     }
     }
    }
  }
}

eqtl.coloc.df <- do.call(rbind.data.frame, milo.coloc.list)

ofile <- paste0(opt$output, q.data, "_UKBBcoloc.tsv")
message("Writing coloc results to: ", ofile)
write.table(eqtl.coloc.df, file=ofile,
            sep="\t", quote=FALSE, row.names=FALSE)

message("All done")






