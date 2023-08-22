#! /usr/bin/env Rscript

## convert genotypes in to plink .ped format
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-g", "--genotype"), type="character",
       	             help="A table of variants X samples with additive genotypes")

parser <- add_option(parser, c("-m", "--snps"), type="character",
                     help="SNP position file - 1 SNP per row")

parser <- add_option(parser, c("-c", "--chrome"), type="character",
                     help="Chromosome to subset to")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output file prefix for DA testing results")

opt <- parse_args(parser)

message("Reading in genotypes: ", opt$genotype)
geno <- read.table(opt$genotype, sep="\t", row.names=1, stringsAsFactors=FALSE)

message("Reading in SNP info: ", opt$snps)
snps <- read.table(opt$snps, row.names=1, stringsAsFactors=FALSE)

print(head(snps))

# make the first 6 columns - the same as the .fam file
fam.mat <- data.frame("FAMID"=colnames(geno), "IID"=colnames(geno), "FID"=rep('0', ncol(geno)), "MID"=rep('0', ncol(geno)), "SEX"=rep('1', ncol(geno)),
	              "PHENO"=rep('-9', ncol(geno)))
# and the .map file
map.mat <- data.frame("CHR"=snps$chr, "SNP"=snps$snp, "BP"=snps$pos)

print(head(map.mat))

message("Subsetting to chromosome Chr", opt$chrome)
keep.snps <- map.mat$SNP[map.mat$CHR %in% paste0("chr", opt$chrome)]

message("Keeping ", length(keep.snps), " SNPs")
geno <- geno[rownames(geno) %in% keep.snps, ]

## first we need to split the data into 2 calls per genotype
ped.matrix <- matrix(NA, nrow=ncol(geno), ncol=nrow(geno))

# emit .ped files for each chromosome
for(x in seq_along(rownames(geno))){
    x.snp <- rownames(geno)[x]
    x.geno <- as.vector(geno[x.snp, ]) # vector of 0,1,2 (-9) genotypes
    x.a1 <- gsub(x.snp, pattern="(\\S+)_([0-9]+)_([ATCG]+)_([ATCG]+)", replacement="\\3")
    x.a2 <- gsub(x.snp, pattern="(\\S+)_([0-9]+)_([ATCG]+)_([ATCG]+)", replacement="\\4")
    x.chr <- gsub(x.snp, pattern="(\\S+)_([0-9]+)_([ATCG]+)_([ATCG]+)", replacement="\\1")

    x.left <- character(length(x.geno))
    x.left[x.geno == "-9"] <- "0"
    x.left[x.geno == "0"] <- x.a1
    x.left[x.geno == "1"] <- x.a1
    x.left[x.geno == "2"] <- x.a2

    x.right <- character(length(x.geno))
    x.right[x.geno == "-9"] <- "0"
    x.right[x.geno == "0"] <- x.a1
    x.right[x.geno == "1"] <- x.a2
    x.right[x.geno == "2"] <- x.a2

    ped.matrix[ , x] <- paste0(x.left, x.right)
    if(x %% 10000 == 0){
       message("Processing SNP ", x, " of ", nrow(geno))
     }
}

# make the first 6 columns
fam.mat <- data.frame("FAMID"=colnames(geno), "IID"=colnames(geno), "FID"=rep('0', ncol(geno)), "MID"=rep('0', ncol(geno)), "SEX"=rep('1', ncol(geno)),
	              "PHENO"=rep('-9', ncol(geno)))
fam.out <- paste0(opt$output, "-Chr", opt$chrome, ".fam")
message("Writing .fam file to: ", fam.out)
write.table(fam.mat, file=fam.out, col.names=FALSE, sep=" ", quote=FALSE, row.names=FALSE)

map.mat <- data.frame("CHR"=snps$chr, "SNP"=snps$snp, "BP"=snps$pos)
map.out <- paste0(opt$output, "-Chr", opt$chrome, ".map")
message("Writing .map file to: ", map.out)
write.table(map.mat[map.mat$SNP %in% rownames(geno),], file=map.out, col.names=FALSE, sep=" ", quote=FALSE, row.names=FALSE)

ped.out <- paste0(opt$output, "-Chr", opt$chrome,  ".ped")
ped.mat <- do.call(cbind.data.frame, list(fam.mat, ped.matrix))

message("Writing full .ped file to: ", ped.out)
write.table(ped.mat, file=ped.out, col.names=FALSE, row.names=FALSE, quote=FALSE, sep=" ")

message("All done")