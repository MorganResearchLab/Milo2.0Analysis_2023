#! /usr/bin/bash

## This needs to be a job array

## 1) Chromosome ID
## 2) PLINK file prefix
## 3) GRM file prefix
## 4) Milo object
## 5) Output file name prefix

## Loop over the SNPs on the chromsome - do this in chunks?
BASE_DIR="/nfs/research/marioni/mdmorgan/Randolph.dir"
SRC=$(echo $BASE_DIR"/src")
module load r/4.2.2

## do a quick count of the number of SNPs to get ~100 SNPs per block.
A=$(cat $2.map | wc -l)
SNP=$(expr $A / 10000) # number of SNP chunks - this is how many array jobs are needed.

FIXED="SOC_genetic_ancestry"
SOLVER="HE-NNLS"
SAMPID="SOC_indiv_ID" # sample ID column - this should match the column names of the nhoods count matrix
RAND="NA" # this should be NA as the GRM is the only random effect
SOLVER="HE-NNLS" # GLMM Variance solver
REDIM="Cosine.MNN"
MAF=0.15

## LSB_JOBINDEX is accessed internally in the script - I couldn't figure out how to pass it in this command without it either defaulting to 0 or being passed as a literal string
bsub -q standard -M 40000 -R "rusage[mem=38000]" -T 120 -n 4 -J $1-Milo_GWAS_[1-$SNP] -e $BASE_DIR/logs/$1-Milo_GWAS_%I.err -o $BASE_DIR/logs/$1-Milo_GWAS_%I.out Rscript $SRC/milo_gwas.R --Milo $4 --fixed $FIXED --sampleID $SAMPID --reddim $REDIM  --GRM $3 --plink $2 --solver $SOLVER --maf $MAF --chunk 10000 --glm --output $5

#echo $JOB
#eval $JOB
