#! /usr/bin/bash

## Extract genotypes from plink file given a set of coordinates
## Input is a file with the following columns:
## 1) Input bed prefix
## 2) Chromosome
## 3) Start
## 4) End
## 5) Output file prefix

FCONT=$(cat $1-$LSB_JOBINDEX".tsv")
BFILE=$(echo $FCONT | cut -d " " -f1)
CHR=$(echo $FCONT | cut -d " " -f2)
START=$(echo $FCONT | cut -d " " -f3)
END=$(echo $FCONT | cut -d " " -f4)
OUT=$(echo $FCONT | cut -d " " -f5)

plink --bfile $BFILE  --recode A-transpose --maf 0.15 --memory 11000  --chr $CHR --from-bp $START --to-bp $END  --out $OUT
