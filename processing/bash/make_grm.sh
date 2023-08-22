#! /usr/bin/bash

# 1) input plink files
# 2) Output file name

gcta64 --mbfile $1 --maf 0.01 --make-grm --threads 10  --out $2


