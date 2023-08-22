#! /usr/bin/bash

## 1) Input gzipped tsv file
## 2) Position column number
## 3) Output file name
## 4) tmp file location

zcat $1 | grep -v "SNP" | sort -k$2n -T $4 | awk 'BEGIN{printf("CHR_A\tBP_A\tSNP_A\tCHR_B\tBP_B\tSNP_B\tR2\n")}{print $0}' |  bgzip > $3;
tabix -f -S1 -b 2 -e 2 -s 1 $3
