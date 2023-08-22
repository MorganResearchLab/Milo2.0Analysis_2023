#! /usr/bin/bash

## 1) input file
## 2) threshold
## 3) output file

zcat $1 | awk -v thresh=$2 '{if($7 <= thresh) {print $0}}'  |  sort -k2n | awk 'BEGIN {printf("SNP\tBP\tCHR\tLFC\tSE\tGeneticVariance\tSpatialFDR\tNhood\n")} {print $0}' | bgzip > $3
tabix -S1  -b 2 -e 2 -s 3  $3
