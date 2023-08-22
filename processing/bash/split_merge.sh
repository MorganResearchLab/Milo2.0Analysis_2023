#! /usr/bin/bash

INFILE=$1

cat $INFILE | grep -v "logFC" | awk -F '\t|_' '{print $0 >> "/nfs/research/marioni/mdmorgan/Randolph.dir/split_gwas.dir/Chr"$12"_mergedGWAS.tsv"}'
