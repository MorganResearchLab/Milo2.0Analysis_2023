#! /usr/bin/bash

## merge all the GWAS results
## only need one header
# # argv=( "$@" )
# # OUT=${argv[-1]}
# # unset 'argv[-1]'
CHROME=$1
DIR=$2
OUT=$3

echo $CHROME
echo $DIR
echo $OUT
all_res=$(find $DIR"/" -name "*GLMM_results.tsv" | grep "${CHROME}")
#echo ${all_res[1]}

# write the headerline first
printf "logFC\tlogCPM\tSE\ttvalue\tPValue\tGenetic variance\tConverged\tDispersion\tLogliklihood\tNhood\tSpatialFDR\tSNP\tN\n" | gzip > $OUT
#zcat $OUT

# make a regex variable for awk to select the correct
CHREGEX=$(echo "^"$CHROME"_[0-9]+_[ATCG]_[ATCG]$")
echo $CHREGEX

for res in ${all_res[@]};
do
    cat $res | grep -v "logFC" | gzip >> $OUT
done


