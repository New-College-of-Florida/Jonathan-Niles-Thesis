#!/bin/sh

BEDFILE=$1
OUTFILE=$(echo $BEDFILE | sed 's/\.bed$/.peaks.bed/')

bedtools closest \
    -a ../$BEDFILE \
    -b ../$BEDFILE \
    -io -d | awk -F '\t' \
    'BEGIN {OFS="\t"} { if ($15 == 1)  print }' \
    > $OUTFILE
