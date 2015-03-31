#!/bin/bash

FNAME="hESC-200kb.boundaries"
FPATH="/home/jniles/thesis/sync/data/domains/boundaries/"
HG19="/home/jniles/data/dna/hg19/hg19.chrom.sizes"

# slop the boundaries in question at 10kb
bedtools slop \
         -b 1000 \
         -i $FPATH$FNAME.bed \
         -g $HG19 \
> $FNAME.plusminus.1kb.bed

# make windows for the new file
bedtools makewindows \
         -b $FNAME.plusminus.1kb.bed \
         -w 5 \
         -i srcwinnum \
| sort -k1,1 -k2,2n \
| tr "_" "\t" \
> $FNAME.plusminus.1kb.5bp.windows.bed

# we now have our base windows, time to map the dna binding protein

BASE=$FNAME.plusminus.1kb.5bp.windows.bed
TF="/home/jniles/data/IMR90/epi/ctcf/GSM935404_hg19_wgEncodeSydhTfbsImr90CtcfbIggrabSig.bedg"
TFNAME="ctcf"

# -c 4 -o mean: get the mean of the coverage
# -null 0: if no overlap with bigwig, set to zero
bedtools map \
         -a $BASE \
         -b $TF \
         -c 4 \
         -o mean \
         -null 0 \
> $TFNAME.window.coverage.bedg

# sort by the window number
# -t$'\t' to specify that TABS should 
# be used as the delimiter
sort -t$'\t' -k5,5n $TFNAME.window.coverage.bedg \
| bedtools groupby \
           -i - \
           -g 5 \
           -c 6 \
           -o sum \
> $TFNAME.window.counts.txt

