#!/bin/bash

# check if element in an array
containsElement () {
    local e
    x=1
    for e in "${@:2}"; do [[ "$e" == "$1" ]] && x=0; done
    echo $x
}

# make sure the user supplies a resolution
if (( $# != 1 ))
then
    echo "Error! No resolution specified!"
    echo "Usage: script.sh [RESOLUTION]"
    exit 1
fi

# the resolution is argument 1
RES=$1

# array of possible resolutions
RESOLUTIONS=("100kb" "200kb" "400kb" "800kb" "1000kb")

f=$(containsElement $RES "${RESOLUTIONS[@]}")

if (("$f" > 0))
then
    echo "Unsupported resolution."
    exit 1
fi

FNAME="IMR90-$RES.peaks"
FPATH="/home/jniles/thesis/sync/data/domains/boundaries/"
HG19="/home/jniles/data/dna/hg19/hg19.chrom.sizes"
WORKING="working/"

# slop the boundaries in question at 10kb
bedtools slop \
         -b 500000 \
         -i $FPATH$FNAME.bed \
         -g $HG19 \
> $WORKING$FNAME.plusminus.500kb.bed

# make windows for the new file
bedtools makewindows \
         -b $WORKING$FNAME.plusminus.500kb.bed \
         -w 100 \
         -i srcwinnum \
| sort -k1,1 -k2,2n \
| tr "_" "\t" \
> $WORKING$FNAME.plusminus.500kb.100bp.windows.bed

# we now have our base windows, time to map the dna binding protein

BASE=$WORKING$FNAME.plusminus.500kb.100bp.windows.bed
TFNAME="pol2"
TF="/home/jniles/data/IMR90/epi/$TFNAME/GSM935513_hg19_wgEncodeSydhTfbsImr90Pol2IggrabSig.bedg"

# -c 4 -o mean: get the mean of the coverage
# -null 0: if no overlap with bigwig, set to zero
bedtools map \
         -a $BASE \
         -b $TF \
         -c 4 \
         -o mean \
         -null 0 \
> $WORKING$TFNAME.window.coverage.bedg

# sort by the window number
# -t$'\t' to specify that TABS should 
# be used as the delimiter
sort -t$'\t' -k5,5n $WORKING$TFNAME.window.coverage.bedg \
| bedtools groupby \
           -i - \
           -g 5 \
           -c 6 \
           -o sum \
> $TFNAME.$RES.window.counts.txt
