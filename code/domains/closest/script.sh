#!/bin/sh

RES=10kb
BEDPATH=../bed

# compare across all replicates
bedtools closest \
    -a $BEDPATH/IMR90-R1-$RES.sorted.bed \
    -b $BEDPATH/IMR90-R2-$RES.sorted.bed \
       $BEDPATH/IMR90-R3-$RES.sorted.bed \
       $BEDPATH/IMR90-R4-$RES.sorted.bed \
       $BEDPATH/IMR90-R5-$RES.sorted.bed \
       $BEDPATH/IMR90-R6-$RES.sorted.bed \
    -mdb each \
    > IMR90-$RES.all.closest.bed

# compare to self
bedtools closest \
    -a $BEDPATH/IMR90-R1-$RES.sorted.bed \
    -b $BEDPATH/IMR90-R1-$RES.sorted.bed \
    -io \
    > IMR90-$RES.self.closest.bed

# Everyday I'm shuffling ...
bedtools shuffle \
    -i $BEDPATH/IMR90-R1-$RES.sorted.bed \
    -g /home/jniles/data/dna/hg19/hg19.chrom.sizes \
    > IMR90-R1-$RES.shuffled.bed

# sort
sort -k1,1 -k2,2n IMR90-R1-$RES.shuffled.bed > IMR90-R1-$RES.shuffled.sorted.bed

# compare to shuffled
bedtools closest \
    -a $BEDPATH/IMR90-R1-$RES.sorted.bed \
    -b IMR90-R1-$RES.shuffled.sorted.bed \
    -io \
    > IMR90-$RES.shuffle.closest.bed
