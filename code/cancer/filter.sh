#!/bin/bash

# Filter out unique values from boundaries
FNAME=$1
sort $FNAME.boundaries.bed \
    | uniq \
    | sort -k 1,1 -k2,2n \
> $FNAME.unique.boundaries.bed
