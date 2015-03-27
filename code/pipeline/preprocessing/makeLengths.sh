#!/bin/bash
# Author : jniles

# This script makes a folder which keeps lengths of all SRA files,
# for each replicate.
# Note that .sra file names arunique

# folder where the .sra file lengths with be stored
if [ ! -d "lengths" ]; then
    mkdir lengths
fi

celltype=$1

# Directory where the replicate samples are stored
DIR="/home/jniles/data/$celltype/sra/"

# For each, use fastq-dump and wc to find the length of each file
for i in `ls | grep r`; do
  for sra in `ls $DIR$i --hide=EXPERIMENT`; do
    fastq-dump -Z $DIR$i/$sra | head -n 2 | tail -n 1 | wc -c > lengths/$sra;
  done
done
