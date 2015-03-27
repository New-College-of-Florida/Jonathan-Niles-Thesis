#!/usr/bin/env python
# plot a checkerboard pattern for a single chromosome

from __future__ import print_function

import os
import sys
import nutils
import matplotlib.pyplot as plt
import numpy as np
from hiclib.binnedData import binnedData

# Constants
# change this to set the resolution
path = "/home/jniles/data/{0}/ic/{0}-{1}-HindIII-200k.hm"

# useful parameters
genome = "/home/jniles/data/dna/hg19"
cellType = "hESC"
replicates = ["all"]

def getFileName(replicate):
    """
    return the heat map file name for a given replicate and resolution
    """
    assert replicate in replicates
    return path.format(cellType, replicate)

def plotChromosomeHeatMap(replicate, chromosome):
    """
    Plots a heatmap for the IMR90 cell line
    at a given resolution and repliate number
    """
    print("Loading Heatmap Data for chromosome ", chromosome)

    # decriment, to match params
    chromosome -= 1

    inFile = getFileName(replicate)
    key = "%s-%s" % (replicate, "200k")

    Tanay = binnedDataAnalysis(200000, genome=genome)
    Tanay.simpleLoad(inFile, key)

    print("Filtering data...")
    data = filterByChromosome(Tanay, key, chromosome)

    print("Plotting figure ", chromosome)
    plt.figure()

    fig, ax = plt.subplot()

    ax.imshow(np.log(data), cmap=plt.cm.Reds, interpolation='None')

    # Chromosome index needs to be human readable!
    label = chromosome + 1

    plt.title('Chromosome {0} Heatmap (Resolution: {1})'.format(label, resolution))
    plt.xlabel('Chromosome {0}'.format(label))
    plt.ylabel('Chromosome {0}'.format(label))

    savePath = "pngs/heatmap-{0}-{1}-chr{2}.png".format(replicate, resolution, label)
    print("Saving figure to", savePath)
    plt.savefig(savePath, dpi=900)

    # explicitly close the figure to conserve memory
    plt.close()
    return

def filterByChromosome(data, key, chromosome):
    """
    Filters a binnedDataSet by the chromosome number
    provided.
    """
    assert chromosome in data.chromosomeIndex
    # get a boolean mask for bins that correspond to
    # the desired chromosome
    mask = data.chromosomeIndex == chromosome

    # Pull out the raw (binned) matrix
    rawMatrix = data.dataDict[key]

    return squareMask(rawMatrix, mask)

def squareMask(data, mask):
    """
    Returns the square mask of a nxn matrix
    by slicing out the rows and columns.
    """
    return data[mask].T[mask]

# This is will run the entire module
def main(chromosome):
    for i in range(1,20):
        plotChromosomeHeatMap("R1", i)

if __name__ == "__main__":
    main(int(sys.argv[1]))
