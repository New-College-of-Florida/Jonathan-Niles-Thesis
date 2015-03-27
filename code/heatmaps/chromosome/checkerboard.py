#!/usr/bin/env python
# plot a checkerboard pattern for a single chromosome

from __future__ import print_function

import os
import nutils
import numpy as np
import nutils as nu
import matplotlib.pyplot as plt
from hiclib.binnedData import binnedData
from hiclib.highResBinnedData import HiResHiC

replicates = {
    "IMR90" : ["R1", "R2", "R3", "R4", "R5", "R6"],
    "hESC" : ["R1", "R2"]
}

# Constants
path = "/home/jniles/data/{0}/ic/"
# useful parameters
genome = "/home/jniles/data/dna/hg19"
resolutions = ["2000k", "1000k", "500k", "200k"]

# set up global loaders dictionary to stuff loaders into
Loaders = {}

def getFileName(cellType, replicate, resolution):
    """
    return the heat map file name for a given
    replicate and resolution.
    """
    assert resolution in resolutions
    tmpl = os.path.join(path, "{0}-{1}-HindIII-{2}.hm")
    fname = tmpl.format(cellType, replicate, resolution)
    return fname

def getLoader(res):
    """
    Memoize the loader so that we are not doing many read/writes
    """
    if not Loaders.has_key(res):
        print("Initializing new loader for {0}.".format(res))
        Loaders[res] = binnedData(res, genome=genome)
    return Loaders[res]

def adjustSpines(ax):
    """pushes spines away from the plot axes"""

    # Move left and bottom spines outward by 10 points
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))

    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    return

def plotCheckerboardHeatMap(cellType, replicate, resolution, chromosome):
    """
    Plots a heatmap for the IMR90 cell line
    at a given resolution and repliate number
    """
    print("Loading Heatmap Data for chromosome %i...", chromosome)

    res = nutils.strToResolution(resolution)
    loader = getLoader(res)

    # decriment, to match params
    chromosome -= 1

    inFile = getFileName(cellType, replicate, resolution)
    key = "%s-%s" % (replicate, resolution)

    loader.simpleLoad(inFile, key)

    data = nutils.filterByChromosome(loader, chromosome, key=key)

    print("Plotting Figure %i...", chromosome)
    fig, ax = plt.subplots()

    # Calculate the checkerboard pattern by dividing
    # through by the average.
    avg = np.average(data)
    normed = np.log(np.divide(data, avg))

    ax.imshow(normed, cmap=plt.cm.seismic)

    # Chromosome index needs to be human readable!
    label = chromosome + 1

    ax.set_title('Chromosome {0} Checkerboard ({1})'.format(label, resolution))
    ax.set_xlabel('Chromosome {0}'.format(label))
    ax.set_ylabel('Chromosome {0}'.format(label))

    savePath = "pngs/br-{0}-{1}-{2}-chr{3}.png".format(cellType, replicate, resolution, label)
    print("Saving figure to ", savePath)
    fig.savefig(savePath, dpi=850)

    # explicitly close the figure to conserve memory
    plt.close()
    return

def main(cellType, resolution):
    """
    Runs the main module code
    """
    # validity checks
    assert resolution in resolutions
    try:
        experiments = replicates[cellType]
    except:
        raise Exception("No experiments for cell type %s" % cellType)

    # plot heatmaps
    for replicate in experiments:
        for chrom in range(1, 5):
            plotCheckerboardHeatMap(cellType, replicate, resolution, chrom)

if __name__ == "__main__":
    main("IMR90", "1000k")
    main("hESC", "1000k")
