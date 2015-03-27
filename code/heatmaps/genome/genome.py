#!/usr/bin/env python

import os
import sys
import nutils
import numpy as np
import matplotlib.pyplot as plt
from hiclib.binnedData import binnedData

# cellType, {raw | ic}, cellType, replicate, resolution
path = "/home/jniles/data/{0}/{1}/{0}-{2}-HindIII-{3}.hm"
genome = "/home/jniles/data/dna/hg19/"

# useful parameters
folder = "correct"
genome = "/home/jniles/data/dna/hg19"
dataSets = {
    "IMR90" : ["R1", "R2", "R3", "R4", "R5", "R6"],
    "hESC" : ["R1", "R2"]
}


def addChromosomeTickMarks(bd, ax):
    """
    removes old ticks and puts tickmarks at the center of chromosomes
    """
    # Finds the midpoint of each chromosome in the dataset
    # Lambda finds the midpoint [(y-x)/2.] and adds it to the previous
    # chromosomeStarts position
    tickPos = map((lambda x,y: (y-x)/2.+x), bd.chromosomeStarts, bd.chromosomeEnds)

    # create chromosome labels on the fly
    tickLabels = map(lambda i: "chr%i" % i, xrange(1,22))
    tickLabels.append('chrX')

    # add the first and last labels to the plot
    ax.set_xticks([tickPos[0], tickPos[-1]])
    ax.set_yticks([tickPos[0], tickPos[-1]])

    ax.set_xticklabels([tickLabels[0], tickLabels[-1]])
    ax.set_yticklabels([tickLabels[0], tickLabels[-1]])
    return

def plotGenomeHeatmap(cellType, replicate, resolution, red=False, dType="raw", cTicks=True):
    """
    plots a genome heatmap
    """
    print "{0} {1} {2} {3}".format(cellType, dType, replicate, resolution)
    res = nutils.strToResolution(resolution)
    bd = binnedData(res, genome)
    
    inFile = path.format(cellType, dType, replicate, resolution)
    bd.simpleLoad(inFile, "key")
    data = bd.dataDict["key"]

    # plot
    fig, ax = plt.subplots()
    if red:
        ax.imshow(np.log(data), cmap=plt.cm.Reds, interpolation='None')
    else:
        ax.imshow(np.log(data))

    # add tick marks if required
    if cTicks: 
        addChromosomeTickMarks(bd, ax)

    # labeling
    ax.set_title("{0} {1} {2} {3}".format(cellType, dType, replicate, resolution))

    savePath = nutils.chkdir("{0}/{1}/".format(dType, cellType))
    outPng = os.path.join(savePath, "{0}-{1}.png".format(replicate, resolution))
    fig.savefig(outPng)

    # explicitly close the figure to conserve memory
    plt.close()
    return

# This is will run the entire module
if __name__ == "__main__":
    for cellType in ["hESC"]:
        for rep in dataSets[cellType]:
            for res in ["1000k", "200k"]:
                plotGenomeHeatmap(cellType, rep, res)
                plotGenomeHeatmap(cellType, rep, res, dType="ic")
