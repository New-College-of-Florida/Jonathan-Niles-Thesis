#!/usr/bin/env python
"""
All heatmap type plots
"""

from __future__ import print_function

import numpy as np
import nutils as nu
import matplotlib.pyplot as plt
from hiclib.binnedData import binnedData
from hiclib.highResBinnedData import HiResHiC

outDir = nu.chkdir(nu.join(nu.sync, "heatmaps/"))

# Constants
cpath = "/home/jniles/data/{0}/ic/"
rpath = "/home/jniles/data/{0}/raw/"
genome = "/home/jniles/data/dna/hg19"
#resolutions = ["2000k", "1000k", "500k", "200k"]
resolutions = ["1000k", "200k"]

# Data loader
Loaders = {}
HiResLoader = None

def lowResData(res):
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

"""
High Res Heatmap Code
"""

def loadHiResData(cellType, chrom, rep="R1" ):
    """ This funcellTypeion caches the matrix to avoid lots of i/o
    
    :param chrom: chromosome number
    :param cellType: cell type name
    
    """
    global HiResLoader
    #resolution = 10000
    resolution = 20000
    #fname = "/home/jniles/data/{0}/ic/{0}-{1}-HindIII-10k_HighRes.byChr"
    fname = "/home/jniles/data/{0}/ic/{0}-{1}-HindIII-20k_HighRes.byChr"
    try:
        data = HiResLoader.data[(chrom,chrom)].getData()
    except:
        HiResLoader = HiResHiC(genome, resolution)
        HiResLoader.loadData(fname.format(cellType, rep), mode='cis')
        data = HiResLoader.data[(chrom, chrom)].getData()
    return data

def highResHM(cellType, chrom, rep="R1"):
    """plots a chromosome heatmap at high resolution (10kb)

    :param chrom: chromosome number to plot
    :param cellType: cell type to load
    :returns: None
    """

    saveDir = nu.chkdir(nu.join(outDir, 'highres/{0}/'.format(cellType)))

    # make sure chromosomes are labeled properly
    if chrom == 23:
        chromLabel = "ChrX"
    else:
        chromLabel = "Chr{0}".format(chrom+1)

    data = loadHiResData(cellType, chrom, rep=rep)

    fig, ax = plt.subplots()

    # plot the matrix
    print("Plotting high resolution heatmap {0} {1}".format(cellType, chromLabel))
    ax.imshow(np.log(data), interpolation='None', cmap=plt.cm.Reds)

    fig.suptitle("{0} {1}".format(cellType, rep))
    #ax.set_title('{0} 10kb Resolution'.format(chromLabel))
    ax.set_title('{0} 20kb Resolution'.format(chromLabel))
    ax.set_xlabel(chromLabel)
    ax.set_ylabel(chromLabel)

    # save the figure
    #figPath = "{0}-{1}-{2}.png".format(cellType, rep, chromLabel)
    figPath = "{0}-{1}-{2}-20kb.png".format(cellType, rep, chromLabel)
    fig.savefig(nu.join(saveDir, figPath), dpi=800)
    plt.close()
    return

def plotAllHighResHeatmaps():
    """
    plots all the high res heatmaps you could ask for
    """
    for cellType in nu.datasets.keys():
        for chrom in xrange(2,20):
            highResHM(cellType, chrom)

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

"""
Whole Genome HeatMap
"""

def plotGenomeHeatmap(cellType, replicate, resolution, red=False, src="raw"):
    """plots a genome heatmap"""
    print("Plotting whole genome heatmap: {0} {1} {2} {3}".format(cellType, src, replicate, resolution))

    res = nu.strToResolution(resolution)
    bd = binnedData(res, genome)
   
    path = "/home/jniles/data/{0}/{1}/{0}-{2}-HindIII-{3}.hm"
    inFile = path.format(cellType, src, replicate, resolution)
    bd.simpleLoad(inFile, "key")
    data = np.log(bd.dataDict["key"])

    # plot
    fig, ax = plt.subplots()
    ax.set_title("{0} {1} {2} {3}".format(cellType, src, replicate, resolution))

    if red:
        ax.imshow(data, cmap=plt.cm.Reds, interpolation='None')
    else:
        ax.imshow(data, interpolation='None')

    addChromosomeTickMarks(bd, ax)

    # save the figure
    savePath = nu.chkdir(nu.join(outDir, "genome/{0}/".format(src)))
    outPng = nu.join(savePath, "{0}-{1}-{2}.png".format(cellType, replicate, resolution))
    fig.savefig(outPng, dpi=850)

    # explicitly close the figure to conserve memory
    plt.close()
    return

def plotGenomeHeatmaps():
    for cellType in nu.datasets.keys():
        for rep in nu.datasets[cellType]:
            for res in resolutions:
                plotGenomeHeatmap(cellType, rep, res, src="raw")
                plotGenomeHeatmap(cellType, rep, res, src="ic")

"""
Checkerboard Heat Map
"""

def plotCheckerboardHeatMap(cellType, replicate, resolution, chrom, src="raw"):
    """
    Plots a heatmap for cell line at a given resolution and replicate number
    """
    print("Loading data for chromosome {0}...".format(chrom))

    if src == "raw":
        path = nu.join(rpath, "{0}-{1}-HindIII-{2}.hm")
    else:
        path = nu.join(cpath, "{0}-{1}-HindIII-{2}.hm")

    # convert string resolution if necessary
    if type(resolution) == type("S"):
        res = nu.strToResolution(resolution)
    else:
        res = resolution

    loader = lowResData(res)
    inFile = path.format(cellType, replicate, resolution)
    key = "%s-%s" % (replicate, resolution)
    loader.simpleLoad(inFile, key)

    data = nu.filterByChromosome(loader, chrom - 1, key=key)

    label = "Chr{0}".format(chrom)

    print("Plotting checkerboard heatmap for {0} {1} {2}...".format(cellType, replicate, label))
    fig, ax = plt.subplots()

    # Calculate the checkerboard pattern by dividing
    # through by the average.
    avg = np.average(data)
    data = np.log(np.divide(data, avg))

    ax.imshow(data, interpolation='None', cmap=plt.cm.seismic)

    ax.set_title('{0} Checkerboard ({1})'.format(label, resolution))
    ax.set_xlabel(label)
    ax.set_ylabel(label)

    saveDir = nu.chkdir(nu.join(outDir, "checkerboard/"))
    savePath = nu.join(saveDir, "{0}-{1}-{2}-{3}.png".format(cellType, replicate, resolution, label))

    fig.savefig(savePath, dpi=850)

    # explicitly close the figure to conserve memory
    plt.close()
    return

def plotCheckboardHeatmaps():
    """plots all checkerboard heatmaps"""
    for cellType in nu.datasets.keys():
        for rep in nu.datasets[cellType]:
            for res in resolutions:
                for i in xrange(1,20):
                    plotCheckerboardHeatMap(cellType, rep, res, i)
                    plotCheckerboardHeatMap(cellType, rep, res, i, src="ic")
    return

if __name__ == "__main__":
    plotAllHighResHeatmaps()
    #plotGenomeHeatmaps()
    #plotCheckboardHeatmaps()
