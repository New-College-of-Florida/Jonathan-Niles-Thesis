#!/usr/bin/env python
"""
All heatmap type plots
"""

from __future__ import print_function

import numpy as np
import nutils as nu
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from hiclib.binnedData import binnedData

# Constants
genome = "/home/jniles/data/dna/hg19"
resolutions = ["1000k", "500k", "200k"]
gpath = "/home/jniles/data/{0}/{1}/{0}-{2}-HindIII-{3}.hm"

def plotCorrelationsByResolution():
    """Plots the correlations between heat maps by resolution
    for hESC data"""

    cellType = "hESC"
    replicates = nu.datasets[cellType]
    rawCorrelations = []

    # load raw heatmaps first
    for res in resolutions:
        print("Initializing new raw loader for {0}".format(res))
        r = nu.strToResolution(res)
        loader = binnedData(r, genome=genome)

        for i, rep in enumerate(replicates):
            path = gpath.format(cellType, "raw", rep, res)
            loader.simpleLoad(path, i)

        matrices = loader.dataDict.values()
        rho, pval= spearmanr(matrices[0].ravel(), matrices[1].ravel())
        rawCorrelations.append(rho)
        matrices = loader.dataDict.values()
        print("Raw Correlations calculated: rho=", rho, "p=", pval)

    correctedCorrelations = []

    # load corrected heatmaps next
    for res in resolutions:
        print("Initializing new loader for {0}.".format(res))
        r = nu.strToResolution(res)
        loader = binnedData(r, genome=genome)

        for i, rep in enumerate(replicates):
            path = gpath.format(cellType, "ic", rep, res)
            loader.simpleLoad(path, i)

        matrices = loader.dataDict.values()
        rho, pval= spearmanr(matrices[0].ravel(), matrices[1].ravel())
        correctedCorrelations.append(rho)
        print("Corrected Correlations calculated: rho=", rho, "p=", pval)

   
    fig, ax = plt.subplots()
    N = len(resolutions) + 1
    x = np.arange(1,N)
    ax.plot(x, rawCorrelations, "y-o", label="raw")
    ax.plot(x, correctedCorrelations, "r-o", label="corrected")
    ax.grid()
   
    ax.xlim(0, N)
    ax.set_xlabel("Replicates")
    ax.set_ylabel("Spearman's Correlation")

    outFig = nu.join(nu.sync, "plots/heatmaps/correlations.png")
    print("Saving figure to", outFig)
    fig.savefig(outFig, dpi=700)
    return

if __name__ == "__main__":
    plotCorrelationsByResolution()
