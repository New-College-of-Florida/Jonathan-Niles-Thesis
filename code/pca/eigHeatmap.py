#!/usr/bin/env python

import os
import nutils
import matplotlib.pyplot as plt
import numpy as np
from hiclib.binnedData import binnedData

# path template
tmpl = "/home/jniles/data/{0}/ic/{0}-all-HindIII-200k.hm"
genome = "/home/jniles/data/dna/hg19"
chroms = range(12)

# datasets
dataSets = {
    "IMR90": {
        "evectors" : "IMR90-200k.eigenvectors",
        "evalues" : "IMR90-200k.eigenvalues",
        "hm" : tmpl.format("IMR90"),
        "save" : "heatmap-chr{0}-imr90-200k.png"
    },
    "hESC" : {
        "evectors" : "hESC-200k.eigenvectors",
        "evalues" : "hESC-200k.eigenvalues",
        "hm" : tmpl.format("hESC"),
        "save" : "heatmap-chr{0}-hESC-200k.png"
    }
}

def plotAll(evector, data, save, name=""):
    """Plots a form of scree plot for eigenvalues"""

    # heatmap and eigen axes
    axh = plt.subplot2grid((5,1), (0,0), rowspan=4)
    axe = plt.subplot2grid((5,1), (4,0), sharex=axh)
    fig = plt.gcf()

    # the eigenvalues can be negative, so plot the magnitude of the eigenvalue
    axe.plot(np.arange(len(evector)), evector)

    # pretty axis limits
    axe.set_xlabel("bin")
    axe.set_ylabel("$\\vec{v}$")

    # heatmap plotting
    axh.imshow(np.log(data), cmap=plt.cm.Reds, interpolation="none")
    axh.set_adjustable('box-forced')
    axh.set_title(name)

    plt.tight_layout()

    fig.suptitle("Comparison of 1st Eigenvector and Compartments")

    fig.savefig(os.path.join("plots", save), dpi=600)

    plt.close()
    return

if __name__ == "__main__":
    print "Plotting Eigenvector vs. Heatmap Plots"

    for key, sets in dataSets.items():
        vectors = np.loadtxt(sets["evectors"])
        values = np.loadtxt(sets["evalues"])

        # load data
        BD = binnedData(200000, genome=genome)
        BD.simpleLoad(sets["hm"], key)

        # get the specific chromosome
        for ch in chroms:
            data = nutils.filterByChromosome(BD, key, ch)
            n = "{0}-chr{1}".format(key, ch+1)
            plotAll(vectors[0], data, sets["save"].format(ch+1), name=n)
