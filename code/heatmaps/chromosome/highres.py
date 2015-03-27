#!/usr/bin/env python
# plot a high-res plot of a single chromosome

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import nutils as nu
from hiclib.highResBinnedData import HiResHiC

HG19 = '~/data/dna/hg19/'
RESOLUTION = 10000

# To be formatted
FILE = '~/data/{0}/ic/{0}-all-HindIII-10k_HighRes.byChr'

Loader = None

def parseChromosomeNumber(chromStr):
    """return the tuple idx for the chromosome"""
    i = int(chromStr.strip('chr')) - 1
    return (i, i)

def loadHiResData(chromIdx, ct):
    """ This function caches the matrix to avoid lots of i/o """
    try:
        data = Loader.data[chromIdx].getData()
        print("Using loaded data...")
    except:
        print("Loading fresh data...")
        global Loader
        Loader = HiResHiC(HG19, RESOLUTION)
        Loader.loadData(FILE.format(ct, ct), mode='cis')
        print "Loaded data. Getting chromosome matrix."
        data = Loader.data[chromIdx].getData()
    return data

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

# This is will run the entire module
def main(chromosome, ct):
    """plot a chromosome at 10kb resolution"""

    caps = lambda x : x.capitalize()

    # parse the chromosome number
    chromIdx = parseChromosomeNumber(chromosome)
    print("Calling loadHighResData with", chromIdx)
    data = loadHiResData(chromIdx, ct)
    logData = np.log(data)

    print("Plotting the matrix...")
    fig, ax = plt.subplots()

    # plot the matrix
    ax.imshow(logData, interpolation='None', cmap=plt.cm.Reds)

    fig.suptitle(caps(ct))
    ax.set_title('Chromosome {0} highres'.format(chromIdx[0]))
    ax.set_xlabel(caps(chromosome))
    ax.set_ylabel(caps(chromosome))

    # make sure the directory exists
    savePath = os.path.join(nu.chkdir('highres/%s/'%ct), '%s-high-res.png'%chromosome)

    # save the figure
    fig.savefig(savePath, dpi=800)
    plt.close()
    return

if __name__ == "__main__":
    for ct in ["IMR90", "hESC"]:
        for chrom in ['chr%i' % x for x in xrange(1,4)]:
            main(chrom, ct)
