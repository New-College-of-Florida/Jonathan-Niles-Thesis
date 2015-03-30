#!/usr/bin/env python

import os
import sys

import matplotlib
matplotlib.use('agg')

import numpy as np
import matplotlib.pyplot as plt

from mirnylib import h5dict

# Constants
MODULE = "TRANS-HEATMAP"

dataPath = "/home/jniles/data/imr90/raw"
genome = "/home/jniles/data/dna/hg19"
cellType = "IMR90"
enzyme = "HindIII"
resolutions = ["2000k", "1000k", "500k", "200k"]
replicates = ["all"]

def fmt(msg, *args):
    print "[%s] %s" % (MODULE, msg % args)

def getFileName(replicate, resolution):
    """
    return the heat map file name for a given
    replicate and resolution
    """
    assert replicate in replicates
    assert resolution in resolutions
    fname = "{0}-{1}-{2}-{3}.hm".format(cellType, replicate, enzyme, resolution)
    return os.path.join(dataPath, fname)

def getChromMask(data, chrm):
    return np.array([i == chrm for i in data['chromosomeIndex']])

def getChromSlice(data, chrm):
    mask = getChromMask(data, chrm)
    return data['heatmap'][mask]

def getChromSquareSlice(data, chrm):
    mask = getChromMask(data, chrm)
    return data['heatmap'][mask].T[mask]

def getChromTransInteractions(data, chrm1, chrm2):
    m1 = getChromMask(data, chrm1)
    m2 = getChromMask(data, chrm2)
    return data['heatmap'][m1].T[m2]

def makeChromosomeInteractionMatrix(data):
    matrix = np.zeros((23,23))
    for y in xrange(23):
        for x in xrange(23):
            if x != y:
                matrix[x][y] = np.sum(getChromTransInteractions(data, x, y))

    return matrix

def getSizeNormalization(data, chrN, chrM):
    sizes = data['chromosomeIndex']
    n = len(sizes[sizes == chrN])
    m = len(sizes[sizes == chrM])
    return n*m

def arrayCount(value, array):
    return np.sum(map(lambda x: int(x == value), array))

def plotGenomeTransInteractions(data, saveFig=False):
    # `data` is the hdf5 file
    # This should be makeChromosomeInteractionMatrix() if need be
    matrix = makeChromosomeInteractionMatrix(data)

    # Strip out the cis interactions
    for x in xrange(23):
        for y in xrange(23):
            if x != y:
                matrix[x][y] = matrix[x][y] / getSizeNormalization(data, x, y)
            else:
                matrix[x][y] = 0

    plt.matshow(np.log(matrix), cmap=plt.get_cmap('seismic'))
    plt.colorbar()

    labels = map(lambda i: str(i), range(1,23))
    labels.append('X')

    ax = plt.gca()
    plt.xticks(xrange(23))
    plt.yticks(xrange(23))
    ax.set_xticklabels([])
    ax.set_yticklabels(labels)

    plt.title("Trans Interaction Enrichment")
    plt.xlabel("Chromosomes")
    plt.ylabel("Chromosomes")

    if saveFig:
        plt.savefig('TransSizeLogNormalization.png', dpi=900, figsize=(25,25))

    return

def main():
    data = h5dict.h5dict(getFileName("all", "1000k"))
    plotGenomeTransInteractions(data, saveFig=True)

if __name__ == "__main__":
    main()
