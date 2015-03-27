#!usr/bin/env python
"""
plot.py
This script plots number of interaction probes as a function of distance.

Copyright (C) 2015 Jonathan Niles

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import pickle
import nutils as nu
import numpy as np
import matplotlib.pyplot as plt
from hiclib.binnedData import binnedData

# switch this from ic to raw if desired
etype = "ic"
path = "/home/jniles/data/{0}/%s/{0}-{1}-HindIII-1000k.hm" % etype
savePath = nu.chkdir(nu.join(nu.sync, "data/probes/"))
figPath = nu.chkdir(nu.join(nu.sync, "plots/probes/"))
genome = "/home/jniles/data/dna/hg19/"

replicates = {
    "IMR90" : ["R1", "R2", "R3", "R4", "R5", "R6"],
    "hESC" : ["R1", "R2"]
}

def loadData(cellType):
    """
    Loads all the data for a given cell type into a binnedData instance for
    analysis.
    """
    BD = binnedData(1000000, genome=genome)
    datasets = replicates[cellType]
    for ds in datasets:
        BD.simpleLoad(path.format(cellType, ds), ds)
    return BD

def calculateScaling(cellType, loader):
    """
    Plots scaling for all chromosome keys
    """
    # output dictionary decribing scaling factors
    scaling = {}

    for key, data in loader.dataDict.items():
        scaling[key] = {}
        for idx in xrange(23):
            wdata = nu.filterByChromosome(loader, idx, key=key)
            length = len(wdata)
            contacts = np.zeros(length)

            # O(n**2)
            # scan the entire binned chromosome
            # for each bin, calculate the distance between
            # it and all other bins, and add probes to a global
            # probe count for each distance.
            # The matrix is symmetric, so we only need to scan
            # the bottom or top diagonal
            for i in xrange(length):
                for j in xrange(i):
                    # absolute distance in bins
                    distance = np.abs(i - j)
                    # if we are not in the same bin,
                    # add the number of probes in that bin
                    # to the total number of probes at that distance
                    if distance != 0:
                        contacts[distance] += wdata[i, j]
            scaling[key][idx] = contacts

    outFile = nu.join(savePath, '{0}-{1}-{2}.pkl'.format(etype, cellType, key))
    return saveScalingData(scaling, outFile)

def saveScalingData(scale, fname):
    """exports the scaling dictionary to a python pickle"""
    with open(fname, 'w') as outfile:
        pickle.dump(scale, outfile)
    return scale

def loadScalingData(cellType):
    """
    Loads scaling by replicate from preprocessed pickle data
    """
    datasets = replicates[cellType]
    dataDict = {}
    for rep in datasets:
        inFile = nu.join(savePath, '{0}-{1}-{2}.pkl'.format(etype, cellType, rep))
        with open(inFile, 'r') as f:
            dataDict[rep] = pickle.load(f)
    return dataDict

def plotScaling(cellType, dataDict, percent=False):
    """
    Plots the scaling plot for each chromosome by replicate
    """
    for rep in dataDict.keys():
        r = dataDict[rep]
        fig, ax = plt.subplots()

        # loop through chromosomes, plotting on the same plot
        for chrom in r.keys():
            contacts = r[chrom]
            length = len(contacts)

            # get an aggregate array of the contacts
            norm = np.cumsum(contacts) / np.sum(contacts)

            x = np.arange(length)
            if percent:
                x = np.linspace(0, 1, length)

            # plot on the same axis
            ax.plot(x, norm, label="Chr{0}".format(chrom))

        # legends
        ax.set_title('{0} {1} Contact Scaling'.format(cellType, rep))

        # labels
        if percent:
            ax.set_xlabel('Percentage of Chromosome')
        else: 
            ax.set_xlabel('Distance (Mb)')
        ax.set_ylabel('Percentage of Interations')
        ax.set_ylim(0, 1.1)
        ax.grid()
        ax.minorticks_on()

        # make sure the legend is properly placed
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)

        # save and close
        extension = '{0}-{1}-{2}.png'.format(etype,cellType, rep)

        # make sure the dir exists
        if percent:
            folder = nu.chkdir(nu.join(figPath, "percent/"))
        else:
            folder = nu.chkdir(nu.join(figPath, "distance/"))

        outFig = nu.join(folder, extension)

        fig.savefig(outFig, dpi=500)
        plt.close()
    return

def plotAll():
    """
    Runs the module
    """
    for cellType in replicates.keys():

        # make the save directory if it doesn't exist
        mtx = loadData(cellType)
        scaling = calculateScaling(cellType, mtx)
        plotScaling(cellType, scaling, percent=True)
        plotScaling(cellType, scaling, percent=False)
    return

if __name__ == "__main__":
    plotAll()
