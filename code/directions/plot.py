#!/usr/bin/env python
"""
plot.py
plots all directionality index for a given chromosome

Copyright (C) 2015  Jonathan Niles

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

from __future__ import print_function

import itertools
import numpy as np
import nutils as nu
import matplotlib.pyplot as plt
from matplotlib.mlab import csv2rec

# data directories
bedDir = nu.chkdir(nu.join(nu.sync, "data/domains/bed"))
diDir = nu.chkdir(nu.join(nu.sync, "data/di/"))
pltDomDir = nu.chkdir(nu.join(nu.sync, "plots/domains/hist/"))
pltDIDir = nu.chkdir(nu.join(nu.sync, "plots/di/index/"))

# templates and paths
btmpl = "{0}-{1}-{2}kb.bed"
dtmpl = "{0}-{1}-{2}kb.npy"
columns = ['chrom', 'start', 'end', 'value', 'sign', 'length', 'uuid']
windows = [100, 200, 400, 800, 1000]
#windows = [200, 400, 800, 1000]

# flattens arrays of arrays to a single array
squash = lambda arr: np.array(list(itertools.chain(*arr)))

def parseBed(cellType, replicate, windowSize):
    """
    parses a bed file given in the template
    """
    fname = nu.join(bedDir, btmpl.format(cellType, replicate, windowSize))
    return csv2rec(fname, delimiter='\t', names=columns)

def loadDI(cellType, replicate, windowSize):
    """
    loads the directionality index for a given cellType, replicate at
    a certain resolution
    """
    fname = nu.join(diDir, dtmpl.format(cellType, replicate, windowSize))
    return np.load(fname)

def plotDomainSizeHistogram(cellType, replicate, windowSize, bySign=False):
    """
    plots a histogram of the domain lengths for a replicate
    """
    print("Plotting domain size histogram for %s %s" % (cellType, replicate))

    domains = parseBed(cellType, replicate, windowSize)

    fig, ax = plt.subplots()

    # seperate by domain directionality
    if bySign:
        posDomains = domains[domains['sign'] > 0]['length'] / 1000
        negDomains = domains[domains['sign'] <= 0]['length'] / 1000
        ax.hist(negDomains, color='r', alpha=0.5, bins=50, label='Negative Domains')
        ax.hist(posDomains, color='g', alpha=0.5, bins=50, label='Positive Domains')
    else:
        ax.hist(domains['length'] / 1000, color='b', alpha=0.5, bins=50, label='Domains')

    median = np.median(domains['length'])

    label = 'Median: ~{0:.2f}kb'.format(median / 1000.0) # format to kb
    fig.text(0.75, 0.75, label, horizontalalignment='center', verticalalignment='center')

    # make sure the axis goes to zero
    x1,x2 = ax.get_xlim()
    ax.set_xlim((0, x2))
    ax.set_ylim((0, 1000)) # limit to 1000

    ax.set_title('{0} {1} {2}kb Windows @ {0}'.format(cellType, replicate, windowSize, "10kb"))
    ax.set_xlabel('Domain Size (kb)')
    ax.set_ylabel('# of Domains')

    # make sure everything is visible
    plt.tight_layout()
    plt.legend()

    # save figure
    if bySign:
        outDir = nu.chkdir(nu.join(pltDomDir, "signed/"))
    else:
        outDir = pltDomDir

    fname = nu.join(outDir, "{0}-{1}-{2}kb.png".format(cellType, replicate, windowSize))
    fig.savefig(fname, dpi=600)
    plt.close()
    return

def plotDirectionalityIndex(cellType, replicate, windowSize, norm=False):
    """
    Plots the directionality index for a cell type and replicate
    TODO: Think about making the x labels based on chromosome rather than on position
    """

    print("Plotting directionality index for %s %s" %(cellType, replicate))

    # load in the data and flatten it into a single array
    data = squash(loadDI(cellType, replicate, windowSize))
    length = len(data)

    # normalize the data if required
    if norm:
        data = data / np.max(np.abs(data))

    fig, ax = plt.subplots(figsize=(18, 4))

    x = np.arange(len(data)) * 10
    ax.plot(x, data, 'k-', linewidth=0.5)

    plt.xlabel('hg19 (10kb)')
    plt.ylabel('DI')
    plt.title('{0} {1} Directionality Index @ {2}kb'.format(cellType, replicate, windowSize))

    # format the axes nicely
    ax.minorticks_on()
    ylimit = np.max(np.abs(ax.get_ylim()))
    ax.set_ylim((-ylimit, ylimit)) # make y axis symmetric
    ax.set_xlim((0, length))

    # make sure everything is visible
    plt.tight_layout()

    # save the figure
    fname = nu.join(pltDIDir, "{0}-{1}-{2}kb.png".format(cellType, replicate, windowSize))
    fig.savefig(fname, dpi=800)

    plt.close()
    return

def plotOverlapHistogram():
    """
    plots overlaps
    """
    cols = ['chromA', 'startA', 'endA', 'valueA', 'signA', 'lengthA', 'idA', 'chromB',
            'startB', 'endB', 'valueB', 'signB', 'lengthB', 'idB', 'distance']

    fig, ax = plt.subplots()

    array = csv2rec('./hESC-20kb.bed', delimiter='\t', names=cols)
    distances = array['distance']
    ax.hist(distances, bins=50)
    
    median = np.median(distances)
    label = 'Median: ~{0:.2f}bp'.format(median) # format to kb
    fig.text(0.75, 0.75, label, horizontalalignment='center', verticalalignment='center')

    # format axes
    ax.set_title("Domain Separation hESC @ 10Kb")
    ax.set_ylabel("Counts")
    ax.set_xlabel("Distance (kb)")
    ax.set_xticklabels(ax.get_xticks() / 1000)

    fig.savefig(dpi=750)
    return

# TODO
def plotDomainBoundaries(data, domains, chrom=0):
    """
    Plot these regions on a graph beautifully.

    :param data: the directionality index
    :param domains: a list of domains
    :param chrom: chromosome number
    :returns: None
    """
    fig, ax = plt.subplots()

    ax.plot(data, '-k', linewidth=0.5)

    for x, y in domains:
        ax.axvline(x[0], linewidth=0.25, color='r')
        ax.axvline(x[2], linewidth=0.25, color='r')

    ax.set_title('Domain Boundaries')
    ax.set_xlabel('')

    fig.savefig("out.png", dpi=600)
    plt.close()
    return

def plotAll():
    for win in windows:
        for cellType in nu.datasets.keys():
            for rep in nu.datasets[cellType]:
                plotDomainSizeHistogram(cellType, rep, win, bySign=False)
                plotDirectionalityIndex(cellType, rep, win, norm=True)

def plotCorrelationsByWindowSize():
    """plots the mean and min correlations by window size"""
   
    print("Plotting correlations by Window Size")
    dataDir = nu.join(nu.sync, "data/di/corr/")
    dpath = "{0}kb.window.npytxt"
    data = []
    labels = []
    for win in windows:
        data.append(np.loadtxt(nu.join(dataDir, dpath.format(win))))
        labels.append("{0}kb".format(win))

    mins = map(np.min, data)
    means = map(np.mean, data)

    fig, ax = plt.subplots()

    x = np.arange(len(windows))
    ax.plot(x, mins, "o-g", label="min")
    ax.plot(x, means, "o-y", label="mean")
    
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Window Sizes")
    ax.set_ylabel("Spearman's $\\rho$")
    ax.set_title("DI Correlation by Window Size")
    ax.set_ylim(0, 1)
    ax.grid()

    ax.legend(loc='best')

    outFig = nu.join(nu.sync, "plots/di/correlationsByWindow.png")
    print("Saving figure to", outFig)
    fig.savefig(outFig, dpi=450)
    plt.close()
    return


if __name__ == "__main__":
    #plotAll()
    plotCorrelationsByWindowSize()
