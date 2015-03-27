#!/usr/bin/env python
"""
plot.py
plots the probes/genes expression values from the Affy dataset

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

from __future__ import print_function

import itertools
import numpy as np
import nutils as nu
from matplotlib import mlab
import matplotlib.pyplot as plt
from matplotlib.mlab import csv2rec
from scipy.stats import spearmanr

tmpl = nu.join(nu.sync, "data/components/{0}-{1}-{2}.eigenvectors")
volcanoDir = nu.chkdir(nu.join(nu.sync, "plots/genes/"))
outDir = nu.chkdir(nu.join(nu.sync, "plots/genes/compartment/"))
histDir = nu.chkdir(nu.join(nu.sync, "plots/genes/histograms/"))
resMap = {"200k":200000, "1000k":1000000}
squash = lambda arr: np.array(list(itertools.chain(*arr)))

def loadRawArray():
    """loads raw expression array"""
    fname = nu.join(nu.sync, 'data/genes/expression.bed.gz')
    return csv2rec(fname, delimiter='\t')

def probesByCompartment(cellType, rep, res, n=0):
    """plots probe expression values and their PCA compartment value"""
    data = loadRawArray()
    components = np.loadtxt(tmpl.format(cellType, rep, res))
    pc = components[n] # load the nth component
    
    # create columns for quick iteration
    imr90Cols = ["imr90_{0}".format(i) for i in xrange(1,3)]
    hescCols = ["hesc_{0}".format(i) for i in xrange(1,3)]
    imr90Means = np.array([np.mean([np.log(i),np.log(j)]) for i,j in data[imr90Cols]])
    hescMeans= np.array([np.mean([np.log(i),np.log(j)]) for i,j in data[hescCols]])
    delta = imr90Means - hescMeans

    # get the numerical resolution
    nres = resMap[res]
    starts = data['start'] / nres
    ends = (data['end'] / nres) + 1

    compartments = np.array([np.mean(pc[start:end]) for start, end in zip(starts, ends)])
    
    fig, ax = plt.subplots()
    ax.plot(compartments, delta, 'og', alpha=0.3, markersize=0.3)

    plt.suptitle("{0} {1} {2}".format(cellType, rep, res))

    ax.set_title("Compartment Expression (pc={0})".format(n+1))
    ax.set_xlabel("Compartment Character")
    ax.set_ylabel("Probe Log2 Expression Change")
    ax.set_ylim((-8, 8))

    corr = spearmanr(compartments, delta)[0]

    s = "Spearman = {0:0.3f}".format(corr)
    ax.text(0.9, 0.9, s,
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax.transAxes)

    ax.legend()
    ax.minorticks_on()

    fname = nu.join(outDir, "{0}-{1}-{2}.png".format(cellType, rep, res))
    print("Saving figure as...", fname)
    fig.savefig(fname, dpi=800)
    plt.close()
    return

def cellTypeExpressionHistogram(cellType):
    """plots the gene expression for a given cell type"""
    data = loadRawArray()
    cols = ["{0}_{1}".format(cellType.lower(), i) for i in xrange(1,3)]
    means = np.array([np.mean([np.log2(i),np.log2(j)]) for i,j in data[cols]])

    print("Plotting expression histogram for", cellType)

    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(means, bins=250, histtype='stepfilled',
            alpha=0.75, label="Genes", normed=True)

    # plot fitted normal distribution
    mean, stddev = np.mean(means), np.std(means)
    text = "N($\mu=%.2f$, $\sigma=%.2f$)" % (mean, stddev)
    normal = mlab.normpdf(bins, mean, stddev)
    ax.plot(bins, normal, 'k--', label=text)

    ax.set_title("H1/IMR90 Gene Expression Changes")
    ax.set_xlabel('Expression Level Log2')
    ax.set_ylabel('Frequency')

    ax.legend()
    ax.minorticks_on()

    fname = nu.join(histDir, "{0}-expr.png".format(cellType))
    print("Saving to", fname)
    fig.savefig(fname, dpi=600)
    plt.close()
    return

def expressionChangeHistogram(): 
    """plots the gene expression changes between cell types"""
    # we are going to do the hESC experiment plots first
    data = loadRawArray()
    
    # create columns for quick iteration
    imr90Cols = ["imr90_{0}".format(i) for i in xrange(1,3)]
    hescCols = ["hesc_{0}".format(i) for i in xrange(1,3)]
    imr90Means = np.array([np.mean([np.log2(i),np.log2(j)]) for i,j in data[imr90Cols]])
    hescMeans= np.array([np.mean([np.log2(i),np.log2(j)]) for i,j in data[hescCols]])
    delta = imr90Means - hescMeans

    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(delta, bins=250, histtype='stepfilled',
            alpha=0.75, label="Expression Change")

    # plot fitted normal distribution
    mean, stddev = np.mean(delta), np.std(delta)
    text = "N($\mu=%.2f$, $\sigma=%.2f$)" % (mean, stddev)
    normal = mlab.normpdf(bins, mean, stddev)
    ax.plot(bins, normal, 'k--', label=text)

    ax.set_title("hESC/IMR90 Gene Expression Changes")
    ax.set_xlabel('Expression Log2 Fold Change')
    ax.set_ylabel('Number of Genes')

    fname = nu.join(histDir, 'both.png')
    print("Saving to", fname)
    fig.savefig(fname, dpi=600)
    plt.close()
    return

def probeChangesByCompartmentChanges(n=0, res="200k"):
    """Plots a volcano plot of probe changes by compartment change"""
    data = loadRawArray()
    icomponents = np.loadtxt(tmpl.format("IMR90", "R1", res))
    hcomponents = np.loadtxt(tmpl.format("hESC", "R1", res))
    ipc = icomponents[n] # component n
    hpc = hcomponents[n]

    # calculate the change in principal components, biasing for large changes
    # this is the "compartment change"
    pc = np.abs(ipc - hpc)
    
    # create columns for quick iteration
    iCols = ["imr90_{0}".format(i) for i in xrange(1,3)]
    hCols = ["hesc_{0}".format(i) for i in xrange(1,3)]
    iMeans = np.array([np.mean([np.log(i),np.log(j)]) for i,j in data[iCols]])
    hMeans= np.array([np.mean([np.log(i),np.log(j)]) for i,j in data[hCols]])
    delta = iMeans - hMeans

    # get the numerical resolution
    nres = resMap[res]
    starts = data['start'] / nres
    ends = (data['end'] / nres) + 1

    compartments = np.array([np.mean(pc[start:end]) for start, end in zip(starts, ends)])
   
    fig, ax = plt.subplots() 

    # plot compartment change versus expression change
    ax.plot(delta, compartments, '.k', alpha=0.5, markersize=0.5)
    #ax.set_ylim(0, np.max(delta)+1)
    #ax.set_xticks([])
    ax.set_ylabel("Compartment Change")
    ax.set_xlabel("Expression Change")

    fname = nu.join(volcanoDir, "volcano-pc{0}.png".format(n+1))
    print("Saving file to", fname)
    fig.savefig(fname, dpi=600)
    plt.close()
    return

def plotAll():
    for cellType in nu.datasets.keys():
        cellTypeExpressionHistogram(cellType)
        for rep in nu.datasets[cellType]:
            for res in ["200k", "1000k"]:
                probesByCompartment(cellType, rep, res)

    expressionChangeHistogram()
    return

if __name__ == "__main__":
    #plotAll()
    for i in range(3):
        probeChangesByCompartmentChanges(n=i)
