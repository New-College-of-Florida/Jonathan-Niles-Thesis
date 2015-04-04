#!/usr/bin/env python
"""
script.py
plots a histogram of the compartment character in the fragile site

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

import numpy as np
import nutils as nu
import matplotlib.pyplot as plt 
from matplotlib.mlab import csv2rec
from hiclib.binnedData import binnedData

path = nu.join(nu.sync, "data/components/{0}-{1}-{2}.eigenvectors")
cfsPath = "/home/jniles/data/annotations/cfs.csv"
outDir = nu.chkdir(nu.join(nu.sync, "plots/fragile/"))
genome = "/home/jniles/data/dna/hg19/"

def loadComponent(cellType, rep, res):
    """loads the PC for a given celltype"""
    fname = path.format(cellType, rep, res)
    print("Loading component file", fname)
    return np.loadtxt(fname)

def plotCfsVsComponent(cfs, pc, res, spath):
    """plots the cfs versus components"""

    # convert resolution if we need to
    if type(res) == type("s"):
        res = nu.strToResolution(res)

    bd = binnedData(res, genome=genome)

    compartments = []
    for site in cfs:

        # strip off all but the numbers/X
        chromosome = site['chromosome'].strip()[3:]
        if chromosome == 'X':
            chrom = 23
        else:
            chrom = int(chromosome)
        
        sbin = site['start'] / res
        ebin = site['stop'] / res

        component = pc[bd.chromosomeIndex == chrom]
        compartments.append(np.sum(component[sbin:ebin]))

    fig, ax = plt.subplots()

    ax.hist(compartments, bins=25, normed=True)

    plt.title("Compartment Character in Common Fragile Sites")
    plt.xlabel("Compartment Character")
    plt.ylabel("Frequency")
    xl = np.max(np.abs(ax.get_xlim()))
    ax.set_xlim(-xl, xl)

    outFig = nu.join(outDir, "{0}-histogram.png".format(spath))
    print("Saving to", outFig)
    fig.savefig(outFig)
    plt.close()
    return

def main(res="200k", n=0):
    """plots the fragile eigenvectors along the fragile sites"""
    
    # load chromosome fragile sites
    cfs = csv2rec(cfsPath)

    # load each component
    for cellType in nu.datasets.keys():
        for rep in nu.datasets[cellType]:
            pc = loadComponent(cellType, rep, res)[n]
            spath = "{0}-{1}".format(cellType, rep)
            plotCfsVsComponent(cfs, pc, res, spath)

if __name__ == "__main__":
    main()
