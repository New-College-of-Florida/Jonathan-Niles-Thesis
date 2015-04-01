#!/usr/bin/env python

import numpy as np
import nutils as nu
from matplotlib.mlab import csv2rec

path = nu.join(nu.sync, "/data/components/{0}-{1}-{2}.eigenvectors")
cfsPath = "/home/jniles/data/annotations/cfs.csv"

def loadComponent(cellType, rep, res):
    """loads the PC for a given celltype"""
    fname = path.format(cellType, rep, res)
    print("Loading component file", fname)
    return np.loadtxt(fname)

def plotCfsVsComponent(cfs, pc):
    """plots the cfs versus components"""
    return

def main(res="200k"):
    """plots the fragile eigenvectors along the fragile sites"""
    
    # load chromosome fragile sites
    cfs = csv2rec(cfsPath)

    # load each component
    for cellType in nu.datasets.keys():
        for rep in nu.datasets[cellType]:
            pc = loadComponent(cellType, rep, res)
            plotCfsVsComponent(cfs, pc)



if __name__ == "__main__":
    main()

