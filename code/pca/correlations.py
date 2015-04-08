#!/usr/bin/env python
"""
correlations.py
calculates correlations for the principal components

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
from scipy.stats import pearsonr, spearmanr

# template for loading data
tmpl = nu.join(nu.sync, "data/components/{0}-{1}-{2}.eigenvectors")
outDir = nu.chkdir(nu.join(nu.sync, "data/components/correlations/"))

def computeCorrelations(setA, setB, numComponents):
    """
    computes correlations between the `numComponents` rows of setA and setB
    """
    if numComponents == 1:
        return [spearmanr(setA[i], setB[i]) for i in xrange(numComponents)][0]
    else:
        return [spearmanr(setA[i], setB[i]) for i in xrange(numComponents)]

def betweenReplicates(cellType, resolution, numComponents=1, saveCSV=True):
    """
    computes the correlations for all the first principal components
    of a cell type.
    """
    data = []
    for rep in nu.datasets[cellType]:
        data.append(np.loadtxt(tmpl.format(cellType, rep, resolution)))

    l = len(data)
    mtx = np.ndarray((l,l), dtype=tuple)
    for i, sA in enumerate(data):
        for j, sB in enumerate(data):
            mtx[i,j] = computeCorrelations(sA, sB, numComponents)

    outFile = nu.join(outDir, "{0}-{1}-pc{2}.npy".format(cellType, resolution, numComponents))

    # save a matrix copy
    np.save(outFile, mtx)

    # save a csv copy
    if saveCSV:
        np.savetxt(outFile+"txt", np.asarray(mtx), delimiter=",", fmt="%s")
    return

def betweenCellTypes(resolution, numComponents=1, saveCSV=True):
    """
    computes the correlations between PCs of different cell types
    """
    i90 = [np.loadtxt(tmpl.format("IMR90", v, resolution)) for v in nu.datasets["IMR90"]]
    h1 = [np.loadtxt(tmpl.format("hESC", v, resolution)) for v in nu.datasets["hESC"]]

    mtx = np.ndarray((len(i90),len(h1)), dtype=tuple)
    for i, repI in enumerate(i90):
        for j, repH in enumerate(h1):
            mtx[i, j] = computeCorrelations(repI, repH, numComponents)

    outFile = nu.join(outDir, "cross-{0}-pc{1}.npy".format(resolution, numComponents))

    # save an array copy
    np.save(outFile, mtx)

    # save a csv copy
    if saveCSV:
        np.savetxt(outFile+"txt", np.asarray(mtx), delimiter=",", fmt="%s")
    return

def computeAllCellTypes():
    """runs betweenCellTypes() for each resolution computed"""
    for res in ["200k", "1000k"]:
        print("Computing", res)
        betweenCellTypes(res)
    return

def computeAllReplicates():
    """runs betweenReplicates() for all resolutions computed"""
    for res in ["200k", "1000k"]:
        for cellType in nu.datasets.keys():
            print("Computing ", cellType, res)
            betweenReplicates(cellType, res)
    return

if __name__ == "__main__":
    computeAllCellTypes()
    computeAllReplicates()
