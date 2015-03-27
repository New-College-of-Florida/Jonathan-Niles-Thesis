#!/usr/bin/env python
"""
script.py
Calculates the directionality index for a given (high-res) chromosome
as described in Dixon et al.  Writes the result to a .npy array.

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

import numpy as np
import nutils as nu
from matplotlib import pyplot as plt
from hiclib.highResBinnedData import HiResHiC

resolution = 10000
genome = "/home/jniles/data/dna/hg19/"
tmpl = "/home/jniles/data/{0}/ic/{0}-{1}-HindIII-10k_HighRes.byChr"
windows = [5, 10, 20, 30, 40, 50]

Loaders = {}

# ensure that the data directory exists
outDir = nu.chkdir(nu.join(nu.sync,"data/di/"))

def di(subarray):
    """calculates the directionality of a particular interaction window

    :param subarray: a numpy array
    :returns index: a +/- float

    S  -- window size
    A  -- sum of the upstream interactions
    B  -- sum of the downstream interactions
    E  -- null hypothesis
    """
    s = subarray.shape[0]
    half = s / 2
    A = subarray[:half, :half].sum() # downstream
    B = subarray[half:,half:].sum() # upstream
    E = (A + B) / 2
    if (B - A) == 0:  # avoid division by zero
        return 0
    return (float((B - A)) / np.abs(B - A))*(pow(A-E,2) / E + pow(B-E,2) / E)

def runDI(data, window=10):
    """Calculates the DI for a region"""
    # for each position along the diagonal, calculate the DI
    # with bins of 1MB
    N = data.shape[0]
    directions = []
    for i in xrange(0, N):
        # make sure that indices are bounded by the array
        up = i + window
        if up > N:
            up = N
        down = i - window
        if down < 0:
            down = 0

        subarray = data[down:up,down:up]
        directions.append(di(subarray))
    return np.array(directions)

def export(cellType, rep, array, windowSize):
    """writes the directionality index to the sync directory"""
    res = windowSize * 2 * 10000
    kb = 1000
    suffix = "{0}kb".format(res / kb)

    fname = nu.join(outDir, "{0}-{1}-{2}.npy".format(cellType, rep, suffix))
    print("Writing index to ", fname)
    np.save(open(fname, "w"), array)
    return

def getLoader(cellType, rep):
    """Memoize the loader so that we don't do lots of i/o"""
    print("Retrieving data for {0} {1}".format(cellType, rep))
    global Loaders
    if not Loaders.get((cellType, rep)):
        loader = HiResHiC(genome, resolution)
        loader.loadData(tmpl.format(cellType, rep), mode="cis")
        Loaders[(cellType, rep)] = loader
    return Loaders[(cellType, rep)]

def calculateDI(cellType, rep, window):
    """calculates the directionality for the entire genome"""

    print("Calculating DI for {0} {1}".format(cellType, rep))
    print("Window size: {0}".format(window))

    genomeIdx = []
    loader = getLoader(cellType, rep)
    for chrom in xrange(0, 23):
        chromIdx = (chrom, chrom)
        data = loader.data[chromIdx].getData()

        # for each position along the diagonal, calculate the DI
        # with bins of 1MB
        idx = runDI(data, window=window)
        genomeIdx.append(idx)

    array = np.array(genomeIdx)
    export(cellType, rep, array, window)
    return

def calculateAll():
    """runs DI for all cell lines"""
    for win in windows:
        for cellType in nu.datasets.keys():
            for rep in nu.datasets[cellType]:
                calculateDI(cellType, rep, win)

if __name__ == "__main__":
    calculateAll()
