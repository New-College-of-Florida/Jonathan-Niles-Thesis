#!/usr/bin/env python
"""
corr.py
Calculates the spearman's correlation for the directionality indices
generated by the previous scripts.

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
from scipy.stats import spearmanr

inDir = nu.join(nu.sync, "data/di/")
outDir = nu.chkdir(nu.join(nu.sync, "data/di/corr/"))
template = "{0}-{1}-{2}kb.npy"
windows = [100, 200, 400, 600, 800, 1000]

def loadAll(cellType, replicate, resolution):
    print("Loading data for", cellType, replicate, resolution)
    fname = template.format(cellType, replicate, resolution)
    return np.load(nu.join(inDir, fname))

def calculateGenomeCorrelation(window):
    """
    calculates spearman's correlation for the entire genome between
    all cell types for a given window size
    """

    genomes = []
    for key in nu.datasets.keys():
        for v in nu.datasets[key]:
            genome = list(itertools.chain(*loadAll(key, v, window)))
            genomes.append(genome)

    length = len(genomes)
    mtx = np.zeros((length, length))
    for i in xrange(length):
        for j in xrange(length):
            mtx[i,j] = spearmanr(genomes[i], genomes[j])[0]

    fname = nu.join(outDir, '{0}kb.window.npytxt'.format(window))
    np.savetxt(fname, mtx)
    return

def byChrom(chrom, window):
    """
    computes Spearmans' r by chromosome across all cellTypes
    """
    
    print("Calculating spearmans for chromosome %i %ikb" % (chrom+1, window))

    sets = []
    for key in nu.datasets.keys():
        for value in nu.datasets[key]:
            sets.append(loadAll(key, value, window)[chrom])

    length = len(sets)
    mtx = np.zeros((length, length))
    for i in xrange(length):
        for j in xrange(length):
            mtx[i,j] = spearmanr(sets[i], sets[j])[0]

    fname = nu.join(outDir,'chr{0}-{1}kb.npytxt'.format(chrom, window))
    np.savetxt(fname, mtx)
    return

def calculateAll():
    """runs the module"""
    for win in windows:
       calculateGenomeCorrelation(win)

if __name__ == "__main__":
    calculateAll()
