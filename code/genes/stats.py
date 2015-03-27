#!/usr/bin/env python
"""
stats.py
kolmogorov-smirnov test to see if the samples came from the same distribution

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
from matplotlib.mlab import csv2rec
from scipy.stats import spearmanr, ks_2samp as ks

def loadRawArray():
    """loads raw expression array"""
    fname = nu.join(nu.sync, 'data/genes/expression.bed.gz')
    return csv2rec(fname, delimiter='\t')

def computeKSStatistics(cellType):
    """compute the KS test"""
    
    data = loadRawArray()

    columns = ["{0}_{1}".format(cellType.lower(), i) for i in xrange(1,3)]

    rep1 = np.array([np.log2(i) for i in data[columns[0]]])
    rep2 = np.array([np.log2(i) for i in data[columns[1]]])

    N = 250
    setA = np.histogram(rep1, bins=N)
    setB = np.histogram(rep2, bins=N)

    # run Kolmogorov-Smirnov test
    result = ks(setA, setB)
    print("Results of K-S test for {0}".format(cellType), result)
    return

def runAllTests():
    for cellType in nu.datasets.keys():
        computeKSStatistics(cellType)

if __name__ == "__main__":
    runAllTests()
    
