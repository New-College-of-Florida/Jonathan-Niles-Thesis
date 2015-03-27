#!/usr/bin/env python
"""
export.py
export the directionality index to a .bed file at 10kb resolution

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

import csv
import itertools
import numpy as np
import nutils as nu
from collections import namedtuple

Domain = namedtuple('Domain', ['chrom', 'start', 'end', 'value', 'sign', 'length', 'id'])

resolution = 10000
template = "{0}-{1}-{2}kb.npy"
diDir = nu.chkdir(nu.join(nu.sync, "data/di/"))
bedDir = nu.chkdir(nu.join(nu.sync, "data/domains/bed/"))
windows =  [100, 200, 400, 800, 1000]

def loadAll(cellType, replicate, resolution):
    """loads the directionality index into a numpy array"""
    fname = nu.join(diDir, template.format(cellType, replicate, resolution))
    print("Loading data from", fname)
    return np.load(fname)

def domains2bed(domains, fname):
    """
    export the given file to a .bed file at 10kb res
    """
    print("Writing %i domains to %s" % (len(domains), fname))

    with open(fname, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for domain in domains:
            writer.writerow(list(domain))
    return

def findRegions(data):
    """
    Finds where the directionality index crosses the zero axis
    If zeros = True, will return individual indexes rather than tuples of
    grouped indexes
    """
    diffs, = np.where(np.diff(np.sign(data))) # detect sign changes
    regions = collapseConsecutive(diffs)
    return regions

def collapseConsecutive(z):
    """
    DEPRECATED
    collapses consecutive numbers in to their start and end positions.
    """
    # duplicate and shift
    x = np.zeros(z.shape)
    x[1:] = z[:-1]
    arr = []
    for i, j in zip(x,z):
        if j - i != 1:
            arr.append(j)
    return np.array(arr)

def identifyDomains(data, boundaries, chrom=0):
    """
    Takes in an array of positions where the line crosses and calls domains
    based on those positions.
    """
    # pad the boundaries with a zero start, and the max length
    offset = np.zeros(len(boundaries) + 2, dtype=np.int)
    offset[1:-1] = boundaries
    offset[-1] = len(data) - 1
    
    # the entire domain must sum to greater than the mean value
    maxd = np.max(np.abs(data))
    #m = np.mean(np.abs(data))
    m = 0.9*maxd

    if chrom == 23:
        chrom = "X"
        
    domains = zip(offset, offset[1:])
    discovered = []
    i = 0
    for start, end in domains:
        v = np.sum(np.abs(data[start:end]))
        if v > m: # we will keep this!  Passes domain test
            sign = np.sign(np.sum(data[start:end])) # with any luck, this will not be zero

            dstart = start * resolution
            dend = (end * resolution)
            length = (end - start) * resolution
            chrm = 'chr{0}'.format(chrom)
            i += 1
            id = "Dom{0}_Chr{1}".format(i, chrom)
            dom = Domain(chrom=chrm, start=dstart, end=dend, value=v, sign=sign, length=length, id=id)
            discovered.append(dom)
    return discovered

def concatDirectionality(cellType, rep, res):
    """
    Takes files from the indices/ folder and makes on large matrix 
    out of them, writing them back to the domains/
    """
    print("Compiling files from %s %s %s" % (cellType, rep, res))
    inFile = 'indices/{0}-{1}-chr{2}-{3}bp.npy'
    data = np.asarray([np.load(inFile.format(cellType, rep, i, res)) for i in xrange(23)])
    pRes = res / resolution
    fname = "byReplicate/{0}-{1}-{2}kb.npy".format(cellType, rep, pRes)
    np.save(fname, data)
    return

def exportBedFromDI(cellType, rep, res):
    """
    Runs the module
    """
    domains = []
    data = loadAll(cellType, rep, res)
    for i, indices in enumerate(data):
        regions = findRegions(indices)
        chromDomains = identifyDomains(indices, regions, chrom=i+1)
        domains.extend(chromDomains)

    savePath = nu.join(bedDir, '{0}-{1}-{2}kb.bed'.format(cellType, rep, res))
    domains2bed(domains, savePath)
    return

def exportAll():
    """computes di for all datasets"""
    for size in windows:
        for cellType in nu.datasets.keys():
            for rep in nu.datasets[cellType]:
                exportBedFromDI(cellType, rep, size)

if __name__ == "__main__":
    exportAll()
