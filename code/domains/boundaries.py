#!/usr/bin/env python
"""
boundaries.py
Gets the boundaries of domains

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
from pybedtools import BedTool
import matplotlib.pyplot as plt
from matplotlib.mlab import csv2rec
from collections import defaultdict

# templates and paths
conservedDir = nu.join(nu.sync, "data/domains/conserved/")
boundDir = nu.chkdir(nu.join(nu.sync, "data/domains/boundaries/"))
tmpl = "{0}.{1}kb.conserved.bed"
windows = [100, 200, 400, 800, 1000]

def loadBed(cellType, window):
    """loads a bed file in for work"""
    fname = nu.join(conservedDir, tmpl.format(cellType, window))
    print("Loading .bed record", fname)
    return BedTool(fname)

def getBoundaries(cellType, window):
    """Exports the boundaries from conserved domains"""
    conserved = loadBed(cellType, window)
    boundaries = []
    for domain in conserved:
        index = float(domain[4]) # directionality index
        if index < 1:
            field = 'end'
        else:
            field = 'start'
        pos = int(domain[field]) # get the bp position of boundary
        bound = (domain.chrom, pos, pos+1)
        boundaries.append(bound)
    beds = BedTool(boundaries)
    fname = nu.join(boundDir, "{0}-{1}kb.boundaries.bed".format(cellType, window))
    print("Saving boundaries to", fname)
    beds.saveas(fname)
    return

def getPeaks(cellType, window):
    """exports the peaks from conserved domains"""
    conserved = loadBed(cellType, window)
    boundaries = []
    for domain in conserved:
        index = float(domain[4]) # directionality index
        if index > 1:
            field = 'end'
        else:
            field = 'start'
        pos = int(domain[field]) # get the bp position of boundary
        bound = (domain.chrom, pos, pos+1)
        boundaries.append(bound)
    beds = BedTool(boundaries)
    fname = nu.join(boundDir, "{0}-{1}kb.peaks.bed".format(cellType, window))
    print("Saving boundaries to", fname)
    beds.saveas(fname)
    return

if __name__ == "__main__":
    for cellType in nu.datasets.keys():
        for window in windows:
            getBoundaries(cellType, window)
            getPeaks(cellType, window)
