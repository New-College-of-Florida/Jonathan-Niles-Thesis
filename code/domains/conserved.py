#!/usr/bin/env python
"""
conserved.py
takes intersections of all replicates to find conserved domains.

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
from pybedtools import BedTool

dataDir = nu.join(nu.sync, "data/domains/bed/")
exportDir = nu.chkdir(nu.join(nu.sync, "data/domains/conserved/"))

# templates and paths
tmpl = "{0}-{1}-{2}kb.bed"
windows = [100, 200, 400, 800, 1000]

def loadBed(cellType, rep, window):
    """loads a bed file in for work"""
    fname = nu.join(dataDir, tmpl.format(cellType, rep, window))
    print("Loading .bed record", fname)
    return BedTool(fname)

def filterConserved(cellType, window=200):
    """this will do things"""
    beds = [loadBed(cellType, rep, window) for rep in nu.datasets[cellType]]
    f = lambda x,y : x+y
    conserved = reduce(f, beds)
    print("Found", len(conserved), "conserved domains.")
    fname = nu.join(exportDir, "{0}.{1}kb.conserved.bed".format(cellType, window))
    print("Saving as", fname)
    conserved.saveas(fname) 
    return

def findAllConserved():
    for cellType in nu.datasets:
        for win in windows:
            filterConserved(cellType, window=win)

if __name__ == "__main__":
    findAllConserved()
