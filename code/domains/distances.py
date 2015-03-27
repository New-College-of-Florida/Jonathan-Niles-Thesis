#!/usr/bin/env python
"""
distances.py
Calculates the distances between replicates compared to the distances
between n

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
import matplotlib.pyplot as plt

inDir = nu.join(nu.sync, "data/domains/bed/")
outDir = nu.chkdir(nu.join(nu.sync, "plots/domains/closest/"))
tmpl = "{0}-{1}-{2}kb.bed"
windows = [100, 200, 400, 800, 1000]

def readBed(cellType, rep, win):
    """read in a bed"""
    fname = nu.join(inDir, tmpl.format(cellType, rep, win))
    return BedTool(fname)

def closest(cellType, win, io=True):
    """Calculates the closest peaks between the same replicate
    and all other replicates in the cellType"""
    
    print("Calculating closest domains for", cellType, win)
    reps = map(lambda r: readBed(cellType, r, win), nu.datasets[cellType])
    self = reps[0]
    nearestReplicates = [self.closest(rep, io=io, d=True) for rep in reps[1:]]
    nearestReplicates.insert(0, self.closest(self, io=io, d=True))
    repMeans = [np.mean([int(line[-1]) for line in rep]) for rep in nearestReplicates]

    x = np.arange(1,len(repMeans)+1)

    fig, ax = plt.subplots()
    ax.plot(x, repMeans, '-o', color="k", label="Avg. Distances")
    ax.set_xlabel("Replicate")
    ax.set_ylabel("Distance (bp)")
    ax.set_title("Distance to Nearest Domain From {0}".format(nu.datasets[cellType][0]))

    plt.suptitle("{0} {1}kb".format(cellType, win))
    plt.ylim(0, np.max(repMeans)+0.1*np.max(repMeans))
    plt.xlim(0, len(x)+1)

    ax.set_xticks(x)
    ax.set_xticklabels(nu.datasets[cellType])

    if io:
        fname = nu.join(outDir, "{0}-{1}kb.png".format(cellType, win))
    else:
        fname = nu.join(outDir, "{0}-{1}kb-overlaps.png".format(cellType, win))
    print("Saving figure to", fname)
    plt.savefig(fname)
    plt.close()
    return

if __name__ == "__main__":
    for cellType in nu.datasets:
        for win in windows:
            closest(cellType, win, io=True)
