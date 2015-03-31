#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import nutils as nu
from pybedtools import BedTool 
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2

path = "{0}-{1}-{2}kb.bed"
windows = [100, 200, 400, 800, 1000]
dataDir = nu.join(nu.sync, "data/domains/bed/")
outDir = nu.chkdir(nu.join(nu.sync, "plots/domains/venn/"))
genome = "/home/jniles/data/dna/hg19/"

def twoWaySizeVenn(A=("IMR90", "R1", 100), B=("IMR90", "R2", 200)):
    """Think about domains of different sizes intersecting each other"""

    a = BedTool(nu.join(dataDir, path.format(A[0], A[1], A[2])))
    b = BedTool(nu.join(dataDir, path.format(B[0], B[1], B[2])))

    left = (a - b).count()
    right = (b - a).count()
    center = (a + b).count()

    alabel = "{0} {1} {2}kb".format(A[0], A[1], A[2])
    blabel = "{0} {1} {2}kb".format(B[0], B[1], B[2])

    venn2(subsets=(left, right, center), set_labels=[alabel, blabel])
    plt.title("Domain Size Retention")
    fname = "scaling-{0}.{1}.{2}-vs-{3}.{4}.{5}.png".format(A[0], A[1], A[2], B[0], B[1], B[2])
    outFig = nu.join(outDir, fname)
    print("Saving figure to", outFig)
    plt.savefig(outFig, dpi=600)
    plt.close()
    return

if __name__ == "__main__":
    for cellType in nu.datasets.keys():
        for i in windows[:-1]:
            for j in windows[1:]:
                if i != j:
                    A = (cellType, "R1", i)
                    B = (cellType, "R1", j)
                    twoWaySizeVenn(A=A, B=B)
