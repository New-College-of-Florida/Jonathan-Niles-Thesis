#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import nutils as nu
from pybedtools import BedTool 
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2

path = "{0}-{1}-{2}kb.bed"
windows = [100, 200, 600, 800, 1000]
dataDir = nu.join(nu.sync, "data/domains/bed/")
outDir = nu.chkdir(nu.join(nu.sync, "plots/domains/venn/"))
genome = "/home/jniles/data/dna/hg19/"

# make a 2-way venn diagram of domains
def threeWayVenn(reps=["R1", "R2", "R3"], cellType="IMR90", win="100"):
    """plots the domain preservation between the first three replicates
    of the IMR90 cell line"""

    mkpath = lambda rep: nu.join(dataDir, path.format(cellType, rep, win))
    paths = map(mkpath, reps)
    a,b,c = map(BedTool, paths)

    # calculate venn regions
    au = (a-b-c).count()
    bu = (b-a-c).count()
    cu = (c-b-a).count()
    ab = (a+b-c).count()
    ac = (a+c-b).count()
    bc = (b+c-a).count()
    abc = (c+a+b).count()

    vlabel = lambda rep : "{0} {1}".format(cellType, rep)
    labels = map(vlabel, reps)
    
    venn3({
        '100': au,
        '010': bu,
        '001': cu,
        '110': ab,
        '011': bc,
        '101': ac,
        '111':abc
    }, set_labels=labels)

    plt.title("Domain Discovery by Replicate")
    outFig = nu.join(outDir, "venn3-{0}kb.png".format(win))
    print("Saving to", outFig)
    plt.savefig(outFig, dpi=600)
    plt.close()
    return

"""
Do not worry about this code.  It produces 0 as a result!
"""
def twoWayRandomVenn(rep="R1", cellType="IMR90", win="100"):
    """plots the domain preservation between the first three replicates
    of the IMR90 cell line"""

    domains = BedTool(nu.join(dataDir, path.format(cellType, rep, win)))
    centers = domains.randomintersection(domains, iterations=1000,
            shuffle_kwargs={"chrom" : True, "g" : genome},
            debug=True)

    center = np.mean([v fo v in centers])
    sides = len(domains) - center

    dlabel = "{0} {1}".format(cellType, rep)
    rlabel = "Shuffled Domains"
    
    venn2(subsets=(sides, sides, center), set_labels=[dlabel, rlabel])

    plt.title("Shuffled Domain Overlaps")
    outFig = nu.join(outDir, "venn2shuffled-{0}kb.png".format(win))
    print("Saving to", outFig)
    plt.savefig(outFig, dpi=600)
    plt.close()
    return

if __name__ == "__main__":
    threeWayVenn()
    #twoWayRandomVenn()
