#!/usr/bin/env python
"""
Exports edges from the half-domains constructed earlier
author : jniles
date : March 18th, 2015
"""

import os
import csv
import itertools
import numpy as np
import nutils as nu
from matplotlib.mlab import csv2rec
from collections import namedtuple

# templates and paths
bedDir = os.path.join(nu.sync, "data/domains/beds/")
outDir = nu.chkdir(os.path.join(nu.sync, "data/domains/edges/"))
tmpl = "{0}-{1}-{2}kb.bed"
columns = ['chrom', 'start', 'end', 'value', 'sign', 'length', 'id']

Peak = namedtuple('Peak', ['chrom', 'start', 'end', 'id'])

def parseBed(cellType, replicate, windowSize):
    """
    parses a bed file given in the template
    """
    fname = os.path.join(bedDir, tmpl.format(cellType, replicate, windowSize))
    return csv2rec(fname, delimiter='\t', names=columns)

def writeBed(fname, data):
    """
    Writes a tab-delimited .bed file using the fname and data passed to it
    """
    print "writing domains to %s" % fname
    writer = csv.writer(open(fname, 'w'), delimiter='\t')
    writer.writerows([list(dom) for dom in data])
    return

def findEdges(cellType, replicate, window):
    """
    Reads domains from .bed files and writes the borders/boundaries to a new
    .bed file.

    :param cellType: cell line neame
    :param replicate: replicate label
    :param window: integer window size
    :returns: None
    """
    data = parseBed(cellType, replicate, window)
   
    peaks = []
    for dom in data:
        p1 = Peak(chrom=dom['chrom'], start=dom['start'], end=dom['start']+1, id=dom['id'])
        p2 = Peak(chrom=dom['chrom'], start=dom['end'], end=dom['end']+1, id=dom['id'])
        peaks.extend([p1 p2])
   
    fname = os.path.join(outDir, tmpl.format(cellType, replicate, window))
    writeBed(fname, peaks)
    return

def findAllEdges():
    """runs find edges for all celltypes/replicates"""
    for ct in nu.datasets.keys():
        for rep in nu.datasets[ct]:
            for size in [10, 20, 40, 80, 100]:
                findEdges(ct, rep, size)

if __name__ == "__main__":
    findAllEdges()
