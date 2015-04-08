#!/usr/bin/env python
"""
plot.py
plots domains for a cell type

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

domBedDir = nu.join(nu.sync, "data/domains/bed/")
domPltDir = nu.chkdir(nu.join(nu.sync, "plots/domains/bar/"))

# templates and paths
tmpl = "{0}-{1}-{2}kb.bed"
columns = ['chrom', 'start', 'end', 'value', 'sign', 'length', 'id']
windows = [100, 200, 400, 800, 1000]

def parseBed(cellType, replicate, windowSize):
    """
    parses a bed file given in the template
    """
    fname = nu.join(domBedDir, tmpl.format(cellType, replicate, windowSize))
    return csv2rec(fname, delimiter='\t', names=columns)

def group(datasets):
    """
    use itertools.groupby to group datasets by chromosome
    """
    groups = defaultdict(list)
    for ds in datasets:
        ds.sort() # must be sorted for itertools
        for k, g in itertools.groupby(ds, lambda x: x['chrom']):
            groups[k].append(len(list(g)))

    # make sure to convert chromosomes
    groups['chr23'] = groups['chrX']
    groups.pop('chrX', None)
    return groups

def getMeans(cellType, window):
    """
    """
    datasets = [parseBed(cellType, rep, window) for rep in nu.datasets[cellType]]
    groups = group(datasets)

    # sort the (text-labeled) chromosomes
    keys = sorted(groups.keys(), cmp=lambda x,y: cmp(int(x[3:]), int(y[3:]))) 

    n = len(datasets)
    means = map(lambda k: np.mean(groups[k]), keys)
    stderr = map(lambda k: np.std(groups[k]) / np.sqrt(n), keys)
    return (means, stderr , keys)

def autolabel(axis, rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        axis.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                ha='center', va='bottom', size='5')

def domainbar(cellType, labelBar=False):
    """
    Draws a bar graph of the average domain size between technical replicates
    with error bars.
    """
    sizes = []
    for win in windows:
        sizes.append(getMeans(cellType, win))

    # width of the bars and number of groups
    width = 0.15
    ind = np.arange(23)

    fig, ax = plt.subplots()
    ax.set_title("{0} Domains by Chromosome".format(cellType))
    ax.set_xlabel("Chromosomes")
    ax.set_ylabel("Domains")

    # get colors from colormap
    e = plt.cm.gist_earth
    colors = map(e, np.linspace(0,0.9, len(sizes)))

    # this is the most pythonic-code ever
    for i in xrange(len(sizes)):
        means, stderr, keys = sizes[i]
        rects = ax.bar(ind+(i*width), means, width, yerr=stderr,
                color=colors[i], label="{0}kb".format(windows[i]))
        if labelBar:
            autolabel(ax, rects)

    xlabels = map(lambda x: "chr{0}".format(x) if x != 23 else "chrX", ind+1)

    ax.legend(loc="best")
    ax.set_xticks(ind + (len(sizes)*width / 2.))
    ax.set_xticklabels(xlabels, rotation=45)
    plt.tight_layout()
    
    # save figure
    fname = nu.join(domPltDir, "{0}-bar.png".format(cellType))
    print("Saving bar chart to", fname)
    fig.savefig(fname, dpi=500)
    plt.close()
    return

def domainSizeByWindowSize():
    """plots domain sizes as they scale for each cell type/replicate
    by the window size"""

 
    data = np.arange(1,len(windows)+1)
    for j, win in enumerate(windows):
        sizes = []
        for cellType in nu.datasets.keys():
            for rep in nu.datasets[cellType]:
                domains = BedTool(nu.join(domBedDir, tmpl.format(cellType, rep, win)))
                sizes.append(np.mean(map(len, domains)) / 1000.0)
        data[j] = np.mean(sizes)

    fig, ax = plt.subplots()

    ax.plot(np.arange(1, len(windows)+1), data, label="Domain Size")
    ax.set_xlim(0, len(windows) + 1)
    ax.grid()
    
    ax.set_title("Domain Size Scaling vs Window Size")
    wins =  windows[:]
    wins.insert(0,0)
    ax.set_xticklabels(map(lambda i: "{0}kb".format(i), wins))
    ax.set_xlabel("Window Size")
    ax.set_ylabel("Domain Size (kb)")

    fname = nu.join(nu.sync, "plots/domains/sizes.png")
    print("Saving to", fname)
    fig.savefig(fname)
    plt.close()
    return

def plotAllBars():
    for cellType in nu.datasets.keys():
        domainbar(cellType)

if __name__ == "__main__":
    #domainSizeByWindowSize()
    plotAllBars()
