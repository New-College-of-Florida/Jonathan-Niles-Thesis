#!/usr/bin/env python
"""
script.py
Calculates the distances between genetic lesions and domains

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
from itertools import cycle
from pybedtools import BedTool
import matplotlib.pyplot as plt


inDir = nu.join(nu.sync, "data/domains/bed/")
tmpl = "{0}-{1}-{2}kb.bed"
windows = [100, 200, 400, 800, 1000]

def mutationTypePieChart():
    lesions = BedTool('tcga.bed')
    types = [l[-1].strip() for l in lesions]
    utypes = np.unique(types)
    counts = [types.count(utype) for utype in utypes]

    fig, ax = plt.subplots(figsize=[10,10])

    cmap = plt.cm.prism
    colors = cmap(np.linspace(0., 1., len(counts)))

    patches = ax.pie(counts, labels=utypes, colors=colors, labeldistance=1.05)
    for wedge in patches[0]:
        wedge.set_edgecolor('white')

    ax.set_title("Mutation Types");

    ax.legend(patches, utypes, loc='left center', bbox_to_anchor=(-0.1, 1.),
                       fontsize=8)

    fname = nu.join(nu.sync, "plots/genes/mutationTypePieChart.png")
    print("Saving mutation types to", fname)
    plt.savefig(fname, bbox_inches='tight', dpi=600)
    return

def distanceToLesions(cellType="IMR90", rep="R1"):
    lesions = BedTool('tcga.bed')

    data = []
    for win in windows:
        fname = nu.join(inDir, tmpl.format(cellType, rep, win))
        print("Loading data from", fname)
        domains = BedTool(fname)
        
        # mutations in the domain
        lesionsInDomains = lesions.intersect(domains)

        # compare to random
        shuffles = []
        for i in xrange(250):
            rand = domains.shuffle(genome='hg19', chrom=True)
            shuffles.append(len(lesions.intersect(rand)))
        randCounts = np.mean(shuffles)
        data.append((len(lesionsInDomains), randCounts))

    print(data)
    fname = "{0}-{1}-shuffle.npytxt".format(cellType,rep)
    print("Saving data to ", fname)
    np.savetxt(fname, data)
    #ax, fig = plt.subplots()
    return

def distanceToMutations(cellType="IMR90", rep="R1"):
    """"""
    lesions = BedTool('tcga.bed')

    data = []
    for win in windows:
        fname = nu.join(inDir, tmpl.format(cellType, rep, win))
        print("Loading data from", fname)
        domains = BedTool(fname)
        
        # mutations in the domain
        lesionsInDomains = lesions.intersect(domains)

        # compare to random
        distribution = []
        for i in xrange(250):
            rand = domains.shuffle(genome='hg19', chrom=True)
            distribution.append(len(lesions.intersect(rand)))

        # add in data
        data.append((len(lesionsInDomains), distribution))

    np.savetxt("histogram.npytxt", np.array(data))
    return data

if __name__ == "__main__":
    mutationTypePieChart()
    #distanceToLesions() 
