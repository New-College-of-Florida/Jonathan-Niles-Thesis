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
import scipy.stats as stats
from matplotlib import mlab
from pybedtools import BedTool
import matplotlib.pyplot as plt

bedDir = nu.join(nu.sync, "data/domains/bed/")
conservedDir = nu.join(nu.sync, "data/domains/conserved/")
boundDir = nu.join(nu.sync, "data/domains/boundaries/")
tmpl = "{0}-{1}-{2}kb.bed"
windows = [100, 200, 400, 800, 1000]
genome = "/home/jniles/data/dna/hg19/"

def mutationTypePieChart():
    lesions = BedTool('tcga.bed')
    types = [l[-1].strip() for l in lesions]
    utypes = np.unique(types)
    counts = [types.count(utype) for utype in utypes]

    fig, ax = plt.subplots(figsize=[10,10])

    cmap = plt.cm.brg
    colors = cmap(np.linspace(0., 1., len(counts)))

    patches, texts = ax.pie(counts, colors=colors, labeldistance=1.05, startangle=90)
    for wedge in patches:
        wedge.set_edgecolor('white')

    ax.set_title("Mutation Types");
    ax.legend(patches, utypes, loc='best', fontsize=14)
    plt.axis('equal')
    plt.tight_layout()

    fname = nu.join(nu.sync, "plots/genes/mutationTypePieChart.png")
    print("Saving mutation types to", fname)
    plt.savefig(fname, bbox_inches='tight', dpi=500)
    plt.close()
    return

def mutationsInDomains(cellType="IMR90", rep="R1", conserved=False):
    """distance between domains and cancerous lesions"""
    lesions = BedTool('tcga.bed')
    saveDir = nu.chkdir(nu.join(nu.sync, "plots/domains/shuffled/"))

    data = []
    for win in windows:
        if conserved:
            fname = nu.join(conservedDir,"{0}.{1}kb.conserved.bed".format(cellType, win))
        else:
            fname = nu.join(bedDir, tmpl.format(cellType, rep, win))

        print("Loading data from", fname)
        domains = BedTool(fname)

        # mutations in the domain
        lesionsInDomains = lesions.intersect(domains)

        # compare to random
        distribution = []
        for i in xrange(500):
            rand = domains.shuffle(genome='hg19', chrom=True)
            distribution.append(len(lesions.intersect(rand)))

        # plot distribution
        fig, ax = plt.subplots()

        if conserved:
            plt.suptitle("{0}".format(cellType))
        else:
            plt.suptitle("{0} {1}".format(cellType, rep))

        ax.set_title("Resampled Null Distribution {0}kb".format(win))
        ax.set_xlabel("Number of Lesions in Domains")
        ax.set_ylabel("Frequency")

        print("Plotting histogram for", win, "kb")
        n, bins, patches = ax.hist(distribution, bins=25, histtype="step",
                color="b", normed=True, label="Shuffled Domains")

        # calculating normal distribution, p value
        mean, stddev = np.mean(distribution), np.std(distribution)
        ndist = mlab.normpdf(bins, mean, stddev)
        text = "N($\mu={0:0.2f}$, $\sigma={1:0.2f}$)".format(mean, stddev)
        ax.plot(bins, ndist, "r--", label=text)

        # plot the actual mean
        v = len(lesionsInDomains)
        ax.axvline(v, color="k", label="Observed Number")
        pValue = 1 - stats.norm(mean, stddev).cdf(v)
        pText = "$p$-value = {0:0.03f}".format(pValue)

        m = np.max(ax.get_ylim())
        ax.set_ylim(0, m+(m*0.1))
        ax.legend(loc='best', title=pText, fancybox=True)
        ax.minorticks_on()

        # save stuff
        if conserved:
            fname = nu.join(saveDir, "window-conserved-{0}kb.png".format(win))
        else:
            fname = nu.join(saveDir, "window-{0}kb.png".format(win))

        print("Saving to", fname)
        fig.savefig(fname, dpi=500)
        plt.close()
    return

def mutationsAtBoundaries(cellType="IMR90", size=5000, iterations=25):
    """resampling the regions around domain boundaries"""
    saveDir = nu.chkdir(nu.join(nu.sync, "plots/domains/shuffled/"))
    lesions = BedTool('tcga.sorted.bed')

    for win in windows:
        fname = nu.join(boundDir, "{0}-{1}kb.unique.boundaries.bed".format(cellType, win))

        # create slop region
        bounds = BedTool(fname).slop(g=genome, b=size)

        overlap = len(bounds.intersect(lesions))

        # resampling
        print("Resampling", iterations, "times for window", win)
        results = bounds.randomintersection(lesions, iterations=iterations,
                shuffle_kwargs={'chrom': True, 'genome' : "hg19"}, debug=False)

        # convert generator to list
        distribution = list(results)
        print("Finished resampling. Plotting...")

        # plot distribution
        fig, ax = plt.subplots()
        plt.suptitle("{0}".format(cellType))

        ax.set_title("Resampling {0}kb Window Sizes".format(win))
        ax.set_xlabel("Number of Lesions in {0}kb of Boundaries".format(size/1000))
        ax.set_ylabel("Frequency")

        print("Plotting histogram for", win, "kb")
        n, bins, patches = ax.hist(distribution, bins=50, histtype="stepfilled",
                color="g", normed=True, label="Shuffled Domains")
        plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

        # calculating normal distribution, p value
        mean, stddev = np.mean(distribution), np.std(distribution)
        ndist = mlab.normpdf(bins, mean, stddev)
        text = "N($\mu={0:0.2f}$, $\sigma={1:0.2f}$)".format(mean, stddev)
        ax.plot(bins, ndist, "r--", label=text)

        # plot the actual mean
        ax.axvline(overlap, color="k", label="Observed Number")
        pValue = 1 - stats.norm(mean, stddev).cdf(overlap)
        pText = "$p$-value = {0:0.03f}".format(pValue)

        # formatting
        m = np.max(ndist)
        ax.set_ylim(0, m + (m*0.2))
        ax.legend(loc='best', title=pText, fancybox=True)
        ax.minorticks_on()

        # save
        fname = nu.join(saveDir, "{0}.boundaries.{1}kb.windows.{2}kb.slop.png".format(cellType, win, size))
        print("Saving to", fname)
        fig.savefig(fname, dpi=500)
        plt.close()

    return

if __name__ == "__main__":
    #mutationTypePieChart()
    #mutationsInDomains(conserved=True)
    mutationsAtBoundaries(size=5000, iterations=1000)
