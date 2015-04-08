#!/usr/bin/env python
"""
com.py
This script

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

import pickle
import nutils as nu
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# constants
etype = "ic"
path = "/home/jniles/data/{0}/%s/{0}-{1}-HindIII-1000k.hm" % etype
dataPath = nu.chkdir(nu.join(nu.sync, "data/probes/"))
figPath = nu.chkdir(nu.join(nu.sync, "plots/probes/com/"))

def loadScalingData(cellType):
    """
    Loads scaling by replicate from preprocessed pickle data
    """
    inFile = nu.join(dataPath, '{0}-{1}.pkl'.format(etype, cellType))
    with open(inFile, 'r') as f:
        dataDict = pickle.load(f)
    return dataDict

def centerOfMass(array):
    total = np.sum(array)
    center = total / 2.
    pos = 0
    n = 0
    for i in xrange(len(array)):
        if n > center:
            pos = i
            break
        n += array[i]
    return pos

def plotComScatter():
    """plots a scatter plot of the center of mass for each dataset"""

    data = loadScalingData("hESC")

    # Generate some random colors
    colors = [np.random.rand(3,1) for x in xrange(23)]
    centers = defaultdict(list)

    for rep in data:

        # idx is the chrom number
        for idx, values in data[rep].items():
            # find the position of the center of mass:
            # sums all previous positions
            z = [sum(values[:i+1]) for i in xrange(len(values))]
            ind = np.where(z > (z[-1]/2))[0]
            
            # append the position (bin number) of the center of mass to
            # the centers dictionary
            centers[idx].append(ind[0])
    
    averages, errors = [], []
    for i in xrange(23):

        # retrieve the centers for each chromosome
        avg = np.mean(centers[i])
        err = np.std(centers[i]) / np.sqrt(len(centers[i]))
        
        # store in proper order
        averages.append(avg)
        errors.append(err)

    length = len(averages)
    labels = [str(i) for i in xrange(1,length)]
    labels.append('X')
   
    y = np.arange(length)

    plt.scatter(averages, y, c="k")
    ax = plt.gca()
    ax.set_yticks(range(length))
    ax.set_yticklabels(labels)

    plt.errorbar(averages, y, xerr=errors, linestyle="None", ecolor="k")
    plt.title('Half Max Probe Scaling')
    plt.ylabel('Chromosome')
    plt.xlabel('Distance')
    plt.grid(True)

    fname = nu.join(figPath, "scatter.png")
    print("Saving figure to", fname)
    plt.savefig(fname)
    plt.close()

if __name__ == "__main__":
    plotComScatter()
