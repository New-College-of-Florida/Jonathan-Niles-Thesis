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

import pickle
import nutils as nu
import numpy as np
import matplotlib.pyplot as plt

# constants
etype = "ic"
path = "/home/jniles/data/{0}/%s/{0}-{1}-HindIII-1000k.hm" % etype
dataPath = nu.chkdir(nu.join(nu.sync, "data/probes/"))
figPath = nu.chkdir(nu.join(nu.sync, "plots/probes/com/"))

def loadScalingData(cellType):
    """
    Loads scaling by replicate from preprocessed pickle data
    """
    datasets = nu.datasets[cellType]
    dataDict = {}
    for rep in datasets:
        inFile = nu.join(dataPath, '{0}-{1}-{2}.pkl'.format(etype, cellType, rep))
        with open(inFile, 'r') as f:
            dataDict[rep] = pickle.load(f)
    return dataDict

def plotCoM(data, rep, colors):
    """
    plot the center of mass of each chromosome in a replicate
    
    The "center of mass" is the distance were 50% of the interactions
    are less than that distance and 50% of the interactions are greater.
    """
    logger.stdout("[COM] Beginning analysis of replicate {0}".format(rep))
    
    maxLength = 0
    plt.subplot(111)
    
    # for each chromosome
    for k in data.keys():
        contacts = data[k]['contacts']
        length = data[k]['length']
      
        # in case I didn't calculate the total
        if 'total' in data[k].keys():
          total = data[k]['total']
        else:
          total = sum(contacts)
      
        # normalize to a percentage per bin
        norm = 100 * (contacts[1:] / total)
      
        # find the center of mass
        s, idx = 0, 0
        while s < 50:
            s += norm[idx]
            idx += 1
      
        centerOfMass = idx
        if centerOfMass > maxLength:
            maxLength = centerOfMass
      
        if k == 22:
            label = "chrX"
        else:
            label = "chr{0}".format(k+1)
      
        plt.scatter([centerOfMass], [0], s=35, c=colors[k], label=label)
        plt.plot(xrange(length), np.zeros(length), '-k', linestyle="dashed", linewidth=0.5)
      
    logger.stdout('Plotting figure.')
    
    plt.title('Replicate {0} Center of Mass'.format(rep))

    # remove ugly ticklables
    plt.setp(plt.gca().get_yticklabels(), visible=False)
    plt.ylim([-0.25, 1])
    plt.xlim([0, maxLength + maxLength/10.])
    plt.xlabel('Genomic Distance (Kb)')
    plt.legend(loc=9, ncol=4, mode="expand", borderaxespad=0., scatterpoints=1, scatteryoffsets=[0.5])
    plt.savefig('vis/corrected/com-{0}.png'.format(rep), dpi=300)
    
    # explicity close the figure to conserve memory
    plt.close()
    return

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
    errors = []
    avgs = []

    for chrom in data:
        # find the position of the center of mass:
        z = [sum(chrom[:i+1]) for i in xrange(len(chrom))]
        ind = np.where(z > (z[-1]/2))[0]

        avgs.append(np.mean(values))


    # for each chromosome, calculate the avg center of mass and
    # append it to the averages.  Also, calculate the standard
    # deviations and append it to errors
    for i in xrange(22,-1,-1):
        logger.stdout("Calculating center of mass for chromosome %i"%i)
        values = []
        for p in pickles:
            values.append(centerOfMass(p[i]['contacts'], total=p[i]['total']))
        errors.append(np.std(values))
        avgs.append(np.average(values))

    length = len(avgs)
    labels = [str(i) for i in xrange(1,length)]
    labels.append('X')
   
    y = np.arange(length)

    logger.stdout('Drawing center of mass scatterplot figure.');

    plt.scatter(avgs, y, c="k")
    ax = plt.gca()
    ax.set_yticks(range(length))
    ax.set_yticklabels(labels)
    plt.errorbar(avgs, y, xerr=errors, linestyle="None", ecolor="k")
    plt.title('Chromosomal Probe Densities')
    plt.ylabel('Chromosome')
    plt.xlabel('Distance')
    plt.grid(True)

    logger.stdout("Saving figure to vis/com/stddev.png")
    plt.savefig('vis/com/stddev.png', dpi=600)
    plt.close()

if __name__ == "__main__":
    plotComScatter()
