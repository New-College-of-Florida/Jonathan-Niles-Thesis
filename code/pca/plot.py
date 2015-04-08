#!/usr/bin/env python
"""
plot.py
Plots all related component plots to the sync folder

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

import numpy as np
import nutils as nu
import matplotlib.pyplot as plt

screeSavePath = nu.chkdir(nu.join(nu.sync, "plots/pca/scree/"))
compSavePath =  nu.chkdir(nu.join(nu.sync, "plots/pca/pc/"))
tmpl = nu.join(nu.sync, "data/components/{0}-{1}-200k.{2}")

def loadComponent(cellType, replicate, dtype="eigenvectors"):
    """
    Loads the data for a certain cell type and replicate
    """
    path = tmpl.format(cellType,replicate, dtype)
    print("Loading file:", path)
    return np.loadtxt(path)

def plotPrincipalComponent(cellType, replicate, components, number=0):
    """Plots a the principal component for a celltype, replicate"""

    fig, ax = plt.subplots(figsize=(22, 2))

    N = len(components[number])

    # the eigenvalues can be negative, so plot the magnitude of the eigenvalue
    x = np.arange(N)
    ax.plot(x, components[number], color="k", linewidth=0.75)
    #ax.fill_between(x, 0, components[number], color="grey")

    # pretty axis limits
    ax.set_ylabel("$\\vec{e}$")
    #plt.suptitle("{0} {1}".format(cellType, replicate))
    ax.set_title("$E_{%i}$" % (number+1))

    m = np.max(components[number])
    ax.set_yticks([-(m + 0.1*m), (m + 0.1*m)])
    ax.set_yticklabels(["-", "+"])
    
    p = N / 20 # put ticks at 5% and 95%, respectively
    ax.set_xticks([p, N - p])
    ax.set_xticklabels(["Chromosome 1", "Chromosome X"])

    plt.tight_layout()

    fname = "{0}-{1}.png".format(cellType, replicate)
    outFig = nu.join(compSavePath, fname)
    print("Saving to:", outFig)
    fig.savefig(outFig, dpi=700)
    plt.close()
    return

def plotScree(cellType, replicate, data):
    """Plots a form of scree plot for eigenvalues"""

    fig, ax = plt.subplots()
    length = len(data)

    # the eigenvalues can be negative, so plot the magnitude of the eigenvalue
    ax.plot(np.arange(1, length+1), np.abs(data), 'ok')

    # pretty axis limits
    ax.set_xlim(0, length + 1)
    ax.set_xlabel("PC Number")
    ax.set_ylabel("$|\lambda|$")
    ax.set_title("{0} {1} Scree Plot".format(cellType, replicate))
    ax.set_xticks(np.arange(length + 1))
    ax.set_xticklabels(np.arange(length + 1))

    ax.minorticks_on()

    saveFile = "{0}-{1}.png".format(cellType, replicate)
    print("Saving to:", saveFile)
    fig.savefig(nu.join(screeSavePath, saveFile))
    plt.close()
    return

def plotAllScree():
    """plots all scree plots"""
    for cellType in nu.datasets.keys():
        for rep in nu.datasets[cellType]:
            eigenvalues = loadComponent(cellType, rep, dtype="eigenvalues")
            print("Plotting scree plots for", cellType, rep)
            plotScree(cellType, rep, eigenvalues)
    return

def plotAllComponents():
    """plots all the component plots"""
    for cellType in nu.datasets.keys():
        for rep in nu.datasets[cellType]:
            eigenvectors = loadComponent(cellType, rep)
            print("Plotting PC plots for", cellType, rep)
            plotPrincipalComponent(cellType, rep, eigenvectors)
    return

if __name__ == "__main__":
    #plotAllScree()
    plotAllComponents()
