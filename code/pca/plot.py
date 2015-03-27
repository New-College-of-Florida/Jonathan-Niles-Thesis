#!/usr/bin/env python
"""
Plots all related component plots to the sync folder
"""

import numpy as np
import nutils as nu
import matplotlib.pyplot as plt

screeSavePath = nu.chkdir(nu.join(nu.sync, "plots/pca/scree/"))
compSavePath =  nu.chkdir(nu.join(nu.sync, "plots/pca/pc/"))
tmpl = nu.join(nu.sync, "data/components/{0}-{1}-200k.{2}")
replicates = {
    "IMR90": ["R1", "R2", "R3", "R4", "R5", "R6"],
    "hESC": ["R1", "R2"]
}

def loadComponent(cellType, replicate, dtype="eigenvectors"):
    """
    Loads the data for a certain cell type and replicate
    """
    path = tmpl.format(cellType,replicate, dtype)
    return np.loadtxt(path)

def plotPrincipalComponent(cellType, replicate, components, number=0):
    """Plots a the principal component for a celltype, replicate"""

    print "Plotting PC {0}".format(number + 1)

    fig, ax = plt.subplots()

    # the eigenvalues can be negative, so plot the magnitude of the eigenvalue
    x = np.arange(len(components[number]))
    ax.plot(x, components[number], color="k")
    ax.fill_between(x, 0, components[number], color="grey")

    # pretty axis limits
    ax.set_xlabel("Bin")
    ax.set_ylabel("$\\vec{v}$")
    ax.set_title("{0} {1} PC{2}".format(cellType, replicate, number+1))
    ax.minorticks_on()
 
    fname = "{0}-{1}.png".format(cellType, replicate)
    fig.savefig(nu.join(compSavePath, fname), dpi=800)
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
    fig.savefig(nu.join(screeSavePath, saveFile), dpi=800)
    plt.close()
    return

def plotAllScree():
    """plots all scree plots"""
    for cellType in replicates.keys():
        for rep in replicates[cellType]:
            print "Loading components for %s %s" % (cellType, rep)
            eigenvalues = loadComponent(cellType, rep, dtype="eigenvalues")
            print "Plotting scree plots for %s %s" % (cellType, rep)
            plotScree(cellType, rep, eigenvalues)
    return


def plotAllComponents():
    """plots all the component plots"""
    for cellType in replicates.keys():
        for rep in replicates[cellType]:
            print "Loading components for %s %s" % (cellType, rep)
            eigenvectors = loadComponent(cellType, rep)
            print "Plotting PC plots for %s %s" % (cellType, rep)
            plotPrincipalComponent(cellType, rep, eigenvectors)
    return

if __name__ == "__main__":
    plotAllScree()
    plotAllComponents()
