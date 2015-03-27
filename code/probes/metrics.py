#!usr/bin/env python

"""
This file has a number of metrics to describe the distance of probe scaling.

last modified: Feb 26 2015
author: jniles
"""

import os
import pickle
import nutils
import numpy as np
import matplotlib.pyplot as plt

path = "data/"

def loadScalingData(cellType):
    """
    Loads scaling by replicate from preprocessed pickle data
    """
    datasets = replicates[cellType]
    dataDict = {}
    for rep in datasets:
        inFile = os.path.join(data, '{0}-{1}.pkl'.format(cellType, rep))
        with open(inFile, 'r') as f:
            dataDict[rep] = pickle.load(f)
    return dataDict

def aggregate(array)
    length = len(array)
    agg = np.zeros(length)
    for i in xrange(length):
        if i == 0:
            agg+= array[i]
        else:
            agg[i] += agg[i-1] + array[i]
    return agg

def halfMetrics(cellType, ctData):
    """
    Plots the half-max metrics for each chromosome
    """
    for rep in dataDict.keys():
        r = dataDict[rep][rep]
        contacts = aggregate(r[chrom]['contacts'])

    return

def main(cellType):
    """
    Runs the module
    """
    # make the save directory if it doesn't exist
    nutils.chkdir(savePath)

    # ctData = loadData(cellType)
    # plotChromosomeScaling(cellType, ctData)
    ctData = loadScalingData(cellType)
    plotMetrics(cellType, ctData, percentage=True)
    return

if __name__ == "__main__":
    main("IMR90")
    main("hESC")
