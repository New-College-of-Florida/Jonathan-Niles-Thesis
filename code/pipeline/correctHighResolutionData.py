#!/usr/bin/env python
"""
Iteratively correct high-res heatmaps.

Author: LMirny
Modified by: jniles
Date: Jan 28th, 2015

"""

# use python 3 printing
from __future__ import print_function

import os
import sys
import sh
import perf
import nutils

from copy import copy
from mirnylib.h5dict import h5dict
from mirnylib.numutils import completeIC
import matplotlib.pyplot as plt
import numpy as np

rawFolderTemplate = "/home/jniles/data/{0}/raw/"
icFolderTemplate = "/home/jniles/data/{0}/ic/"

def correct(filename, cellType):
    print("Reading {0}".format(filename))
    dataset = h5dict(filename,'r')  # open in the "read" mode

    # a bit of a weird way to find chromosome number
    keys = dataset.keys()
    cisKeys = [i for i in keys if len(set(i.split()))==1] # extract keys of the type "a a"
    numChromosomes = len(cisKeys)

    for chromosome in range(numChromosomes):
        chromosomeHeatmap = dataset["{0} {0}".format(chromosome)] # extracting cis heatmap

        # This line executes proper Iterative Correction, which accounts for regions with low coverage
        # It only works for cis (symmetric) heatmaps.
        try:
            correctedHeatmap, bias = completeIC(chromosomeHeatmap,returnBias=True)
        except:
            e = sys.exc_info()[0]
            print("Warning: ", e)
            pass

    # save the dataset
    outDir = nutils.chkdir(icFolderTemplate.format(cellType))
    outFile = os.path.join(outDir, os.path.basename(filename))
    outData = h5dict(outFile,'w') # open in "write" mode
    for key in dataset.keys():
        outData[key] = copy(dataset[key])
    outData.update()
    return

def main(celltype):
    # this only works with Hi-C dataset saved by chromosome

    # set up test timer
    timer = perf.Test("High Res Iterative Correction")
    timer.start()

    folder = rawFolderTemplate.format(cellType)

    # get all the files in the data directory.
    # execute $ find - masdepth 1 -type f -name *.byChr to get high-res files
    try:
        files = sh.find(folder, '-maxdepth', 1, '-type', 'f', '-name', '*.byChr')
    except:
        raise "Error: /usr/bin/find raised an error"
    files = map(lambda f : f.strip(), list(files))

    for f in files:
        correct(f, cellType)

    # print resultant times
    timer.stop()
    timer.results()
    return

if __name__ == "__main__":
    cellType = sys.argv[1]
    print("Correcting all high resolution datasets in {0}".format(cellType))

    # run the module code
    main(cellType)
