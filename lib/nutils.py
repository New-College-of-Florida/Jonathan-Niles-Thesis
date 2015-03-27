#!/usr/bin/env python
"""
nutils.py - utilities making analysis easier
copyright (c) 2015  Jonathan Niles

this program is free software: you can redistribute it and/or modify
it under the terms of the gnu general public license as published by
the free software foundation, either version 3 of the license, or
(at your option) any later version.

this program is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  see the
gnu general public license for more details.

you should have received a copy of the gnu general public license
along with this program.  if not, see <http://www.gnu.org/licenses/>.
"""

import os
from collections import namedtuple

Annotation = namedtuple('Annotation', ('chromosome', 'start', 'end'))
sync = "/home/jniles/thesis/sync/"
datasets = {
    "IMR90" : ["R1", "R2", "R3", "R4", "R5", "R6"],
    "hESC" : ["R1", "R2"]
}

def join(*args):
    """alias for os.path.join"""
    return os.path.join(*args)

def chkdir(directory):
    """ensures that a directory exists, creating it if it does not."""
    d = os.path.dirname(directory)
    if not os.path.exists(d):
        print "%s does not exist.  Creating..." % directory
        os.makedirs(d)
    return directory

def strToResolution(resString):
    assert 'k' in resString
    kb = resString.split('k')[0]
    return int(kb)*1000

def sliceUCSCCoords(index, model, **kwargs):
    """
    Converts UCSC genome index notation into indices for the interaction
    matrix, and slices in at the given resolution.

    Arguments
    index -- A UCSC-style string (e.g. "chr6:117,607,316-117,744,804")
    """
    resolution = model.resolution
    region = parseUCSCAnnotation(index)
    # Get the start bins and end bins
    startBin = region.start / resolution
    endBin = region.end / resolution + 1
    # filter the chromosome
    dataArray = filterByChromosome(model, region.chrom, kwargs)
    return dataArray[startBin:endBin].T[startBin:endBin]


def filterByChromosome(data, chromosome, key="key"):
    """filters a binnedDataSet by the chromosome number provided.
    
    :params data: binnedDataSet
    :params chromosome: chromosome number
    :params key: key into binnedDataset
    """
    assert chromosome in data.chromosomeIndex
    # get a boolean mask for bins that correspond to
    # the desired chromosome
    mask = data.chromosomeIndex == chromosome

    # Pull out the raw (binned) matrix
    rawMatrix = data.dataDict[key]

    return squareMask(rawMatrix, mask)

def squareMask(data, mask):
    """
    Returns the square mask of a nxn matrix
    by slicing out the rows and columns.

    Arguments:
    data -- A numpy matrix
    mask -- An array of booleans to mask the data matrix
    """
    return data[mask].T[mask]

def parseUCSCAnnotation(annotation):
    """
    parses a UCSC-style annotation
    """
    try:
        chrom = int(annotation[3:annotation.find(':')])
        startStr = annotation[annotation.find(':') + 1: annotation.find('-')]
        endStr = annotation[annotation.find('-') + 1: ]
        start = int(startStr.replace(',', ''))
        end = int(endStr.replace(',', ''))
    except:
        raise Exception("This doesn't seem like UCSC-style coordinates.")
    return Annotation(chromosome=chrom, start=start, end=end)
