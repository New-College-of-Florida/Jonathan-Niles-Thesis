#!/usr/bin/env python
"""
This is a draft of a script which corrects one or multiple Hi-C datasets. 
It uses BinnedData, meaning that you need to have resolution of 1000kb or more
"""

from __future__ import print_function

import nutils as nu
from hiclib.binnedData import binnedData

# paths
genome = "/home/jniles/data/dna/hg19/"
path = "/home/jniles/data/{0}/raw/{0}-{1}-HindIII-{2}.hm"
resolutions = {
    "2000k" : 2000000,
    "1000k" : 1000000,
    "500k"  : 500000,
    "200k"  : 200000,
}

def correctResolution(paths, resolution):
    """corrects a series of maps at a given resolution"""

    print("Correcting all datasets at ", resolution)
    data = binnedData(resolution, genome)
    names = map(lambda s: s.split("/")[-1], paths)
    exportnames = map(lambda s: s.replace("raw", "ic"), paths)

    for name,fname,exportname in zip(names,paths,exportnames):
        print("Loading ...", name)
        data.simpleLoad(fname, name)

    # we never ever use diagonal
    data.removeDiagonal()

    # new filter: omit all bins with less than 0.5 coverage by sequenced bases
    # (i.e. bases present in the genome)
    data.removeBySequencedCount() 

    data.removePoorRegions(cutoff = 0.5, coverage=True)
    # remove .5% bins with the lowest number of records (i.e. non-zero entrees in the matrix)
    # This filter was updated to remove bins which have zero contacts and one PCR blowout.
    # Those bins would have many reads, but all reads will be with one or few other bins. 
    data.removePoorRegions(cutoff = 0.5, coverage=False)  # standard filter. Cutoff reduced to 0.5 from 2. 

    # remove PCR blowouts from trans data
    data.truncTrans()

    print("Iteratively correcting maps")

    # do iterative correction 
    data.iterativeCorrectWithoutSS() 

    print("Saving files to", exportnames)
    for name, exportname in zip(names, exportnames): 
        print("Saving", exportnames)
        data.export(name, exportname)
    return


def main():
    """corrects all resolutions sequentially"""
    for kres, nres in resolutions.items():
        paths = []
        for cellType in nu.datasets.keys():
            for rep in nu.datasets[cellType]:
                paths.append(path.format(cellType, rep, kres))
        print("Found paths:", paths)
        correctResolution(paths, nres)
    return

if __name__ == "__main__":
    main()
