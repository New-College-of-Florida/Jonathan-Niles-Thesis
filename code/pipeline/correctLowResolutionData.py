#!/usr/bin/env python
"""
This is a draft of a script which corrects one or multiple Hi-C datasets. 
It uses BinnedData, meaning that you need to have resolution of 1000kb or more
"""

from __future__ import print_function

from hiclib.binnedData import binnedData

# names to be mapped
fnames = [
    "/home/jniles/data/IMR90/raw/IMR90-all-HindIII-1000k.hm",
    "/home/jniles/data/IMR90/raw/IMR90-R1-HindIII-1000k.hm",
    "/home/jniles/data/IMR90/raw/IMR90-R2-HindIII-1000k.hm",
    "/home/jniles/data/IMR90/raw/IMR90-R3-HindIII-1000k.hm",
    "/home/jniles/data/IMR90/raw/IMR90-R4-HindIII-1000k.hm",
    "/home/jniles/data/IMR90/raw/IMR90-R5-HindIII-1000k.hm",
    "/home/jniles/data/IMR90/raw/IMR90-R6-HindIII-1000k.hm",
    "/home/jniles/data/hESC/raw/hESC-R1-HindIII-1000k.hm",
    "/home/jniles/data/hESC/raw/hESC-R2-HindIII-1000k.hm",
    "/home/jniles/data/hESC/raw/hESC-all-HindIII-1000k.hm"
]

path = "/home/jniles/data/{0}/raw/{0}-{1}-HindIII-{2}.hm"
resolutions = {
    "2000k" : 2000000,
    "1000k" : 1000000,
    "500K"  : 500000,
    "200K"  : 200000,
}

names = map(lambda s: s.split("/")[-1], fnames)
exportnames = map(lambda s: s.replace("raw", "ic"), fnames)

resolution = 1000000
genome = "/home/jniles/data/dna/hg19/"

def main():
    data = binnedData(resolution, genome)

    for name,fname,exportname in zip(names,fnames,exportnames):
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

    for name, exportname in zip(names, exportnames): 
        print("Saving", exportnames)
        data.export(name, exportname)
    return

if __name__ == "__main__":
    main()
