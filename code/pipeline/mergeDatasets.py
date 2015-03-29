#!/usr/bin/env python
"""
This script does the following:

- Loads individual files from the location defined by source(ID)
- Parses individual files in memory (inMemory = True)
   - If your system does not have enough memory, you might need to switch to
     hdf5 here.
- Merges files corresponding to the same experiment together, on the HDD.
- Filters datasets, builds heatmaps
- Combines multiple replicas of the same experiment together, builds heatmaps

--Locations of individual files are defined by source(ID) function.
--Datasets are defined in the datasets.tsv file
--genome is defined by genomeFolder function, and workingGenome identifyer
--output files are arranged to folders named by their workingGenome IDs

Warnings:
    Running this over NFS might cause unexpected slow-downs because NFS is
    unhappy with repeated read/write cycles to the same file

    You could do entire thing in memory, if you have RAM or your datasets are
    small.  Actually, using HDF5 is then equivalent to storing compressed data
    in RAM, and might be in fact pretty fast.

modified by: jniles
"""

import os
import sys
import nutils as nu
from hiclib.fragmentHiC import HiCdataset

def genomeFolder(name):
    return nu.join("/home/jniles/data/dna", name)

byChromosomeResolutionsKb = [100, 40]
HiResWithOverlapResolutionsKb = [20, 10]
wholeGenomeResolutionsKb = [2000, 1000,500,200]

def refineDataset(filenames, create=True, delete=True, parseInMemory=True):
    """
    Parameters
    ----------

    filenames[0] is a list of filenames of incoming files
    filenames[1] is a folder for outgoing file
    filenames[2] is a working genome, that is output directory
    filenames[3] is an enzyme for a given experiment


    create : bool, optional
        If True, parse each file.
        If False, assume that files were already parsed
        (e.g. if you are just playing around with filtering parameters)
    delete : bool, optional
        If True, delete parsed files after merging.
        Man, these files may be huge... if you don't have a 10TB RAID, this may be useful.
    parseInMemory : bool, optional
        Perform parsing input files in memory.
    """
    in_files = filenames[0]
    out_file = filenames[1]

    statFolder = nu.join("statistics", out_file)

    workingGenome = filenames[2]
    enzyme = filenames[3]


    if create == True:  # if we need to parse the input files (.hdf5 from mapping).
        for onename in in_files:

            # Parsing individual files
            if parseInMemory == True:

                # create dataset in memory, parse and then save to destination
                TR = HiCdataset("bla", genome=genomeFolder(workingGenome),
                                maximumMoleculeLength=500,
                                inMemory=True)  # remove inMemory if you don't have enough RAM


                TR.parseInputData(dictLike=onename, enzymeToFillRsites=enzyme)

                folder = os.path.split(onename)[0]
                print onename
                TR.save(nu.chkdir(onename + "_parsed.frag"))
                folder, fname = os.path.split(onename)
                statSubFolder = os.path.join(statFolder, folder)

                TR.printMetadata(saveTo=nu.chkdir(os.path.join(statSubFolder, fname + ".stat")))
            else:
                #Create dataset at destination, parse on HDD, then no need to save.
                TR = HiCdataset(nu.chkdir(onename + "_parsed.frag"),
                                genome=genomeFolder(workingGenome),
                                maximumMoleculeLength=500, mode='w')
                TR.parseInputData(dictLike=onename, enzymeToFillRsites=enzyme)
                TR.printMetadata(saveTo=nu.chkdir(os.path.join(statFolder, onename + ".stat")))

        "Merging files alltogether, applying filters"
        TR = HiCdataset(nu.chkdir(out_file + "_merged.frag"),
                        genome=genomeFolder(workingGenome),
                        mode="w")
        TR.merge([i + "_parsed.frag" for i in in_files])
            #Merge in all parsed files from one experiment

        if delete == True:  # cleaning up parsed files
            for delFile in [i + "_parsed.frag" for i in in_files]:
                os.remove(delFile)

        "Now opening new dataset for refined data, and performing all the filtering "
        TR = HiCdataset(out_file + "_refined.frag",
                        genome=genomeFolder(workingGenome),
                        mode='w')
        TR.load(out_file + "_merged.frag")

        #----------------------------Set of filters applied -------------
        TR.filterRsiteStart(offset=5)
        TR.filterDuplicates()
        #TR.save(out_file+".dat")
        TR.filterLarge()
        TR.filterExtreme(cutH=0.005, cutL=0)
        TR.writeFilteringStats()
        TR.printMetadata(saveTo=statFolder + ".stat")

        #------------------------End set of filters applied----------

    else:
        #If merging & filters has already been done, just load files
        TR = HiCdataset(out_file + "_working.frag",
                        mode='w', genome=genomeFolder(workingGenome))
        TR.load(out_file + "_refined.frag")
        TR.printMetadata(saveTo=statFolder + ".stat")

    print "-----> Building Raw heatmap at different resolutions"
    TR.printStats()
    for res in wholeGenomeResolutionsKb:   
        TR.saveHeatmap(out_file + "-{0}k.hm".format(res), res*1000)
    for res in byChromosomeResolutionsKb:
        TR.saveByChromosomeHeatmap(out_file + "-{0}k.byChr".format(res), res*1000)
    for res in HiResWithOverlapResolutionsKb:
        TR.saveHiResHeatmapWithOverlaps(out_file + "-{0}k_HighRes.byChr".format(res), res*1000)

def main(filename):
    """runs when called from the command line"""

    # This code is actually parsing datasets.tsv file
    print "[MERGE] Opening file: {0}".format(filename)

    fsplit = os.path.split(filename)
    if len(fsplit[0]) > 0:
        os.chdir(fsplit[0])
    filename = fsplit[1]

    # open and parse the datasets file
    lines = open(filename).readlines()
    lines = [i for i in lines if i[0] != "#"]
    lines = [i.split() for i in lines if (len(i) > 3) and (i[0] != "#")]
    print lines

    for i in lines:
        if len(i) == 4:
            # Now we assume that enzyme is fifth argument in datasets.tsv or second commandline argument
            try:
                i.append(sys.argv[2])
            except:
                err = """
                      Please specify enzyme as a second command line argument,
                      or as a fifth column
                      """
                print err
                raise ValueError(err)

    # Make sure we have an enzyme for each file
    assert False not in [len(i) == 5 for i in lines]

    dataFiles = lines

    #experiment is defined by experiment name, replica name, genome and enzyme
    experimentNames = set((i[1], i[2], i[3], i[4]) for i in dataFiles)
    byExperiment = []
    combinedExperimentNames = []

    for experiment in experimentNames:

        #experiment is IMR90,R1,hg19,HindIII
        workingGenome = experiment[2]
        enzyme = experiment[3]
        filenames = [i[0] for i in dataFiles if (i[1], i[2], i[3], i[4]) == experiment]
        outName = "{0}-{1}-{3}".format(*experiment)

        #Filenames, path to save, genome, enzyme
        byExperiment.append((filenames, os.path.join(workingGenome, outName), workingGenome, enzyme))
        # merged files belonging to one expriment
        combinedExperimentNames.append((experiment[0], os.path.join(workingGenome, outName), workingGenome, enzyme))


    # Now running refineDataset for each experiment
    for i in byExperiment:
        print i
        refineDataset(i, create=True, delete=True)

    # Now merging different experiments alltogether
    # note that the first column is not here, as it is a replica
    experiments = set([(i[0], i[2], i[3]) for i in combinedExperimentNames])
    print experiments

    for experiment in experiments:
        workingGenome = experiment[1]
        myExperimentNames = [i[1] + "_refined.frag" for i in combinedExperimentNames if (i[0], i[2], i[3]) == (experiment[0], experiment[1],experiment[2])]   
        assert len(myExperimentNames) > 0
        if len(myExperimentNames) > 1:
            #If we have more than one experiment (replica) for the same data, we can combine.
            TR = HiCdataset(os.path.join(workingGenome, "%s-all-%s_refined.frag" %
                                         (experiment[0],experiment[2])), genome=genomeFolder(workingGenome))
            statSaveName = os.path.join("statistics", workingGenome, "%s-all-%s_refined.stat" % (experiment[0], experiment[2]))

            TR.merge(myExperimentNames)
            TR.printMetadata(saveTo=statSaveName)
            for res in wholeGenomeResolutionsKb:   
                TR.saveHeatmap(os.path.join(workingGenome, "%s-all-%s-{0}k.hm" % (experiment[0], experiment[2])).format(res), res*1000)
            for res in byChromosomeResolutionsKb:
                TR.saveByChromosomeHeatmap(os.path.join(workingGenome, "%s-all-%s-{0}k.byChr" % (experiment[0], experiment[2])).format(res), res*1000)
            for res in HiResWithOverlapResolutionsKb:
                TR.saveHiResHeatmapWithOverlaps(os.path.join(workingGenome, "%s-all-%s-{0}k_HighRes.byChr" % (experiment[0], experiment[2])).format(res), res*1000)


if __name__ == "__main__":
    try:
        filename = sys.argv[1]
    except:
        filename = "datasets.tsv"
    main(filename)
