"""
iterativeMapping.py
Author : mirnylib
Adapted by: jniles
October 7, 2014
"""

"""
This scripts takes .sra files from fastq directory, maps them to the genome and
saves them to .hdf5 files in a directory "genomeName-mapped".
Please follow comments along the text.
"""

import atexit
import time
import glob
import os
import logging
from hiclib import mapping
from mirnylib import h5dict, genome

import numpy as np
logging.basicConfig(level=logging.DEBUG)

def cleanFile(filename):
    if os.path.exists(filename):
        os.remove(filename)

# cell typs is either hESC or IMR90
cellType = "hESC"
genomeName = "hg19"
threads = 14 # result `nproc`
bowtiePath = "../bin/bowtie2/bowtie2"
if not os.path.exists(bowtiePath): raise

sraDir = "/home/jniles/data/{0}/sra/{1}"
replicates = ["r1", "r2"]


bowtieIndex = "../bin/bowtie2/index/{0}".format(genomeName)
tmpDir = "/tmp"

samFolder = "sams-{0}".format(genomeName)
savePath = "mapped-{0}".format(genomeName)

# Specify location of the genome files here
dnaDb = genome.Genome('/home/jniles/data/dna/{0}'.format(genomeName), readChrms=["#", "X"])

if not os.path.exists(samFolder):
    os.mkdir(samFolder)

if not os.path.exists(savePath):
    os.mkdir(savePath)

def calculateStep(length, minlen, approxStep=10, maxSteps=4):
    """returns minimum length and step based on the
    length of sequence and proposed minimum length"""

    actualDif = length - minlen
    if actualDif < approxStep * 0.6:
        return length, 100

    numIter = np.array(np.around(actualDif / float(approxStep)), dtype=int)
    if numIter == 0:
        numIter = 1
    if numIter > maxSteps:
        numIter = maxSteps
    actualStep = actualDif / numIter

    minlen = length - actualStep * numIter

    return minlen, actualStep

def iterativeMapDir(directory):
  # three times because we want to get sure that if some other process got
  # killed and left an un-mapped file (which we skipped because it had a lock
  # on it), then we would map it as well.
  files = os.listdir(directory)
  print "[ITER] Found %i files " % len(files)

  if files.count("EXPERIMENT") > 0:
    print "[ITER] Removing descriptor file ..."
    files.remove("EXPERIMENT")

  for i in 3 * sorted(files):
      expName = i
      print "[ITER] Mapping file {0}".format(expName)
      file1 = os.path.join(directory, expName)
      if not os.path.exists(file1): raise
      if not os.path.exists(file1): raise
      lengthFile = os.path.join("preprocessing/lengths", expName)
      length = (int(open(lengthFile).readlines()[0]) - 1) / 2
      print "[ITER] Found length {0} for {1}".format(length, expName)
      minlen, step = calculateStep(length, 25)


      finalName = '%s/%s.hdf5' % (savePath, expName.replace(".sra", ""))
      lockName = finalName + ".lock"
      print "[ITER] Writing file to {0}.".format(finalName)

      if os.path.exists(finalName) and not os.path.exists(lockName):
          print "skipping", finalName
          continue

      if os.path.exists(lockName):
          print "someone is working on", finalName
          continue

      lock = open(lockName, "w")
      lock.close()

      atexit.register(cleanFile, lockName)

      os.system("rm -rf {0}/{1}*".format(samFolder, expName.replace(".sra", "")))

  # First step. Map the reads iteratively.
      mapping.iterative_mapping(
          bowtie_path=bowtiePath,
          bowtie_index_path=bowtieIndex,
          fastq_path=file1,
          out_sam_path='{0}/{1}_1.bam'.format(samFolder, expName),
          min_seq_len=minlen,  # for bacteria mimimal mappable length is 15 bp, so I start with something slightly longer
          len_step=step,  # and go with a usualy step
          nthreads=threads,  # on intel corei7 CPUs 4 threads are as fast as
                       # 8, but leave some room for you other applications
          # max_reads_per_chunk = 10000000,  #optional, on low-memory machines
          temp_dir=tmpDir,
          seq_start=0,
          seq_end=length,
          bash_reader="fastq-dump -Z",
          bowtie_flags=" --very-sensitive ",
          )

      mapping.iterative_mapping(
          bowtie_path=bowtiePath,
          bowtie_index_path=bowtieIndex,
          fastq_path=file1,
          out_sam_path='{0}/{1}_2.bam'.format(samFolder, expName),
          min_seq_len=minlen,
          len_step=step,
          nthreads=threads,  # on intel corei7 CPUs 4 threads are as fast as
                       # 8, but leave some room for you other applications
          # max_reads_per_chunk = 10000000,  #optional, on low-memory machines
          temp_dir=tmpDir,
          seq_start=length,
          seq_end=2 * length,
          bash_reader="fastq-dump -Z",
          bowtie_flags=" --very-sensitive ",
          )

      # Second step. Parse the mapped sequences into a Python data structure,
      #    assign the ultra-sonic fragments to restriction fragments.
      mapped_reads = h5dict.h5dict(finalName)

      sam_base1 = '{0}/{1}_1.bam'.format(samFolder, expName)
      sam_base2 = '{0}/{1}_2.bam'.format(samFolder, expName)

      mapping.parse_sam(sam_basename1=sam_base1,
              sam_basename2=sam_base2,
              out_dict=mapped_reads,
              genome_db=dnaDb,
              save_seqs=False)

      os.remove(lockName)
  return

# For each replicate, iteratively map the replicates
print "[INIT] Starting iterative mapping at {0}".format(time.asctime())
for r in replicates:
  inDir = sraDir.format(cellType, r)
  print "[INIT] Iteratively mapping replicate {0}".format(r)
  iterativeMapDir(inDir)

print "[INIT] Ending iterative mapping at {0}".format(time.asctime())
