#!/usr/bin/env python
"""
uses pybedtools to calculate the nearest domains between all replicates in a
cell line, between peaks in the same replicate, and a the same, shuffled
replicate.
"""

import numpy as np
from pybedtools import BedTool

betweenReplicates = BedTool('IMR90-10kb.all.closest.bed')
self = BedTool('IMR90-10kb.self.closest.bed')
shuffled = BedTool('IMR90-10kb.shuffle.closest.bed')

distanceReps = [int(i[-1]) for i in betweenReplicates]
distanceSelf = [int(i[-1]) for i in self]
distanceShuffled = [int(i[-1]) for i in shuffled if i[-1] != '.'] # why does this happen?

print np.median(distanceReps)
print np.median(distanceSelf)
print np.median(distanceShuffled)
