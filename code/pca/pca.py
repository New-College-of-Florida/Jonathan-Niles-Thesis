#!/usr/bin/env python
"""
pca.py
Runs the ICE algorithm as described in Imakeav et al and exports the
component eigenvectors/eigenvalues to numpy arrays.

Copyright (C) 2015  Jonathan Niles

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
from hiclib.binnedData import binnedData

# templates
tmpl = "/home/jniles/data/{0}/ic/{0}-{1}-HindIII-{2}.hm"
genome = "/home/jniles/data/dna/hg19/"
outDir = nu.chkdir(nu.join(nu.sync, "data/components/"))

def savePCAResults(BD, key):
    """saves the results from PCA to some vectors for further analysis"""
    np.savetxt(nu.join(outDir, '{0}.eigenvectors'.format(key), BD.PCDict[key]))
    np.savetxt(nu.join('{0}.eigenvalues'.format(key), BD.PCAEigenvalueDict[key]))
    return

def pca(cellType, resolution):
    """plots the eigenvector plot for a given chromosome"""

    print("Loading", cellType, resolution)

    keys = []
    res = nu.strToResolution(resolution)
    BD = binnedData(res, genome)

    # load each replicate dataset.  They are corrected together
    for rep in nu.datasets[cellType]:
        key = "{0}-{1}-{2}".format(cellType, rep, resolution)
        BD.simpleLoad(tmpl.format(cellType, rep, resolution), key)
        keys.append(key)

    # We never use the diagonal
    BD.removeDiagonal()

    # new filter: omit all bins with less than 0.5 coverage by sequenced bases
    # (i.e. bases present in the genome)
    BD.removeBySequencedCount()

    # remove .5% bins with the lowest number of records
    # (i.e. non-zero entrees in the matrix)
    BD.removePoorRegions(cutoff=0.5, coverage=True)

    # Remove bins with zero counts for PCA analysis
    BD.removeZeros()

    # remove PCR blowouts from trans data
    BD.truncTrans()

    # Fake cis counts. Data gets iteratively corrected during this process...
    BD.fakeCis()

    # Get the first three PCs
    BD.doPCA()

    # Now we restore regions with zero counts.
    # We specify that we fill them in with zeros. By default it goes with NANs.
    for key in keys:
        savePCAResults(BD, key)
    return

if __name__ == "__main__":
    pca("hESC", "200k")
    pca("IMR90", "200k")
