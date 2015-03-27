#!/usr/bin/env python
"""
preprocess.py
preprocess the raw expression values from Affy U+133v2 arrays and aligns them
with the gene dataset from biomart.

Copyright (C) 2015 Jonathan Niles

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
from matplotlib.mlab import csv2rec, rec2csv
from numpy.lib.recfunctions import append_fields

def transformData():
    """preprocessing"""
    names = [
        'name',             # Associated Gene Name
        'strand',
        'end',              # Gene End (bp)
        'start',            # Gene Start (bp)
        'ensemble_id',      # Ensembl Gene ID
        'affy_id',          # Affy HG U133-PLUS-2 probeset
        'chrom'             # Chromosome Name
    ]

    fname = nu.join(nu.sync, 'data/genes/affy_gene.txt')
    print("Loading gene data from", fname)
    genes = csv2rec(fname, delimiter='\t', names=names, skiprows=1)

    # add fields for the new data sets
    shape, = genes.shape
    genes = append_fields(genes, 'hesc_1', np.zeros(shape), usemask=False)
    genes = append_fields(genes, 'hesc_2', np.zeros(shape), usemask=False)
    genes = append_fields(genes, 'imr90_1', np.zeros(shape), usemask=False)
    genes = append_fields(genes, 'imr90_2', np.zeros(shape), usemask=False)

    # load data
    print("Loading cell line probe sets")
    stem = nu.join(nu.sync, 'data/genes/')
    hesc_1 = csv2rec(nu.join(stem, 'hesc/GSM1309417.txt'), delimiter='\t')
    hesc_2 = csv2rec(nu.join(stem, 'hesc/GSM1309418.txt'), delimiter='\t')
    imr90_1 = csv2rec(nu.join(stem, 'imr90/GSM51626.txt'), delimiter='\t')
    imr90_2 = csv2rec(nu.join(stem, 'imr90/GSM51627.txt'), delimiter='\t')

    print("Sorting probe sets")
    id = 'id_ref'
    hesc_1.sort(order=id) # speed gains!
    hesc_2.sort(order=id)
    imr90_1.sort(order=id)
    imr90_2.sort(order=id)

    print("Assigning probe expression values")
    v = 'value' # for conciseness
    for i, probe in enumerate(genes['affy_id']):
        genes[i]['hesc_1'] = hesc_1[hesc_1[id] == probe][v]
        genes[i]['hesc_2'] = hesc_2[hesc_2[id] == probe][v]
        genes[i]['imr90_1'] = imr90_1[imr90_1[id] == probe][v]
        genes[i]['imr90_2'] = imr90_2[imr90_2[id] == probe][v]

    print("Adjusting chromosome labels")
    # makes sure chromosmes are not just number
    dtypes = np.dtype([
        ('name', 'S28'),
        ('strand', '<i8'),
        ('end', '<i8'),
        ('start', '<i8'),
        ('ensemble_id', 'S15'),
        ('affy_id', 'S27'),
        ('chrom', 'S5'),
        ('hesc_1', '<f8'),
        ('hesc_2', '<f8'),
        ('imr90_1', '<f8'),
        ('imr90_2', '<f8')
    ])

    genes = genes.astype(dtypes)
    genes['chrom'] = np.array(map(lambda x: 'chr'+x, genes['chrom']))

    # reorder columns so this fits as a .bed file
    cols = ['chrom', 'start', 'end', 'ensemble_id', 'strand', 'name', 'hesc_1', 'hesc_2', 'imr90_1', 'imr90_2']
    bed = genes[cols]

    fname = nu.join(nu.sync, 'data/genes/expression.bed.gz')
    print("Writing new file to", fname)
    rec2csv(bed, fname, delimiter='\t')
    return

if __name__ == "__main__":
    transformData()
