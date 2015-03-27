#!/usr/bin/env python
"""
Plots domains on human chromosomes
"""

import numpy as np
from matplotlib.mlab import csv2rec
from reportlab.lib.units import cm
from Bio.Graphics import BasicChromosome

# chromosome sizes
entries = [
    ('Chr1', 249250621),
    ('Chr2', 243199373),
    ('Chr3', 198022430),
    ('Chr4', 191154276),
    ('Chr5', 180915260),
    ('Chr6', 171115067),
    ('Chr7', 159138663),
    ('ChrX', 155270560),
    ('Chr8', 146364022),
    ('Chr9', 141213431),
    ('Chr10', 135534747),
    ('Chr11', 135006516),
    ('Chr12', 133851895),
    ('Chr13', 115169878),
    ('Chr14', 107349540),
    ('Chr15', 102531392),
    ('Chr16', 90354753),
    ('Chr17', 81195210),
    ('Chr18', 78077248),
    ('Chr20', 63025520),
    ('Chr19', 59128983),
    ('Chr22', 51304566),
    ('Chr21', 48129895),
]

# templates
bedTmpl = "beds/{0}-{1}-{2}kb.bed"
columns = ['chrom', 'start', 'end', 'value', 'sign', 'length', 'id']

def parseBed(cellType, replicate, windowSize):
    """
    parses a bed file given in the template
    """
    fname = bedTmpl.format(cellType, replicate, windowSize)
    return csv2rec(fname, delimiter='\t', names=columns)

def plotOnChrom(domains, chromosomes):
    chr_diagram = BasicChromosome.Organism()
    chr_diagram.page_size = (29.7*cm, 21*cm) #A4 landscape

    max_len = np.max([c[1] for c in chromosomes])
    telomere_length = 1000000 #For illustration

    for idx, (name, length) in enumerate(chromosomes):
        cur_chromosome = BasicChromosome.Chromosome(name)
        features = domains[domains['chrom'] == name.lower()]

        #Set the scale to the MAXIMUM length plus the two telomeres in bp,
        #want the same scale used on all five chromosomes so they can be
        #compared to each other
        cur_chromosome.scale_num = max_len + 2 * telomere_length

        #Add an opening telomere
        start = BasicChromosome.TelomereSegment()
        start.scale = telomere_length
        cur_chromosome.add(start)

        #Add a body - using bp as the scale length here.
        body = BasicChromosome.AnnotatedChromosomeSegment(length, features)
        body.scale = length
        cur_chromosome.add(body)

        #Add a closing telomere
        end = BasicChromosome.TelomereSegment(inverted=True)
        end.scale = telomere_length
        cur_chromosome.add(end)

        # This chromosome is done
        chr_diagram.add(cur_chromosome)


    return

if __name__ == "__main__":
    cellType = "IMR90"
    replicate = "R1"
    windowSize = 40
    domains = parseBed(cellType, replicate, windowSize)
    plotOnChrom(domains, entries[:5]) # plot the first five chromosomes
