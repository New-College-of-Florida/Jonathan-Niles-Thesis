#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt 
import csv

rand = False 

if rand:
    fname = 'imr90.r1.vs.random.bed'
else:
    fname = 'imr90.r1.vs.r2.bed'

reader = csv.reader(open(fname), delimiter='\t')
data = np.array([row for row in reader])

distances = [np.int(row[-1]) for row in data]

fig, ax = plt.subplots()

ax.hist(distances, bins=500, normed=True)

median = np.median(distances)
label = 'Median: ~{0:.2f}kb'.format(median / 1000.0) # format to kb
fig.text(0.75, 0.75, label, horizontalalignment='center', verticalalignment='center')


# set tick labels to be 1MB intervals
size = 1000000 # 1MB
ax.set_xlim(0, size*10)
ax.set_xticklabels([str(n / 100000) for n in ax.get_xticks()])

if rand:
    ax.set_title('IMR90 R1 vs Random')
else:
    ax.set_title('IMR90 R1 vs R2 ')

ax.set_xlabel('Distance (MB)')
ax.set_ylabel('Frequency')

plt.tight_layout()

if rand:
    fig.savefig('hist.imr90.r1.vs.random.png', dpi=500)
else:
    fig.savefig('hist.imr90.r1.vs.r2.png', dpi=500)
