#!/usr/bin/env python

#pScatter
#Interaction visualizer for pLink data, by Justin Jee (with design by Katelyn McGary Shipper)
#June 2015
#
#input: <fasta protein file> <alphabet of linkable amino acids> <dir1> <dir2> ... <output name>
#where <dirN> is a directory containing all the pLink txt outputs for a single experiment 
#output: summary txt file + heatmaps corresponding to the pairwise prot interactions
#example: python parseraw.py test.fasta test1 test2 testout
#generates: test2.txt, test2.pdf
#
#parses data from column that looks like name1(pos)-name2(pos)
#IMPORTANT NOTE: If saving input txt file from excel it is important to save in Windows txt format
#dependencies: numpy, matplotlib

import loadfiles
import splotch
import sys
import argparse
import numpy as np
from collections import defaultdict, OrderedDict
from bokeh.plotting import figure
from bokeh.io import gridplot, output_file, show

parser = argparse.ArgumentParser()
parser.add_argument('fasta_file', type=str, help='fasta file with protein sequences')
parser.add_argument('data_dirs', nargs='+', help='pLink output directory')
parser.add_argument('-a', '--aminoacids', default='K', help='cross-linkable aminoacids. Defaults to Lysine (K).')
parser.add_argument('-z', '--zoom', help='Prot1-Prot2 only display subplot for proteins Prot1 vs Prot2',
                    type=str)
parser.add_argument('-s', '--scale', help='scale both plot and outputs so that only amino acids of interest are considered', action='store_true')
parser.add_argument('-e', '--evalue', default=1.0, help='e-value cutoff', type=float)
parser.add_argument('-o', '--output', default='output.html', help='output file (HTML) name')
parser.add_argument('-v', '--verbose', help='increase output verbosity', action='store_true')
args = parser.parse_args()

### load fasta file ###
prot2map = loadfiles.loadfasta(args.fasta_file, args.aminoacids)
### load plink file ###
print "fasta loaded. loading pLink files"
allprot = defaultdict(int)
interactions = []
for dirname in args.data_dirs:
    gg2i, allprot = loadfiles.loadplink(dirname, prot2map, allprot, args.scale, args.evalue)
    interactions.append(gg2i)

### output ###
print "generating plot"

#Define basic plot parameters
color = ['blue','red','green']
marker = ['o','o','o']
basesize = 50 
buffer = 10
numprot = len(allprot)
#sort proteins by decreasing size
sallprot = sorted(allprot.items(), key=lambda t:t[1])
sallprot = [k[0] for k in sallprot]
sallprot = sallprot[::-1]

for interaction in interactions:
    print interaction
    print '\n'

#Now do actual plotting
if args.output.endswith('.html'):
    output_file(args.output) 
else:
    output_file(args.output + '.html')

if not args.zoom:
    m = [[None] * (numprot + 1)] * numprot
    for xind,interaction in enumerate(interactions):
        m[0][-1] = splotch.makelegend(m[0][-1], numprot, color[xind], marker[xind], args.data_dirs[xind])
        i =  numprot - 1
        for prot2 in sallprot:
            j = 0
            for prot1 in sallprot:
                key = prot1 + '-' + prot2
                (x, y, r, mc) = loadfiles.writesummary(key, interaction, prot2map, args.scale)
                m[i][j] = splotch.splotch(m[i][j], m[numprot-1][0], x, y, r, basesize, buffer, mc, key, 
                        i==numprot-1, j==0, numprot, color[xind], marker[xind], None, False, args.unjoin)
                j+=1
            i-=1
    p = gridplot(m)
else:
    for xind,interaction in enumerate(interactions):
        (x, y, r, mc) = loadfiles.writesummary(args.zoom, interaction, prot2map, args.scale)
        p = splotch.splotch(None, None, x, y, r, basesize, buffer, mc, args.zoom, True, True, 1,
                color[xind], marker[xind], args.data_dirs[xind], True)

show(p)
