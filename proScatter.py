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
parser.add_argument('-v', '--verbose', help='increase output verbosity', action='store_true')
parser.add_argument('fasta_file', type=str, help='fasta file with protein sequences')
parser.add_argument('data_dirs', nargs='+', help='pLink output directory')
parser.add_argument('-z', '--zoom', help='Prot1-Prot2 only display subplot for proteins Prot1 vs Prot2',
                    type=str)
parser.add_argument('-s', '--scale', help='scale both plot and outputs so that only amino acids of interest are considered', action='store_true')
parser.add_argument('-e', '--evalue', help='e-value cutoff', type-float)
args = parser.parse_args()

#handle user options
args = []
args2 = sys.argv[:]
zoom=False
scale=True #sets default to scaling based on aa only
evalue=1
for i in range(len(args2)):
    a = args2[i]
    if '--zoom' in a:
        zoom=True
        (z,zkey)=a.split('=')
    elif '--scale' in a:
        scale=False
    elif '--evalue' in a:
         (e,evalue)=a.split('=')
         evalue=float(evalue)
    else:
        args.append(a)

output = args[-1]

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

handles = [0] * len(interactions)

#Now do actual plotting
output_file(output+'.html')
if not args.zoom:
    m = [[None for r in range(numprot+1)] for s in range(numprot)]
    for xind in range(len(interactions)):
        m[0][-1]=splotch.makelegend(m[0][-1],numprot,color[xind],marker[xind],dirs[xind])
        i=numprot-1
        for prot2 in sallprot:
            j=0
            for prot1 in sallprot:
                key = prot1 +'-'+ prot2
                gg2i=interactions[xind]
                (x,y,r,mc) = loadfiles.writesummary(key,gg2i,prot2map,scale)
                m[i][j]=splotch.splotch(m[i][j],m[numprot-1][0],x,y,r,basesize,buffer,mc,key,i==numprot-1,j==0,numprot,color[xind],marker[xind],None,False)
                j+=1
            i-=1
    p = gridplot(m)
else:
    for xind in range(len(interactions)):
        gg2i=interactions[xind]
        (x,y,r,mc) = loadfiles.writesummary(zkey,gg2i,prot2map,scale)
        p = splotch.splotch(None,None,x,y,r,basesize,buffer,mc,zkey,True,True,1,color[xind],marker[xind],dirs[xind],True)

show(p)
