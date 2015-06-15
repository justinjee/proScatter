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
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from pylab import *

#handle user options
args = []
args2 = sys.argv[:]
zoom=False
xkcd=False
scale=True #sets default to scaling based on aa only
for i in range(len(args2)):
    a = args2[i]
    if '--zoom' in a:
        zoom=True
        (z,zkey)=a.split('=')
    elif '--xkcd' in a:
        xkcd=True
    elif '--scale' in a:
        scale=False
    else:
        args.append(a)

fasta = args[1]
aa = args[2]
dirs = args[3:-1]
output = args[-1]

### load fasta file ###
prot2map = loadfiles.loadfasta(fasta,aa)

### load plink file ###
print "fasta loaded. loading pLink files"
allprot = defaultdict(int)
interactions = []
for dir in dirs:
    gg2i,allprot = loadfiles.loadplink(dir,prot2map,allprot,scale)
    interactions.append(gg2i)

### output ###
print "generating plot"
o = open(output+'.txt','w')

#Define basic plot parameters
if xkcd:
    plt.xkcd()
color = ['b','r','g']
marker = ['o','x','*']
basesize = 80
buffer=10
numprot = len(allprot)
sallprot = sorted(allprot)

#Now do actual plotting
j=1
if not zoom:
  f, ax = plt.subplots(numprot,numprot,sharex=True,sharey=True)
  for prot1 in sallprot:
    i=numprot-1
    for prot2 in sallprot:
        key = prot1 +'-'+ prot2
        ax = subplot(numprot,numprot,i*numprot+j)
        for xind in range(len(interactions)):
            gg2i=interactions[xind]
            (x,y,r,mc) = loadfiles.writesummary(o,key,gg2i,prot2map,scale)
            splotch.splotch(ax,x,y,r,basesize,buffer,mc,key,j==1,i==numprot-1,numprot,color[xind],marker[xind])
        i-=1
    j+=1
else:
    f, ax = plt.subplots(1,1)
    for xind in range(len(interactions)):
        gg2i=interactions[xind]
        (x,y,r,mc) = loadfiles.writesummary(o,zkey,gg2i,prot2map,scale)
        splotch.splotch(ax,x,y,r,basesize,buffer,mc,zkey,True,True,1,color[xind],marker[xind])

def onpick(event):
    ind = event.ind
    (prot1,prot2)=zkey.split('-')
    print prot1+"="+str(np.take(x, ind))+ " "+prot2+"="+str(np.take(y, ind))
plt.tight_layout(pad=1, w_pad=0, h_pad=0)
f.savefig(output+'.png')
if zoom:
    print "Click features enabled."
    f.canvas.mpl_connect('pick_event', onpick)
plt.show()
