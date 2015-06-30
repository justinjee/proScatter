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
from bokeh.plotting import figure
from bokeh.io import gridplot, output_file, show

#handle user options
args = []
args2 = sys.argv[:]
zoom=False
scale=True #sets default to scaling based on aa only
for i in range(len(args2)):
    a = args2[i]
    if '--zoom' in a:
        zoom=True
        (z,zkey)=a.split('=')
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

#Define basic plot parameters
color = ['blue','red','green']
marker = ['o','*','+']
basesize = 25 
buffer=10
numprot = len(allprot)
sallprot = sorted(allprot)
handles = [0]*len(interactions)

#Now do actual plotting
output_file(output+'.html')
if not zoom:
  m = [[None for r in range(numprot)] for s in range(numprot)]
  i=numprot-1
  for prot2 in sallprot:
      j=0
      for prot1 in sallprot:
          key = prot1 +'-'+ prot2
          for xind in range(len(interactions)):
              gg2i=interactions[xind]
              (x,y,r,mc) = loadfiles.writesummary(key,gg2i,prot2map,scale)
              m[i][j]=splotch.splotch(m[i][j],x,y,r,basesize,buffer,mc,key,i==numprot-1,j==0,numprot,color[xind],marker[xind],dirs[xind],i==0 and j==numprot-1)
          j+=1
      i-=1
  p = gridplot(m)
else:
    p = None
    for xind in range(len(interactions)):
        gg2i=interactions[xind]
        (x,y,r,mc) = loadfiles.writesummary(zkey,gg2i,prot2map,scale)
        p = splotch.splotch(p,x,y,r,basesize,buffer,mc,zkey,True,True,1,color[xind],marker[xind],dirs[xind],True)

#f.legend(tuple(handles),tuple(dirs),'upper right')
show(p)
