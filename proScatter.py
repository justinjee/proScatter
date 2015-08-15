#!/usr/bin/env python

#proScatter
#Interaction visualizer for pLink data
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

class ProScatter(object):

    color = ['blue','red','green']
    marker = ['o','o','o']
    basesize = 50 
    buffer = 10
    
    def __init__(self, **kwargs):
        super(ProScatter, self).__init__()
        output = kwargs.pop('output')
        if output.endswith('.html'):
            self.output = output
        else:
            self.output = output + '.html'
        for key,value in kwargs.iteritems():
            setattr(self, key, value)
        self._interactions = None
        self._prot2map = None
        self._m = None
        self._allprot = None
        self._numprot = 0

    def _reverse_size_sort(self, proteins):
        sprot = sorted(proteins.items(), key=lambda t:t[1])
        return [k[0] for k in sprot][::-1]

    def load_data(self):
        self._prot2map = loadfiles.loadfasta(self.fasta_file, self.aminoacids)
        allprot = defaultdict(int)
        interactions = []
        for dirname in self.data_dirs:
            gg2i, allprot = loadfiles.loadplink(dirname, self._prot2map, allprot, self.scale, self.evalue)
            interactions.append(gg2i)
        self._interactions = interactions
        self._numprot = len(allprot)
        self._allprot = self._reverse_size_sort(allprot)

    def build_plot(self):
        if self.zoom is None:
            m = [[None for r in range(self._numprot+1)] for s in range(self._numprot)]
            for xind,interaction in enumerate(self._interactions):
                m[0][-1] = splotch.makelegend(m[0][-1],
                        self._numprot,
                        self.color[xind],
                        self.marker[xind],
                        self.data_dirs[xind])
                i =  self._numprot - 1
                for prot2 in self._allprot:
                    j = 0
                    for prot1 in self._allprot:
                        key = prot1 + '-' + prot2
                        (x, y, r, mc) = loadfiles.writesummary(key, interaction, self._prot2map, self.scale)
                        m[i][j] = splotch.splotch(m[i][j], m[self._numprot-1][0], x, y, r,
                                self.basesize,
                                self.buffer,
                                mc, key, i==self._numprot-1, j==0,
                                self._numprot,
                                self.color[xind],
                                self.marker[xind], None, False, self.unjoin)
                        j+=1
                    i-=1
            self.proscatter = gridplot(m)
        else:
            for xind,interaction in enumerate(self._interactions):
                (x, y, r, mc) = loadfiles.writesummary(self.zoom, interaction,
                        self._prot2map, self.scale)
                self.proscatter = splotch.splotch(None, None, x, y, r,
                        self.basesize,
                        self.buffer, mc, self.zoom, True, True, 1,
                        self.color[xind],
                        self.marker[xind],
                        self.data_dirs[xind], True)

    def show_scatter(self):
        if self.proscatter is not None:
            output_file(self.output)
            show(self.proscatter)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''
    pScatter
    Interaction visualizer for pLink XLMS data, by Justin Jee (with design by Katelyn McGary Shipper)

    '''
    )
    parser.add_argument('fasta_file', type=str, help='fasta file with protein sequences')
    parser.add_argument('data_dirs', nargs='+', help='pLink output directory')
    parser.add_argument('-a', '--aminoacids', default='K', help='cross-linkable aminoacids. Defaults to Lysine (K).')
    parser.add_argument('-z', '--zoom', help='Prot1-Prot2 only display subplot for proteins Prot1 vs Prot2',
                        type=str)
    parser.add_argument('-s', '--scale', help='scale both plot and outputs so that only amino acids of interest are considered', action='store_false')
    parser.add_argument('-e', '--evalue', default=1.0, help='e-value cutoff', type=float)
    parser.add_argument('-u', '--unjoin', help='unjoin plot axes', action='store_true')
    parser.add_argument('-o', '--output', default='output.html', help='output file (HTML) name')
    parser.add_argument('-v', '--verbose', help='increase output verbosity', action='store_true')
    args = parser.parse_args()

    kwargs = vars(args)
    proscatter = ProScatter(**kwargs)
    proscatter.load_data()
    proscatter.build_plot()
    proscatter.show_scatter() 
