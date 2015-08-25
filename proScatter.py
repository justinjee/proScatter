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

from __future__ import division, print_function
import loadfiles
import math
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
        self._df_sum = None
        self._df_details = None

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

    def load_html(self):
        self._df_sum, self._df_details = loadfiles.load_plink_html(self.html_file)

    def print_summary(self):
        for prot2 in self._allprot:
            for prot1 in self._allprot:
                key = prot1 + '-' + prot2
                print("Summary for {}".format(key))
                for xind,interaction in enumerate(self._interactions):
                    points = set()
                    if key in interaction:
                        points = set(interaction[key])
                    print("\t{0}: {1}".format(splotch.stripfolder(self.data_dirs[xind]),len(points)))
                    if xind==0:
                        commonpoints = points
                    else:
                        commonpoints = commonpoints.intersection(points)
                print("\tIntersection: {}".format(len(commonpoints)))

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

    def make_plot(self):
        TOOLS = "resize,crosshair,pan,wheel_zoom,box_zoom,reset,previewsave"
        rows = []
        row = []
        num_prot = len(set(self._df_sum['prot1']))
        grouped = self._df_sum.groupby(['prot1', 'prot2'])

        for i, (key,group) in enumerate(grouped):
            fig = figure(tools=TOOLS, width=250, height=250, title=None, min_border=10)
            scatter = fig.scatter(group['res2'], group['res1'],
                        size=np.sqrt(self.basesize / num_prot * group['SpecNum']),
                        alpha=0.25)
            fig.xaxis.axis_label = key[1]
            fig.yaxis.axis_label = key[0]
            fig.xaxis.major_label_orientation = math.pi / 4
            fig.xaxis.visible = False
            fig.yaxis.visible = False
            row.append(fig)
            if (i + 1) % num_prot == 0:
                fig.yaxis.visible = True
                rows.append(row[::-1])
                row = []

        for fig in rows[-1]:
            fig.xaxis.visible = True

        self.proscatter = gridplot(rows)

    
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
    parser.add_argument('-h', '--html', type=str, help='pLink output in .html format')
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
    if proscatter.html:
        proscatter.load_html()
    proscatter.print_summary()
    proscatter.build_plot()
    proscatter.show_scatter() 
