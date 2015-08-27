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
import argparse
import numpy as np
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
        self._df_sum = None
        self._df_details = None
        self._df_fasta = None

    def load_data(self):
        self._df_details, self._df_sum = loadfiles.load_plink_html(self.plink)
        self._df_fasta = loadfiles.load_aa_positions(self.fasta_file, self.aminoacids)

    def print_summary(self):
        grouped = self._df_sum.groupby(['prot1', 'prot2'])
        for gname,group in grouped:
            print("Summary for {}".format(gname))
            print("{} interactions".format(group.shape[0]))

    def build_plot(self):
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
    ProScatter
    Interaction visualizer for pLink XLMS data, by Justin Jee (with design by Katelyn McGary Shipper)

    '''
    )
    parser.add_argument('fasta_file', type=str, help='fasta file with protein sequences')
    parser.add_argument('plink', type=str, help='pLink output .html file')
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
    proscatter.print_summary()
    proscatter.build_plot()
    proscatter.show_scatter() 
