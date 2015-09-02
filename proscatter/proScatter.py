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
import pandas as pd
import os
from bokeh.plotting import figure
from bokeh.io import gridplot, output_file, show
from bokeh.models import Range1d

class ProScatter(object):

    def __init__(self, **kwargs):
        super(ProScatter, self).__init__()
        output = kwargs.pop('output')
        if output.endswith('.html'):
            self.output = output
        else:
            self.output = output + '.html'
        for key,value in kwargs.iteritems():
            setattr(self, key, value)
        self._df_sum = pd.DataFrame()
        self._df_details = pd.DataFrame()
        self._df_fasta = None

    def load_data(self):
        for folder in self.plink:
            for htmlfile in os.listdir(folder):
                df_details, df_sum = loadfiles.load_plink_html(folder+'/'+htmlfile, folder)
                self._df_details = pd.concat([self._df_details, df_details])
                self._df_details.reset_index(drop=True, inplace=True)
                self._df_sum = pd.concat([self._df_sum, df_sum])
                self._df_sum.reset_index(drop=True, inplace=True)
        self._df_fasta = loadfiles.load_aa_positions(self.fasta_file, self.aminoacids)

    def print_summary(self):
        grouped = self._df_sum.groupby(['prot1', 'prot2'])
        for gname,group in grouped:
            print("Summary for {}".format(gname))
            print("{} interactions".format(group.shape[0]))

    def build_plot(self):

        prot_list = set(self._df_sum['prot1'])
        sorted_prot_list = sorted(prot_list.intersection(self._df_fasta['prot']))

        numprot = len(sorted_prot_list)
        rows = [[None for r in range(numprot+1)] for s in range(numprot)]
        rows[0][-1]=splotch.makelegend(self._df_sum)
        (xr, yr) = (Range1d(0,self._df_sum['res1'].max()), Range1d(0,self._df_sum['res2'].max()))
        i=numprot-1
        for prot2 in sorted_prot_list:
            j=0
            for prot1 in sorted_prot_list:
                subdf = self._df_sum[(self._df_sum['prot1']==prot1) & (self._df_sum['prot2']==prot2)]
                if self.unjoin:
                    (xr, yr) = ((0,self._df_fasta[self._df_fasta['prot']==prot1]['pos'].max()), 
                                (0,self._df_fasta[self._df_fasta['prot']==prot2]['pos'].max()))
                rows[i][j] = splotch.splotch_df(subdf, prot1, prot2, xr, yr, numprot, i==numprot-1, j==0)
                j+=1
            i-=1

        self.proscatter = gridplot(rows)

    
    def show_scatter(self):
        if self.proscatter is not None:
            output_file(self.output)
            show(self.proscatter)

def main():
    parser = argparse.ArgumentParser(description='''
    ProScatter
    Interaction visualizer for pLink XLMS data, by Justin Jee (with design by Katelyn McGary Shipper)

    '''
    )
    parser.add_argument('fasta_file', type=str, help='fasta file with protein sequences')
    parser.add_argument('-p', '--plink', nargs='+', type=str, help='pLink output .html file')
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

if __name__ == "__main__":
    main()
