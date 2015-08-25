from collections import defaultdict
import pandas as pd
import parse
import os
import re
import numpy as np

RE = re.compile(r'(?P<segid1>\w+)[(](?P<resid1>\d+)[)]-(?P<segid2>\w+)[(](?P<resid2>\d+)[)]')        

def loadfasta(fasta, aa):
    '''
    input: name of fasta file
    create map of ordered amino acids of interest
    '''
    prot2map = {}
    for name,seq in parse.parse_fasta(fasta):
        prot2map[name] = defaultdict(int)
        prot2map[name].update({i:i for i,c in enumerate(seq) if c in aa})
    return prot2map

def load_fasta(fasta):
    prot_map = {}
    for name,seq in parse.parse_fasta(fasta):
        protein = []
        for i,c in enumerate(seq):
            protein.append({'aminoacid': c, 'pos': i+1})
        prot_map[name] = pd.DataFrame.from_records(protein)
    return prot_map


def load_plink_html(filename):
    tp = parse.PLinkParser()
    with open(filename, 'r') as fi:
        for line in fi:
            tp.feed(line)

    def build_xl_df(df, xl_field):
        def extract_xl_coord(rec):
            m = RE.search(rec[xl_field])
            return pd.Series([
                m.group('segid1'), int(m.group('resid1')), m.group('segid2'), int(m.group('resid2'))
                ], index=['prot1', 'res1', 'prot2', 'res2'])

        return df.merge(df.apply(extract_xl_coord, axis=1), left_index=True, right_index=True)
    
    df_details = pd.DataFrame.from_records(tp.details_records).dropna()
    df_sum = pd.DataFrame.from_records(tp.sum_records).dropna()
    return build_xl_df(df_details, 'Proteins'),  build_xl_df(df_sum, 'ProteinAC')


def loadplink(directory, prot2map, allprot, scale, cutoff=1):
    '''
    Input: directory full of plink txt files, the prot2map from loadfasta,
    and the dictionary of all involved proteins thus far creates a dictionary of 
    gene-gene to coordinate tuples and updates allprot
    '''
    gg2i = {} #gene-gene to interaction coordinates

    filenames = [fname for fname in os.listdir(directory) if fname.endswith('.txt')]
    for filename in filenames:
        with open(os.path.join(directory, filename), 'r') as fi:
            for line in fi:
                if (',' in line) and ')-' in line and not ('REVERSE' in line):
                    tokens = (line.rstrip()).split()
                    interaction = tokens[-1]
                    # optional e value cutoff with 2.0E-04 as example
                    e = float(tokens[5])
                    if e < cutoff:
                        subtokens = re.split(r'[)|(|-]',interaction)
                    protloc = [(subtokens[0],subtokens[1]),(subtokens[3],subtokens[4])]
                    protloc = sorted(protloc)
                    prot1 = protloc[0][0]
                    prot2 = protloc[1][0]
                    if prot1 in prot2map and prot2 in prot2map:
                        loc1 = prot2map[prot1][int(protloc[0][1])]
                        loc2 = prot2map[prot2][int(protloc[1][1])]
                        if scale:
                            loc1 = int(protloc[0][1])
                            loc2 = int(protloc[1][1])
                        key = prot1 + '-' + prot2
                        if not key in gg2i:
                            gg2i[key] = defaultdict(int)
                        gg2i[key][(loc1,loc2)] += 1

                        #Now handle symmetry
                        key = prot2 + '-' + prot1
                        if not key in gg2i:
                            gg2i[key] = defaultdict(int)
                        gg2i[key][(loc2,loc1)] += 1

                        #Add protein names to allprot for graphing purposes
                        allprot[prot1]=max(prot2map[prot1])
                        allprot[prot2]=max(prot2map[prot2])
    return (gg2i,allprot)


def writesummary(key, gg2i, prot2map, scale):
    (prot1, prot2) = key.split('-')
    x = np.empty([0])
    y = np.empty([0])
    r = np.empty([0])
    mc = (max(prot2map[prot1].values()),max(prot2map[prot2].values()))
    if scale:
        mc = (max(prot2map[prot1]),max(prot2map[prot2]))
    if key in gg2i:
        for t in gg2i[key]:
            x=np.append(x,t[0])
            y=np.append(y,t[1])
            r=np.append(r,gg2i[key][t])
    return (x, y, r, mc)
