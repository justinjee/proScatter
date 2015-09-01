from collections import defaultdict
import pandas as pd
import parse
import os
import re
import numpy as np

RE = re.compile(r'(?P<segid1>\w+)[(](?P<resid1>\d+)[)]-(?P<segid2>\w+)[(](?P<resid2>\d+)[)]')        

def load_aa_positions(fasta, aa):
    '''
    For every protein sequence in `fasta` builds an index of aminoacids specified in `aa`
    Returns long DataFrame
    '''

    def build_aa_df(name, seq, aa):
        res = []
        for i,c in enumerate(seq):
            if c in aa:
                res.append({'prot': name, 'aa': c, 'pos': i+1})
        return pd.DataFrame.from_records(res)
    
    aa_map = []
    for name,seq in parse.parse_fasta(fasta):
        aa_map.append(build_aa_df(name, seq, aa))
    return pd.concat(aa_map)


def sym_df(df, columns={'res1': 'res2', 'res2': 'res1', 'prot1':'prot2', 'prot2':'prot1',}):
    df_ = df.copy()
    df_.rename(columns=columns, inplace=True)
    df_ = pd.concat([df, df_])
    df_.reset_index(drop=True, inplace=True)
    return df_


def load_plink_html(filename, explabel=None):
    tp = parse.PLinkParser()
    with open(filename, 'r') as fi:
        for line in fi:
            tp.feed(line)

    def build_xl_df(df, xl_field):
        def extract_xl_coord(rec):
            m = RE.search(rec[xl_field])
            return pd.Series([
                m.group('segid1'), int(m.group('resid1')), m.group('segid2'), int(m.group('resid2')), explabel
                ], index=['prot1', 'res1', 'prot2', 'res2', 'label'])

        return df.merge(df.apply(extract_xl_coord, axis=1), left_index=True, right_index=True)
    
    df_details = build_xl_df(pd.DataFrame.from_records(tp.details_records).dropna(), 'Proteins')
    df_sum = build_xl_df(pd.DataFrame.from_records(tp.sum_records).dropna(), 'ProteinAC')
    return  sym_df(df_details), sym_df(df_sum)
