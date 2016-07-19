# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 14:04:10 2016

@author: mmayers
"""

import pandas as pd
import numpy as np
from uuid import uuid4
from io import StringIO


def read_file(fname):
    with open(fname, "r") as fin:
        for line in fin:
            yield line.rstrip("\n")


def get_id():
    return str(uuid4()).split("-")[0]


def read_header(fname):
    for line in read_file(fname):
        if not line.startswith("H"):
            break
        yield line.split('\t')[1:]


def get_colnames(fname):
    # parse out the dataframe column headers
    colnames = dict()
    for head in read_header(fname):
        if head[0] in ["PLINE", "SLINE", "&SLINE"] and len(head) > 2:
            colnames[head[0]] = ["type"] + [h.lower() for h in head[1:]]
            if head[0] == "SLINE":
                colnames["SLINE"] += ["enrichment_mr"]

    return colnames


def get_corrval(fname):
    import math
    for head in read_header(fname):
        if head[0] == 'Correction Factor (Ln)':
            return math.exp(float(head[1]))
    return 1

def parse_file(fname, norm = True):
    """Read through the original file, split into P and S files"""

    import math
    
    colnames = get_colnames(fname)

    p_out, s_out, sand_out = [], [], []
    p_id = get_id()
    make_new_id = False
    for line in read_file(fname):
        if not line.startswith("H"):
            if line.startswith("P"):
                if make_new_id:
                    p_id = get_id()
                    make_new_id = False

                out = p_out
            else:
                # S or &S
                make_new_id = True

                if line.startswith("S"):
                    vals = line.split("\t")
                    line = "\t".join(vals[:-6] + vals[-2:])

                out = sand_out if line.startswith("&") else s_out

            out.append("{}\t{}\n".format(p_id, line))


    # read into pandas
    p_table = (pd
        .read_csv(StringIO(''.join(p_out)), sep = '\t', names = colnames["PLINE"])
        .dropna(axis = 1, how = "all")
        .reset_index()
    )

    s_table = (pd
        .read_csv(StringIO(''.join(s_out)), sep = '\t', names = colnames["SLINE"])
        .reset_index()
    )
    sand_table = (pd
        .read_csv(StringIO(''.join(sand_out)), sep = '\t', names = colnames["&SLINE"])
        .reset_index()
    )

    s_table = pd.concat([s_table, sand_table]).reset_index(drop=True)

    # Add a column 'aa_sequence' without the leading and trailing Residues
    s_table['aa_sequence'] = s_table['sequence'].str.split('.').apply(lambda x: '.'.join(x[1:-1]))

    # Add a boolean to show the source of the protein
    s_table['n15'] = s_table['file_name'].str.contains('n15')
    
    # Get the correction Factor (-median of the log of the peptide ratios)
    corr_factor = -1*(np.log(s_table
                             .drop_duplicates(subset=['cs', 'aa_sequence'])['ratio']    
                             .replace(0, np.nan))
                      .dropna()
                      .median())
    
    # Multiply by the correction factor to get the new 
    s_table['sam_int_corr'] = (s_table['sam_int'] * math.exp(corr_factor)).astype(int)


    #Normilazaion routine
    # ratios are multiplied by math.exp(corr_factor)
    # reverse ratios are divided by math.exp(corr_factor)
    if norm:    
        s_ratios = ['ratio', 'area_ratio']
        s_rev_ratios = ['rev_slope_ratio']
        
        s_table[s_ratios] = s_table[s_ratios].apply(lambda c: c * math.exp(corr_factor))
        s_table[s_rev_ratios] = s_table[s_rev_ratios].apply(lambda c: c / math.exp(corr_factor))
        
        p_ratios = ['average_ratio', 'standard_deviation', 'composite_ratio', 'weighted_average', 'log_inv_average']
        p_rev_ratios = ['average_ratio_rev', 'standard_deviation_rev', 'log_inv_average_rev']

        p_table[p_ratios] = p_table[p_ratios].apply(lambda c: c * math.exp(corr_factor))        
        p_table[p_rev_ratios] = p_table[p_rev_ratios].apply(lambda c: c / math.exp(corr_factor))
        
    
    for col in ['sam_int', 'ref_int', 'peak_int']:
        s_table[col] = s_table[col].astype(int)

    for integration in ['sam_int_corr', 'ref_int', 'peak_int']:
        s_table['log_{}'.format(integration)] = np.log2(s_table[integration])            


    return (p_table, s_table)