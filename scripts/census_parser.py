# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 14:04:10 2016

@author: b34k
"""

import os
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

def parse_file(fname, heavy = False, norm = False):
    """Read through the original file, split into P and S files"""

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

    s_table = pd.concat([s_table, sand_table])

    if heavy:
        s_table['source'] = 'n15'
    else:
        s_table['source'] = 'n14'

    s_table['sam_int_corr'] = (s_table['sam_int'] * get_corrval(fname)).astype(int)
    for col in ['sam_int', 'ref_int', 'peak_int']:
        s_table[col] = s_table[col].astype(int)

    for integration in ['sam_int_corr', 'ref_int', 'peak_int']:
        if norm:
            s_table['log_{}'.format(integration)] = norm_integrations(s_table[integration])            
        else:
            s_table['log_{}'.format(integration)] = np.log2(s_table[integration])            

    return (p_table, s_table)


def norm_integrations(col):
    mu = np.log2(col).replace(np.log2(0), np.nan).dropna().mean()
    sig = np.log2(col).replace(np.log2(0), np.nan).dropna().std()

    log_col = np.log2(col).replace(np.log2(0), np.nan)

    log_col = (log_col - mu) / sig
    return log_col

def join_ptables(light, heavy):

    keys = ['locus', 'description']
    return light.merge(heavy, on=keys, how='outer', suffixes = ['_l', '_h'])


def join_stables(light, heavy):

    return pd.concat([light, heavy])


def parse_l_h_together(light_file, heavy_file, join_s = False, only_s = False, norm = False):
    '''
    Takes census output files generated from light and heavy search on same
    sample and retruns two dataframes, one for proteins, one for peptides.
    Df's are linked by index_l, index_h in protein df to index in the peptide
    df.
    '''

    prot_l, pep_l = parse_file(light_file, norm = norm)
    prot_h, pep_h = parse_file(heavy_file, heavy = True, norm = norm)

    if join_s:
        return (join_ptables(prot_l, prot_h), join_stables(pep_l, pep_h))
    elif only_s:
        return (pep_l, pep_h)
    else:
        return (join_ptables(prot_l, prot_h), join_stables(pep_l, pep_h))


def prep_for_pca(df_list, name = 'Unenr', clean = True):
    """
    Takes a list of paired peptide dfs ex:
    [[df_light1, df_heavy1], [df_light2, df_heavy2], [df_light3, df_heavy3]]
    Can add a name that can be used for groups in PCA
    returns a single DF with peptides each column as a peptide,
    each row as a sample, and the values are the peak integrations
    """


    dats = []

    for i, dfs in enumerate(df_list):
        dats.append((dfs[0]
        .groupby(['sequence', 'cs'])
        .mean()[['log_sam_int_corr', 'log_ref_int']]
        .combine_first(dfs[1].groupby(['sequence', 'cs'])
        .mean()[['log_sam_int_corr', 'log_ref_int']])
        ).rename(columns={'log_sam_int_corr':'{}_n14_{}'.format(name, i+1), 'log_ref_int': '{}_n15_{}'.format(name, i+1)})
        .T)

    if clean:
        return clean_pca_df(pd.concat(dats))
    return pd.concat(dats)



def prep_for_pca2(df_list, name = 'Unenr', clean = False):
    """
    Returns same DF as prep for pca, however
    for Singletons, instead of taking the integration
    value for both peptides, it only takes it for the
    relevant peptides, and places a 0 (log(1)) in the
    for the integration of the other peptide.
    """


    dats = []

    for i, dfs in enumerate(df_list):
        dats.append(((dfs[0].groupby(['type']).get_group('S'))
        .groupby(['sequence', 'cs'])
        .mean()[['log_sam_int_corr', 'log_ref_int']]
        .combine_first(dfs[1]
            .groupby(['type']).get_group('S')
            .groupby(['sequence', 'cs'])
            .mean()[['log_sam_int_corr', 'log_ref_int']])
        .combine_first(pd.DataFrame(dfs[0]
                .groupby(['type']).get_group('&S')
                .groupby(['sequence', 'cs'])
                .mean()['log_sam_int_corr'])
                      )
        .combine_first(pd.DataFrame(dfs[1]
                .groupby(['type']).get_group('&S')
                .groupby(['sequence','cs'])
                .mean()['log_ref_int'])
                      )
        )
        .rename(columns={'log_sam_int_corr':'{}_n14_{}'.format(name, i+1), 'log_ref_int': '{}_n15_{}'.format(name, i+1)})
        .T
        .fillna(0))

    if clean:
        return clean_pca_df(pd.concat(dats))
    return pd.concat(dats)


def plot_pca(df, loc = 'best', labels = False, save = False, out = 'pca_plot'):
    """
    Given a df with samples
    """


    from sklearn import preprocessing
    from sklearn.decomposition import PCA
    import matplotlib.pyplot as plt
    BASE = os.getcwd()

    metadata = pd.DataFrame(['_'.join(s.split('_')[:2]) for s in list(df.index)], columns = ['sample_type'])
    metadata.index = list(df.T.columns)
    metadata = metadata.T


    GROUP = 'sample_type'
    PLOT_PCS = tuple([1, 2])
    NUM_PCS = max(PLOT_PCS)
    X_scaled = pd.DataFrame(preprocessing.scale(df), index=df.index, columns = df.columns)
    pca = PCA(n_components=NUM_PCS)
    X_r = pd.DataFrame(pca.fit_transform(X_scaled), index=df.index, columns = list(range(1,NUM_PCS+1)))

    print('explained variance ratios: %s' % str(pca.explained_variance_ratio_))
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    groups = np.unique(list(metadata.loc[GROUP]))
    num_colors = len(groups)
    cm = plt.get_cmap("Paired")
    colors = [cm(1.*i/num_colors) for i in range(num_colors)]

    for i,group in enumerate(groups):
        samples = list(metadata.columns[metadata.loc[GROUP]==group])
        data = X_r[X_r.index.isin(samples)]
        ax.scatter(data[PLOT_PCS[0]], data[PLOT_PCS[1]], c=[colors[i]]*len(df), label=group, s=50)
        if labels:
            for i, samp in enumerate(samples):
                ax.annotate(samp.split('_')[-1], xy = (data[PLOT_PCS[0]][i], data[PLOT_PCS[1]][i]), xytext = (0, 10),
                    textcoords = 'offset points', ha = 'right', va = 'bottom',
                    arrowprops = dict(arrowstyle = '-', connectionstyle = 'arc3,rad=0'))

    plt.legend(loc=loc,title=GROUP)
    plt.title('PCA of Protein Clusters')
    plt.xlabel("PC{0}: {1:.2f}%".format(PLOT_PCS[0], pca.explained_variance_ratio_[PLOT_PCS[0]-1]*100))
    plt.ylabel("PC{0}: {1:.2f}%".format(PLOT_PCS[1], pca.explained_variance_ratio_[PLOT_PCS[1]-1]*100))

    fig.set_tight_layout(False)
    if save:
        fig.savefig(os.path.join(BASE, out + ".png"))
        fig.savefig(os.path.join(BASE, out + ".pdf"))

def clean_pca_df(df, verbose = True):
    """
    Removes peptides from the DF that are not quantified in ALL samples
    """
    if verbose:
        print('Starting peptides: ', len(df.T))
    df = df.T.reset_index().T[2:].replace(np.log2(0), np.nan).dropna(axis = 1)
    if verbose:
        print('Peptides after Cleaning: ', len(df.T))
    return df
