# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 15:47:55 2016

@author: mmayers
"""

import sys
sys.path.append('..')
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns

from scripts.analysis import build_loci


def clean_df(df):
    # You lose about 90% of peptides if you just drop empty ones.
    # Set them to 1 and then do the log tranform
    df.replace(np.nan, 0, inplace=True)
    df = np.log(df+1)
    return df
    


def plot_pca(df, metadata, group_name = 'sample_type', loc=3, title='PCA of Protein Clusters', out = 'pca', save = False, label_points = False):
    """
    Needs dataframe with samples as row and features as columns and corresponding
    metadata with samples as row and column with categorical variable for grouping.
    Groups is the categorical variable to group by for analysis.
    """    
        
    sns.set_style('white')
   
    from sklearn import preprocessing
    from sklearn.decomposition import PCA

    BASE = os.getcwd()
    GROUP = group_name
    PLOT_PCS = tuple([1, 2])        
    NUM_PCS = max(PLOT_PCS)

    metadata = metadata.T
    df = clean_df(df)
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
        samples = list(df.loc[metadata.loc[GROUP]==group].index)
        data = X_r[X_r.index.isin(samples)]
        ax.scatter(data[PLOT_PCS[0]], data[PLOT_PCS[1]], c=[colors[i]]*len(df), label=group, s=50)
        if label_points:
            for i, samp in enumerate(samples):
                ax.annotate(i+1, xy = (data[PLOT_PCS[0]][i], data[PLOT_PCS[1]][i]), xytext = (0, 10),
                    textcoords = 'offset points', ha = 'right', va = 'bottom',
                    arrowprops = dict(arrowstyle = '-', connectionstyle = 'arc3,rad=0'))

    plt.legend(loc=loc,title=GROUP)
    plt.title(title, size=20)
    plt.xlabel("PC{0}: {1:.2f}%".format(PLOT_PCS[0], pca.explained_variance_ratio_[PLOT_PCS[0]-1]*100), size=16)
    plt.ylabel("PC{0}: {1:.2f}%".format(PLOT_PCS[1], pca.explained_variance_ratio_[PLOT_PCS[1]-1]*100), size=16)

    fig.set_tight_layout(False)
    if save:
        fig.savefig(os.path.join(BASE, out + ".png"))
        fig.savefig(os.path.join(BASE, out + ".pdf"))


def plot_volcano(grouped_loci, fc_cutoff=4, p_val_cutoff=.05, labels=False, title='Volcano Plot'):
    sns.set_style('whitegrid')    

    # This code should probably be somewhere else
    def get_gene_name(r):
        if r['name'].find('GN=') > 0:
            return r['name'][r['name'].index('GN=')+3:].split(' ')[0]
        else:
            return ''
    
    df = build_loci.get_sig_df(grouped_loci, fc_cutoff=fc_cutoff, p_val_cutoff=p_val_cutoff)
    df['gene_name'] = df.apply(get_gene_name, axis=1)
    
    X = df['logfc'].values
    Y = df['-logp'].values
    
    # Plot the points
    fig = plt.figure(figsize=(15, 9))
    fig.set_tight_layout(False)
    plt.scatter(X, Y, c = df['passes'], s=50, cmap='bwr')
    
    if labels:
        for name, x, y in zip(df['gene_name'],X, Y):
            if abs(x) > np.log2(fc_cutoff) and y > -1*np.log10(p_val_cutoff):
                plt.annotate(name, xy = (x, y), xytext = (5, 5), textcoords = 'offset points', arrowprops = dict(arrowstyle = '-', connectionstyle = 'arc3,rad=0'))

    
    # Plot asjustments
    plt.title(title, size=24)
    plt.xlabel('log2 Fold Change', size=18)
    plt.ylabel('-log10 p-value', size=18)
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.ylim(bottom=-0.2);