BASE = "../data"

import sys
sys.path.append('..')

from scripts import utils
import os
import pandas as pd
from itertools import chain
import shelve
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from sklearn import preprocessing
from collections import Counter
from collections import defaultdict

from scripts.analysis import build_loci
from scripts.analysis.DBInfo import DBInfo

db_info = DBInfo("compil_mgm")
metadata = build_loci.read_metadata(os.path.join(BASE,"metadata.csv"))

#%% Parse samples
samples = shelve.open(os.path.join(BASE,"samples.shelve"))
already_done_samples = set(list(samples.keys()))
for sample_name, sample_info in tqdm(list(metadata.iteritems())):
    if sample_name not in already_done_samples:
        sample = build_loci.Sample(sample_name, sample_info.path, db_info, sample_info)
        samples[sample.sample_name] = sample
#%%
protein_clusters = shelve.open(os.path.join(BASE,"protein_clusters.shelve"))
for name,sample in samples.items():
    protein_clusters[name] = sample.build_protein_clusters()

"""
grouped_loci = build_loci.group_across_samples(list(chain(*protein_clusters.values())), db_info)
for locus in tqdm(grouped_loci):
    locus.annotate()
utils.save(grouped_loci, os.path.join(BASE,"grouped_loci.pkl.gz"))
#grouped_loci = utils.load(os.path.join(BASE,"grouped_loci.pkl.gz"))

#%% filtering
grouped_loci = [x for x in grouped_loci if x.passes_thresh(metadata, min_samples_per_group=0)]
utils.save(grouped_loci, os.path.join(BASE,"grouped_loci_filt.pkl.gz"))
build_loci.to_json(grouped_loci, list(samples.keys()), "/mongoc/mmayers/N15-Samples/N15_Samples_.json")
#%% Normalize
nf = build_loci.yates_normalization(samples)
for locus in grouped_loci:
    locus.normalize(nf)
utils.save(grouped_loci, os.path.join(BASE,"grouped_loci_filt_norm.pkl.gz"))
#grouped_loci = utils.load(os.path.join(BASE,"grouped_loci_filt_norm.pkl.gz"))
#%%
df = build_loci.to_df(grouped_loci)
# keep only some samples
samples_keep = {'RTcoh_unenr1','RTcoh_unenr2','RTcoh_unenr3','RTsep_unenr1','RTsep_unenr2','RTsep_unenr3','Ragsep_unenr1','Ragsep_unenr2','Ragsep_unenr3','WTsep_unenr1','WTsep_unenr2','WTsep_unenr3'}
df = df.iloc[df.index.isin(samples_keep)]

# Log transform
df = np.log(1+df)
# Scale
X = preprocessing.scale(df)
"""
