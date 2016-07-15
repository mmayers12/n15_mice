import os
from itertools import chain
from collections import defaultdict
from pymongo import MongoClient
import pandas as pd
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from tqdm import tqdm
rcParams.update({'figure.autolayout': True})



if False:
    #%%  test
    from ... import utils
    import os
    import pandas as pd
    from itertools import chain
    import shelve
    from collections import Counter
    from .. import build_loci
    from .. import DBInfo
    BASE = "/home/gstupp/projects/Wolan/wolan/UC12_02042016"
    db_info = DBInfo("compil_mgm")
    metadata = build_loci.read_metadata(os.path.join(BASE,"metadata.csv"))
    grouped_loci = utils.load(os.path.join(BASE,"grouped_loci.pkl.gz"))

#%%


def annotate(grouped_loci, db_info):
    # Grouping together some mongodb queries to make this faster. 10-20x faster
    all_cluster_prot_ids = set(chain(*[locus.cluster_prot_ids for locus in grouped_loci]))
    prot_result = {x['_id']:x for x in db_info.protDB.find({'_id':{'$in':list(all_cluster_prot_ids)}})}
    taxIDs_doc = {x['_id']: x for x in db_info.taxDB.find({'_id':{'$in':list(all_cluster_prot_ids)}},{'_id':True,'taxid': True})}
    for locus in tqdm(grouped_loci):
        locus.lookup_name(prot_result)
        locus.lookup_function()
        locus.get_go()
        tax_result = [taxIDs_doc[x] for x in locus.cluster_prot_ids]
        locus.tax_id = [x['taxid'] for x in tax_result]
        locus.lca = locus.ncbi_taxonomy.LCA(locus.tax_id)
        locus.lca_organism = locus.ncbi_taxonomy.taxid_to_taxonomy(locus.lca)['scientific_name'] if locus.lca else ''


def yates_normalization(samples):
    
    # See if these are n15 quantified
    quant_key = {sample.sample_name: sample['quant'] for sample in samples.values()}
    
    toppep = pd.DataFrame()
    for sample_name in quant_key:
        if quant_key[sample_name]:
            sample = samples[sample_name]

            # do we want n14 or n15 specs?
            if not sample['n15']:
                spec = 'l_spec'
            else:
                spec = 'h_spec'

            # get 500 top peptides in each sample
            top_500 = (sample.pep_quant
                             .T
                             .sort_values(spec, ascending=False)
                             .iloc[:500][spec]
                             .rename(sample_name))
        else:
            top_500 = pd.Series(samples[sample_name].pep_quant).sort_values(ascending=False).iloc[:500]
        toppep = toppep.join(top_500, how='outer')
    

    # Count how many samples each peptide appears in
    toppep['sample_counts'] = toppep.T.count().T

    # Total the peptides and get the mean
    toppep['avg'] = toppep.apply(lambda r: r.drop('sample_counts').mean(), axis=1)
    
    # Only get those in more than half the samples
    common_pep = toppep.query('sample_counts >= {}'.format(len(samples)/2))
    print("{} common peptides for Normalizaion".format(common_pep.shape[0]))
    
    return (common_pep.sum() / common_pep.sum().mean()).to_dict()


def deseq_normalization(grouped_loci, show_plot=False):
    # operates on a list of MultiSampleProteinCluster
    # return norm_factors and calls MultiSampleProteinCluster.normalize on each locus
    # http://pubs.acs.org/doi/full/10.1021/pr401249b
    """ DESeq, a ratio is calculated for each protein by dividing the counts of a 
    protein in a given sample by the geometric mean of counts for that protein across all samples. 
    Finally, each count is corrected by dividing it by the median of all ratios determined in the
     corresponding sample."""
    from scipy.stats import mstats
    for locus in grouped_loci:        
        # Add 1 to counts so no divide by zero errors later on
        counts = [x['counts'] for x in locus.quantification.values() if x['counts'] > 0]
        if len(counts) > 0:
            this_gmean = mstats.gmean(counts)
        else:
            this_gmean=1
        locus.locus_norm = {sample: value['counts']/this_gmean for sample,value in locus.quantification.items()}
    ratios = dict()
    for locus in grouped_loci:
        for key,value in locus.locus_norm.items():
            ratios.setdefault(key,[]).append(value)
    plt.boxplot(list(ratios.values()), labels=list(ratios.keys()))
    norm_factors = {sample:np.median(value) for sample, value in ratios.items()}
    for locus in grouped_loci:
        locus.normalize(norm_factors)
    return norm_factors

make_absolute = lambda path,root_dir: path if os.path.isabs(path) else os.path.join(root_dir, path) if path else ''

def read_metadata(metadata_path):
    metadata = pd.read_csv(metadata_path, index_col=0).T
    metadata = metadata.replace("TRUE",True)
    metadata = metadata.replace("FALSE",False)
    metadata.loc['path'].fillna('',inplace=True)

    # Make paths in 'path' field absolute
    metadata_path = os.path.abspath(metadata_path)
    root_dir = os.path.dirname(metadata_path)
    metadata.loc['path'] = metadata.loc['path'].apply(make_absolute, root_dir=root_dir)
    
    return metadata

def format_qdict(qdict):
    reform = defaultdict(dict)
    for sampName, sampDict in qdict.items():
        for seqName, seqDict in sampDict.items():
            for qType, qVal in seqDict.items():
                reform['\n'.join([sampName, qType])].update({seqName:qVal})
    return pd.DataFrame(reform).fillna(0)

# maybe to_json and to_df should go into a `DataSet` class? Which is a list of MultiSampleProteinClusters?
def to_json(protein_clusters, samples, json_filename, functionizer=None, norm=True):
    """
    List of MultiSampleProteinCluster to json format needed by datatables dta-select project
    """
    import json
    import numpy
    class SetEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, set):
                return list(obj)
            if isinstance(obj, numpy.integer):
                return int(obj)
            return json.JSONEncoder.default(self, obj)

    if functionizer:
        all_pfam = set(chain(*[x.annotations.get("Pfam",[]) for x in protein_clusters]))
        Pfam_info = functionizer.pfam_info(all_pfam)
    else:
        Pfam_info = None
        
    data = [x.as_dict() for x in protein_clusters]
    
    for cluster in data:
        pt = format_qdict(cluster['cluster_peptides'])
        cluster['peptide_table'] = pt.to_html()
        cluster['peptides'] = ';' + ';'.join(pt.index) + ';'
        cluster['max_quant'] = round(max([x['counts'] for x in cluster['quantification'].values()]))

        # Fix Ratio values so no np.nan appears
        if np.isnan(cluster['avg_ratio']):
            cluster['avg_ratio'] = 0
        for samp, values in cluster['quantification'].items():
            if np.isnan(values['ratio']):
                cluster['quantification'][samp]['ratio'] = 0
        
        #pt['peptide'] = pt.index        
        records = pt.to_dict('split')
        records['aoColumns'] = [{"sTitle":x} for x in records['columns']]
        records['aaData'] = records['data']
        
        # add the peptide column on at the front
        records['aoColumns'] = [{"sTitle":"peptide"}] + records['aoColumns']
        records['aaData'] = [[records['index'][idx]]+data for idx,data in enumerate(records['data'])]
        cluster['peptide_records'] = records
        
        
    if norm:
        for pc in data:
            for key in pc['quantification'].keys():
                pc['quantification'][key]['counts'] = int(round(pc['norm_quantification'][key]))

    with open(json_filename,'w') as f:
        json.dump({'data':data, "samples":samples, "Pfam_info": Pfam_info}, f, indent=1, cls=SetEncoder)

def to_df(protein_clusters, norm=True):
    """
    List of MultiSampleProteinCluster to pandas dataframe
    If `norm`, use `norm_quantification` field, else use `quantification` field
    """
    if norm:
        return pd.DataFrame({x.cluster_id:x.norm_quantification for x in protein_clusters}).fillna(0)
    else:
        return pd.DataFrame({x.cluster_id:x.quantification for x in protein_clusters}).fillna(0)

def is_good_db(s):
    # ['RefSeq','UniProt*', 'HMP_Reference_Genomes']
    return True if "refseq" in s.lower() or "uniprot" in s.lower() or s.lower()=="hmp_reference_genomes" else False

def get_good_name(protIDs, protDB):
    # from a list of protIDs, return the description for the largest protein from a good db or from anydb 
    p_result = list(protDB.find({'_id':{'$in':list(protIDs)}}))
    good_result = [x for x in p_result if is_good_db(x['r'])]
    if good_result:
        p_result = good_result
    return max([(len(p['s']),p['d']) for p in p_result], key=lambda x:x[0])[1]
    
def build_loci_from_all_peptides(all_peptides, ppp=2, seqDB=None, group_subsets=True, verbose=False):
    """
    From a list of peptides (probably parsed from a dtaselect-filter file), build loci.


    :param all_peptides:
    :param ppp:
    :param seqDB: 
        pymongo.Collection
    :param group_subsets:
    :param verbose:
    :return:

    If group_subsets is False, a defaultdict is returned where
        key: frozenset of peptide sequences
        value: list of prot_ids matching those peptides
    Example: {'DALDDAFFEEGK', 'IPYVSSPR', 'STGIGDTVFVGPEPEFFIFDSVK'}: [23667280, 13978755, 4937114, 14451791, 79216136, 14754104, 18933409]

    If group_subsets if True: Returns a list of locus dicts.
    A `locus` contains the following fields:
        forward_loci: list of prot_ids in this locus NOT including subset loci
        root_protID: shadow of 'forward_loci'
        protID: list of prot_ids in this locus *including* subset loci
        id: string. all peptides sorted and concated together
        peptide_seq: set of peptide sequences
        subset: list of loci that are a subset of this locus
    """

    def group_subset(pep_sets, verbose=False):
        pep_sets = sorted(pep_sets, key=lambda x: len(x))
        new_pep_sets = []
        while len(pep_sets) > 0:
            if verbose:
                if len(pep_sets) % 100 == 0:
                    print(str(len(pep_sets)))
            # size len 2 can't be subset of size len 2
            if len(pep_sets[-1]) == 2:
                new_pep_sets.extend([(p,) for p in pep_sets])
                break
            root_pep = pep_sets.pop()
            subs = []
            for x in pep_sets:
                if x.issubset(root_pep):
                    subs.append(x)
                    pep_sets.remove(x)
            new_pep = (root_pep, subs)
            new_pep_sets.append(new_pep)
        return new_pep_sets

    all_peptides = sorted(list(all_peptides))

    if seqDB == None:
        seqDB = MongoClient("wl-cmadmin", 27018)["SeqDB_072114"]["SeqDB_072114"]
    if verbose: print('Querying seqDB')
    # Get entries in seqDB for each peptide
    # Example seqDB entry: {"_id" : "AAAAAAAAAA", "p" : [ {"r" : "PTG", "d" : true, "i" : 690796, "o" : 475, "l" : "AAK"}] }
    seq_dict = {x['_id']: [int(i['i']) for i in x['p'] if 'd' not in i] for x in
                seqDB.find({'_id': {'$in': all_peptides}})}
    if verbose: print('Counting peptides')
    count = defaultdict(set)
    for peptide, ids in seq_dict.items():
        for id in ids:
            count[id].add(peptide)
    # Frozenset is hashable and can be used as dict keys (below)
    # count = {protID: frozenset of peptide sequences for that ID}
    count = dict((key, frozenset(value)) for key, value in count.items() if
                 len(value) >= ppp)  # keep only those pIDs with 2 or more peptides
    # Organize into loci
    # key: set of peptides, value: list of protDB_IDs that contain these peptides
    if verbose: print('Assembling into loci')
    loci = defaultdict(list)
    for (id, p) in count.items():
        loci[p].append(id)

    """
    # If you pass in sample_peptides, which is a dict where keys are samples, values: set of peptides in that sample, then filter out any
    # loci with less than "ppp" peptides in any one sample
    if sample_peptides:
        if verbose: print('Before subsets, before filtering: ' + str(len(loci)))
        loci = {locus: value for (locus, value) in loci.items() if
                any([sum([p in peptides for p in locus]) >= 2 for peptides in sample_peptides.values()])}
        if verbose: print('Before subsets, after filtering: ' + str(len(loci)))
    else:
        if verbose: print('Warning: not filtering loci with less than ' + str(ppp) + ' peptides in any one sample')
    """

    if not group_subsets:
        return loci

    if verbose: print('Grouping subsets')
    pep_sets = [set(x) for x in loci.keys()]
    pep_sets_grouped = group_subset(pep_sets, verbose=verbose)

    if verbose: print('Cleanup')
    ps = []
    for pep_set in pep_sets_grouped:
        if len(pep_set) == 2:
            x = {'protID': loci[frozenset(pep_set[0])], 'peptide_seq': pep_set[0],
                 'subset': [{'protID': loci[frozenset(x)], 'peptide_seq': x} for x in pep_set[1]]}
        else:
            x = {'protID': loci[frozenset(pep_set[0])], 'peptide_seq': pep_set[0], 'subset': []}
        ps.append(x)

    for p in ps:
        p['root_protID'] = list(p['protID'])
        p['protID'] = list(chain(*[x['protID'] for x in p['subset']])) + list(p['root_protID'])
        p['id'] = ','.join(str(all_peptides.index(peptide)) for peptide in sorted(list(p['peptide_seq'])))

        # for compatability with blazmass_tools.dta_select_parser
        p['forward_loci'] = p['root_protID']

    return ps

