from .common import is_good_db
from ..functional_analysis import Functionizer
from .. import taxonomy
from itertools import chain
from collections import defaultdict
import numpy as np



def group_across_samples(protein_clusters, sample_pep_quant, sample_quant_dict=None, db_info=None):
    """Build loci from clustering groups across all samples
    accepts lists of build_loci.ProteinCluster"""

    # sample_pep_quant: dict of {'sample_name': {'pep': quant}}
    # can be built like: sample_pep_quant = {sample.sample_name:sample.pep_quant for sample in samples.values()}
    # build sample_quant_dict like this: sample_quant_dict = {sample: data['n15'] for sample, data in samples.items()}
    # where samples is a dict of {sample_name: build_loci.Sample}
    
    sample_names = set(sample_pep_quant.keys())
    protein_clusters = list(protein_clusters)

    # test if all of the protein_clusters have the same db_info
    if db_info == None and len(set([tuple(x.db_info.items()) for x in protein_clusters])) == 1:
        # if so, use it. Otherwise, must pass it explicitly
        db_info = protein_clusters[0].db_info
    elif db_info == None:
        raise ValueError("Must pass db_info")

    all_cluster_ids = set([x.cluster_id for x in protein_clusters])  # set of all clusters in this dataset

    if not MultiSampleProteinCluster.inited:
        functionizer = Functionizer(protDB=db_info.protDB, domainDB=db_info.domainDB, hashDB=db_info.hashDB)
        ncbi_taxonomy = taxonomy.Taxonomy(mongo_coll = db_info.ncbi_taxonomy)
        MultiSampleProteinCluster.init(ncbi_taxonomy = ncbi_taxonomy, functionizer = functionizer)

    id_clusters = defaultdict(list)
    for x in protein_clusters:
        id_clusters[x.cluster_id].append(x)

    grouped_loci = []
    for cluster_id in all_cluster_ids:
        # the ProteinClusters across all samples for this cluster_id
        this_clusters = id_clusters[cluster_id]


        #cluster_peptides = {x.sample_name: x.peptide_quant for x in this_clusters}
        #quantification = {x.sample_name:x.quantification for x in this_clusters}


        
        peptide_seq = set(chain(*[x.peptide_seq for x in this_clusters]))
        cluster_peptides = dict()
        for sample_name in sample_names:
            pep_quant = {pep:ratio_pep_quant(sample_pep_quant[sample_name][pep], sample_quant_dict[sample_name]) 
                                             for pep in peptide_seq if pep in sample_pep_quant[sample_name]}
            if pep_quant:
                cluster_peptides[sample_name] = pep_quant
        
        quantification = {sample_name:ratio_prot_quant(pep_quant) for sample_name,pep_quant in cluster_peptides.items()}
        

        cluster_prot_ids = set(chain(*[x.cluster_prot_ids for x in this_clusters]))
        grouped_loci.append(MultiSampleProteinCluster(quantification=quantification, cluster_id=cluster_id,
                             cluster_peptides=cluster_peptides, cluster_prot_ids=cluster_prot_ids, db_info=db_info))
    return grouped_loci


def ratio_pep_quant(pep, n15):
    # Columns with important quant info
    if n15:
        cols = ['rev_slope_ratio', 'regression_factor', 'h_spec']
    else:
        cols = ['ratio', 'regression_factor', 'l_spec']

    # Regularize column names across n14 and n15 samples
    new_cols = ['ratio', 'reg_fact', 'counts']
    quant = pep[cols].rename({o:n for o, n in zip(cols, new_cols)}).to_dict()

    # If there is no ratio do a ratio of specs... but add 5 to both
    # to reduce importance low count ratios e.g 2:1 or 1:0 or 2:0
    if np.isnan(quant['ratio']):
        quant['ratio'] = (5+pep['l_spec']) / (5+pep['h_spec'])
        quant['reg_fact'] = 0
        if n15:
            quant['ratio'] = 1/quant['ratio']
    
    # Back Calc counts from ratio
    # n15 back_counts = n14counts * n15_ratio AND n14_backCounts = n14Counts
    # if no n14 counts, n15 back_counts = n15 coutns, and n14_back_counts = n15Counts * n14_ratio
    if n15:
        quant['back_calc'] = pep['l_spec'] * quant['ratio']
        if quant['back_calc'] == 0:
            quant['back_calc'] = pep['h_spec']
    else:
        quant['back_calc'] = pep['l_spec']
        if quant['back_calc'] == 0:
            quant['back_calc'] = pep['h_spec'] * quant['ratio']
        

    return quant


def ratio_prot_quant(peptide_quant):

    # extract peptide quant values
    ratios = [q['ratio'] for q in peptide_quant.values()]
    counts = [q['counts'] for q in peptide_quant.values()]
    reg_facts = [q['reg_fact'] for q in peptide_quant.values()]
    back_calcs = [q['back_calc'] for q in peptide_quant.values()]

    # Get aggergate quant values
    num = sum([r*f for r, f in zip(ratios, reg_facts)])
    denom = sum(reg_facts)
    tot_counts = sum(counts)
    tot_back_calc = sum(back_calcs)

    # If no census ratios, but spec counts are high, might be able to
    # infer a ratio from these
    if denom == 0:
        if tot_counts > 5:
            # Solves for l+h spec eqn: ratio = this_spec+5 / other_spec+5
            weights = [ (((c+5)/r)-5)+c for c,r in zip(counts, ratios)]
            num = sum([w*r for w,r in zip(weights,ratios)])
            denom = sum(weights)
        else:
            num, denom = (np.nan, np.nan)

    # get weighted ratio for protein
    ratio = num/denom


    return {'ratio' : ratio, 'counts': tot_counts, 'back_calc': tot_back_calc}


class MultiSampleProteinCluster():
    """
    A ProteinCluster is built from group_across_samples based on clustering specified in clusterdb
    cluster_id: Int
    cluster_peptides: Dict[Str:Dict[Str:Int]] -> {sample: {peptide: spectral_count}}
    cluster_prot_ids: protIDs that are supported by peptides in any sample (not necesarily all possible protIDs in this cluster)
    quantification:  {sample: value}, where value is the sum of all peptide spectral counts in that sample for peptides in cluster_peptides

    Example:
    from metaproteomics.analysis.build_loci import ProteinCluster
    ProteinCluster(cluster_id = 1, cluster_prot_ids = [1,2,3],
                   cluster_peptides = {'sampleA':{'aaa': {'counts':100, 'ratio':1.213, 'reg_fact': .8892},
                                                  'bbb': {'counts':23, 'ratio':0.866, 'reg_fact': .9421},
                                       'sampleB':{'aaa': {'counts':101, 'ratio':1.271, 'reg_fact': .8747},
                                                  'ccc': {'counts':123, 'ratio':3.741, 'reg_fact': .9548}},
                   quantification = {'sampleA': {'counts': 123, 'ratio': 1.304},
                                     'sampleB': {'counts': 224, 'ratio': 2.560})

    """
    functionizer = None
    ncbi_taxonomy = None
    go_ontology = None
    inited = False

    # attempting to make whats sets in init pickleable.
    # not quite working.....
    """
    def __getstate__(self):
        state = dict(self.__dict__)
        state['functionizer'] = self.functionizer
        state['go_ontology_name'] = self.obo_parser.__name__
        state['ncbi_tax_host'] =
        state['ncbi_tax_port'] =
        return state

    def __setstate__(self, d):
        self.__dict__.update(d)
        obo_parser = __import__(state['go_ontology_name'])

        MultiSampleProteinCluster.init(ncbi_taxonomy = taxonomy.Taxonomy(
                                                host=state['ncbi_tax_host'], port=state['ncbi_tax_port']),
                                       functionizer = functionizer, go_ontology = obo_parser.GODag())
    """
    @classmethod
    def init(cls, ncbi_taxonomy=None, functionizer=None, go_ontology=None):
        from ..functional_analysis import Functionizer
        from .. import taxonomy
        from ... import obo_parser

        cls.ncbi_taxonomy = ncbi_taxonomy if ncbi_taxonomy else taxonomy.Taxonomy(host='wl-cmadmin', port=27017)
        cls.functionizer = functionizer if functionizer else Functionizer()
        cls.go_ontology = go_ontology if go_ontology else obo_parser.GODag()
        cls.inited = True

    def __init__(self, cluster_id, cluster_peptides, cluster_prot_ids, quantification, db_info):
        if not MultiSampleProteinCluster.inited:
            MultiSampleProteinCluster.init()
        self.cluster_id = cluster_id
        self.cluster_peptides = cluster_peptides
        self.cluster_prot_ids = cluster_prot_ids
        self.quantification = quantification
        self.db_info = db_info

    def annotate(self):
        self.lookup_name()
        self.lookup_function()
        self.get_go()
        self.get_tax_info()
        self.quantify()

    def quantify(self, samples=None):
        """
        Get average ratio and p-value from ratio t-test
        Ratio t-test: log-transform ratios, then t-test against 0
        samples: a subset of samples to quantify across e.g if both heavy and
        light samples included, it only makes sense to average ratios across one
        set or the other.
        """
        from scipy import stats

        if samples:
            ratios = [q['ratio'] for s,q in self.quantification.items() if s in samples]
        else:
            ratios = [q['ratio'] for q in self.quantification.values()]

        # Filter out 0 and NaN
        ratios = np.array([r for r in ratios if r != 0])
        ratios = ratios[~np.isnan(ratios)]

        log_ratios = np.log(ratios)

        t, p = stats.ttest_1samp(log_ratios, 0)

        self.avg_ratio = np.mean(ratios)
        self.p_value = p

    # get a "representative" locus
    def lookup_name(self, all_prot_result=None):

        if all_prot_result:
            prot_result = [all_prot_result[x] for x in self.cluster_prot_ids]
        else:
            prot_result = list(self.db_info.protDB.find({'_id':{'$in':list(self.cluster_prot_ids)}}))

        # pick the largest protein in one of dbs: ['RefSeq','UniProt*', 'HMP_Reference_Genomes']
        self.prot_info = [{'_id':p['_id'], 'r':p['r'], 'd':p['d'], 's_len':len(p['s'])} for p in prot_result]
        good_result = [x for x in self.prot_info if is_good_db(x['r'])]
        proteins = good_result if good_result else self.prot_info
        self.name = max([(p['s_len'],p['d']) for p in proteins], key=lambda x:x[0])[1]
        #Keep the full sequence of the winner
        self.seq = [x['s'] for x in prot_result if x['d'] == self.name][0]

    def lookup_function(self):
        assert self.functionizer, "Must initialize functionizer"
        self.annotations = self.functionizer.get_annotations_from_protIDs(self.cluster_prot_ids)

    def get_go(self):
        if "annotations" not in self.__dict__:
            self.lookup_function()
        assert self.go_ontology, "Must initialize go_ontology"
        self.go = [self.go_ontology[x] for x in self.annotations.get('go',[])]
        self.go_names = [str(x) for x in self.go]
        # These go terms aren't pickleble because of recursion.
        # TODO: fix this somehow?
        del self.go

    def get_tax_info(self):
        # set tax_id, lca, lca_organism
        taxIDs_doc = list(self.db_info.taxDB.aggregate(
            [{'$match': {'_id': {'$in': list(self.cluster_prot_ids)}}},
             {'$group': {'_id': None, 'taxid': {'$addToSet': '$taxid'}}}]))
        self.tax_id = [x for x in taxIDs_doc[0]['taxid'] if x] if taxIDs_doc else []
        self.lca = self.ncbi_taxonomy.LCA(self.tax_id)
        self.lca_organism = self.ncbi_taxonomy.taxid_to_taxonomy(self.lca)['scientific_name'] if self.lca else ''


    def __str__(self):
        return """Cluster ID: {cluster_id}\nQuantification: {q}\nGO Terms: {g}""".format(cluster_id=self.cluster_id, q=self.quantification, g=self.go_names)

    def __repr__(self):
        return self.__dict__.__repr__()

    def passes_thresh(self, metadata, min_quant = 2, ratio_ok = True, ratio_only = False, min_samples = 2, min_samples_per_group = 2, group = "biological"):
        """
        min_quant: Must have a minimum quantification of `min_quant` in at least one sample
        ratio_ok: Allow presence of a ratio, rather than `min_quant`, to allow threshold passing
        min_samples: Must have a minimum quantification of `min_quant` in at least `min_samples` samples
        min_samples_per_group: Must be at least `min_samples_per_group` samples that are in the same group `group`,
                                having a minimum quantification of `min_quant`
        """

        if ratio_only:
            num_ratios = (~np.isnan([x['ratio'] for x in self.quantification.values()])).sum()
            if num_ratios == 0:
                return False
            if len([x['counts'] for x in self.quantification.values() if ~np.isnan(x['ratio'])]) < min_samples:
                return False
            if min_samples_per_group:
                samples_that_pass = [sample for sample,value in self.quantification.items() if  ~np.isnan(value['ratio'])]
                if not any(metadata[samples_that_pass].loc[group].value_counts() >= min_samples_per_group):
                    return False
            return True

        if ratio_ok:
            num_ratios = (~np.isnan([x['ratio'] for x in self.quantification.values()])).sum()
            if num_ratios == 0:
                if not any([x['counts']>=min_quant for x in self.quantification.values()]):
                    return False
            if len([x['counts'] for x in self.quantification.values() if (x['counts']>=min_quant or ~np.isnan(x['ratio']))]) < min_samples:
                return False
            if min_samples_per_group:
                samples_that_pass = [sample for sample,value in self.quantification.items() if (value['counts']>=min_quant or ~np.isnan(value['ratio']))]
                if not any(metadata[samples_that_pass].loc[group].value_counts() >= min_samples_per_group):
                    return False
            return True

        if not any([x['counts']>=min_quant for x in self.quantification.values()]):
            return False
        if len([x['counts'] for x in self.quantification.values() if x['counts']>=min_quant]) < min_samples:
            return False
        if min_samples_per_group:
            samples_that_pass = [sample for sample,value in self.quantification.items() if value['counts']>=min_quant]
            if not any(metadata[samples_that_pass].loc[group].value_counts() >= min_samples_per_group):
                return False
        return True

    def as_dict(self):
        d = {}
        keys = ['prot_info','tax_id','lca','lca_organism','go','go_names','annotations','quantification',
                'cluster_id','cluster_peptides','cluster_prot_ids','name','peptide_table','max_quant',
                'crapome','seq','norm_quantification', 'avg_ratio', 'p_value']
        for key in keys:
            if key in self.__dict__:
                d[key] = self.__dict__[key]
        return d

    def normalize(self, norm_factors, field = "norm_quantification"):
        self.__dict__[field] = {sample: value['counts']/norm_factors[sample] for sample,value in self.quantification.items()}
