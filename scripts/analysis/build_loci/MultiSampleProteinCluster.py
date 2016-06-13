from .common import is_good_db
from ..functional_analysis import Functionizer
from .. import taxonomy
from itertools import chain
from collections import defaultdict


def group_across_samples(protein_clusters, db_info=None):
    """Build loci from clustering groups across all samples
    accepts lists of build_loci.ProteinCluster"""
    
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
        cluster_peptides = {x.sample_name: x.peptide_quant for x in this_clusters}
        quantification = {x.sample_name:x.quantification for x in this_clusters}
        cluster_prot_ids = set(chain(*[x.cluster_prot_ids for x in this_clusters]))
        grouped_loci.append(MultiSampleProteinCluster(quantification=quantification, cluster_id=cluster_id,
                             cluster_peptides=cluster_peptides, cluster_prot_ids=cluster_prot_ids, db_info=db_info))
    return grouped_loci
    
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
                   cluster_peptides = {'sampleA':{'aaa':100,'bbb':23}, 'sampleB':{'aaa':101,'ccc':123}}, 
                   quantification = {'sampleA': 123, 'sampleB': 224})

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
        
    def passes_thresh(self, metadata, min_quant = 2, min_samples = 2, min_samples_per_group = 2, group = "biological"):
        """
        min_quant: Must have a minimum quantification of `min_quant` in at least one sample
        min_samples: Must have a minimum quantification of `min_quant` in at least `min_samples` samples
        min_samples_per_group: Must be at least `min_samples_per_group` samples that are in the same group `group`, 
                                having a minimum quantification of `min_quant`
        """
        if not any([x>=min_quant for x in self.quantification.values()]):
            return False
        if len([x for x in self.quantification.values() if x>=min_quant]) < min_samples:
            return False
        if min_samples_per_group:
            samples_that_pass = [sample for sample,value in self.quantification.items() if value>=min_quant]
            if not any(metadata[samples_that_pass].loc[group].value_counts() >= min_samples_per_group):
                return False
        return True
    
    def as_dict(self):
        d = {}
        keys = ['prot_info','tax_id','lca','lca_organism','go','go_names','annotations','quantification',
                'cluster_id','cluster_peptides','cluster_prot_ids','name','peptide_table','max_quant',
                'crapome','seq','norm_quantification']
        for key in keys:                
            if key in self.__dict__:
                d[key] = self.__dict__[key]
        return d
    
    def normalize(self, norm_factors, field = "norm_quantification"):
        self.__dict__[field] = {sample: value/norm_factors[sample] for sample,value in self.quantification.items()}
