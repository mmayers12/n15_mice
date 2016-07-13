from itertools import chain
from .common import get_good_name

class ProteinCluster():
    """
    # ProteinCluster for one sample:

    cluster_id: Int
        _id for protein cluster in clusterDB
    sample: build_loci.Sample
        Sample class
    cluster_prot_ids: List[Int]
        protIDs in this cluster that are supported by peptides in this sample (not necesarily all possible protIDs in `cluster_id`)
    peptide_seq: Set[Str]
        set of peptides that support `cluster_prot_ids` in this sample
    peptide_quant: Dict[Str:Int] -> {'PEPTIDE': spectral_count}
        spectral counts of peptides in `peptide_seq`
    quantification: Int
        sum of spectral counts of peptides in `peptide_seq`

    In each locus, the quantification is the sum of the
    quantification for the peptides in each sample. So a peptide can "participate"
    in multiple loci
    """

    
    def __init__(self, cluster_doc, sample):
        # we don't need everyhing in sample. Just the metadata
        self.sample = dict(sample.items())
        self.sample_name = sample.sample_name
        self.db_info = sample.db_info
        self.cluster_id = cluster_doc['_id']
        self.all_cluster_prot_ids = set(cluster_doc['pID'])
        self.cluster_prot_ids = self.all_cluster_prot_ids & sample.prot_ids
        self.peptide_seq = set(chain(*[sample.prot_to_pep[prot_id] for prot_id in self.cluster_prot_ids]))
        
        if self.sample['quant']:
            if self.sample['n15']:
                self.peptide_quant = {seq: self.ratio_pep_quant(sample.pep_quant[seq], self.sample['n15']) for seq in self.peptide_seq}
                self.quantification = self.ratio_prot_quant()
        else:    
            self.peptide_quant = {seq: sample.pep_quant[seq] for seq in self.peptide_seq}
            self.quantification = int(sum(list(self.peptide_quant.values())))
            
    def __str__(self):
        return "cluster id: {}".format(self.cluster_id)
    
    def __repr__(self):
        return self.__dict__.__repr__()
    
    def lookup_name(self):
        # pick the largest protein in one of dbs: ['RefSeq','UniProt*', 'HMP_Reference_Genomes']
        return get_good_name(self.cluster_prot_ids, self.sample.db_info.protDB)

    
    def ratio_pep_quant(self, pep, n15):
        
        import numpy as np        

        new_cols = ['ratio', 'reg_fact', 'counts']
        if n15:
            cols = ['rev_slope_ratio', 'regression_factor', 'h_spec']
        else:
            cols = ['ratio', 'regression_factor', 'l_spec']
        
        quant = pep[cols].rename({o:n for o, n in zip(cols, new_cols)}).to_dict()
        
        if np.isnan(quant['ratio']):
            quant['ratio'] = (10+pep['l_spec']) / (10+pep['h_spec'])
            quant['reg_fact'] = 0            
            if n15:
                quant['ratio'] = 1/quant['ratio']
            
        return quant
    
    def ratio_prot_quant(self):

        num = 0
        denom = 0
        counts = 0        
        for seq in self.peptide_quant.values():
            r = seq['ratio']
            c = seq['counts']
            f = seq['reg_fact']
            num += r*f
            denom += f
            counts += c
        
        ratio = num/denom
        
        return {'ratio' : ratio, 'counts': counts}