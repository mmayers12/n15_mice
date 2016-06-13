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
        self.peptide_quant = {seq: sample.pep_quant[seq] for seq in self.peptide_seq}
        self.quantification = int(sum(list(self.peptide_quant.values())))
        
    def __str__(self):
        return "cluster id: {}".format(self.cluster_id)
    
    def __repr__(self):
        return self.__dict__.__repr__()
    
    def lookup_name(self):
        # pick the largest protein in one of dbs: ['RefSeq','UniProt*', 'HMP_Reference_Genomes']
        return get_good_name(self.cluster_prot_ids, self.sample.db_info.protDB)

