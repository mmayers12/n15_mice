import numpy as np
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
            self.peptide_quant = {seq: self.ratio_pep_quant(sample.pep_quant[seq], self.sample['n15']) for seq in self.peptide_seq}
            self.quantification = self.ratio_prot_quant()

        else:
            self.peptide_quant = {seq: {'counts': sample.pep_quant[seq], 'ratio': np.nan, 'reg_fact': np.nan}
                                        for seq in self.peptide_seq}
            self.quantification = {'counts': int(sum(list(x['counts'] for x in self.peptide_quant.values()))), 'ratio': np.nan}

    def __str__(self):
        return "cluster id: {}".format(self.cluster_id)

    def __repr__(self):
        return self.__dict__.__repr__()

    def __eq__(self, other):
        return self.cluster_id == other.cluster_id

    def __ne__(self, other):
        return not self.__eq__(other)

    def lookup_name(self):
        # pick the largest protein in one of dbs: ['RefSeq','UniProt*', 'HMP_Reference_Genomes']
        return get_good_name(self.cluster_prot_ids, self.sample.db_info.protDB)


    def ratio_pep_quant(self, pep, n15):

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

        return quant

    def ratio_prot_quant(self):

        # extract peptide quant values
        ratios = [q['ratio'] for q in self.peptide_quant.values()]
        counts = [q['counts'] for q in self.peptide_quant.values()]
        reg_facts = [q['reg_fact'] for q in self.peptide_quant.values()]

        # Get aggergate quant values
        num = sum([r*f for r, f in zip(ratios, reg_facts)])
        denom = sum(reg_facts)
        tot_counts = sum(counts)

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


        return {'ratio' : ratio, 'counts': tot_counts}
