from itertools import chain, repeat
import pandas as pd
from ... import blazmass_tools
from .common import build_loci_from_all_peptides
from .ProteinCluster import ProteinCluster

class Sample(dict):
    # Stores sample information parsed from dtaselect-filter
    # The peptides are parsed from dtaselect-filter and then re-built into loci using seqDB
    # Also stores sample metadata
    """
    *Constructor Parameters*
    sample_name: Str
        Sample Name
    dtaselect_path: path to DTASelect-filter.txt file
        DTASelect-filter.txt file
    db_info: metaproteomics.analysis.DBInfo
        Class containing mongodb database information
        
    *Optional parameters*
    metadata: Dict or pandas.Series or pandas.DataFrame
        Dictionary containing metadata associated with this sample. pandas.Series will be converted to 
        a dict. If passed a pandas.DataFrame, the column matching `sample_name` will be used.
    ppp: Int
        peptides per protein (default: 2)

        
    *Attributes*
    dta_select: metaproteomics.file_processing.blazmass_tools.dta_select_parser Dict
        Parsed DTASelect-filter.txt file
    peptides: Set[Str]
        Set of peptides parsed from dta_select
    pep_quant: Dict[Str: Int]
        Spectral count of each peptide in the sample. Built from blazmass_tools.build_pep_quant_dict
    loci_pep_prot: Dict[Set[Str]: List[Int]] -> {set of peptides: list of prot_ids}
        loci built dta-select-like. Ex: {{'DALDDAFFEEGK', 'IPYVSSPR'}: [23667280, 13978755, 4937114]}
    prot_ids: Set[Int]
        Set of protIDs in this sample using the info in `loci_pep_prot`
    prot_to_pep: Dict[Int: Set[Str]] -> {prot_id: {'pepA','pepB'}}
        For each protid, we want the peptides that support it
    """
    def __init__(self, sample_name, dtaselect_path, db_info, metadata = {}, ppp = 2):
        self.sample_name = sample_name
        self.db_info = db_info
        self.update(blazmass_tools.dta_select_header(dtaselect_path))
        self.set_metadata(metadata)
        
        """ We don't need this do we????
        dta_len = len(list(blazmass_tools.dta_select_parser(dtaselect_path,small=True, get_tax=False)))
        self.dta_select = list(tqdm(blazmass_tools.dta_select_parser(dtaselect_path, 
                                    get_tax=True, taxDB=db_info.taxDB, protDB=db_info.protDB),total=dta_len, nested=True))"""

        self.dta_select = list(blazmass_tools.dta_select_parser(dtaselect_path,get_tax=False))
        self.pep_quant = blazmass_tools.build_pep_quant_dict(self.dta_select, field='Redundancy')
        self.peptides = set(self.pep_quant.keys())
        
        self.loci_pep_prot = build_loci_from_all_peptides(self.peptides, group_subsets=False, ppp=ppp, seqDB=db_info.seqDB)
        self.prot_ids = set(chain(*self.loci_pep_prot.values()))
        self.prot_to_pep = {protid: set(peptides) for (peptides, protid) in chain(
            *[zip(repeat(peptides), prot_ids) for (peptides, prot_ids) in self.loci_pep_prot.items()])}
    
    @classmethod
    def from_peptides(cls, sample_name, pep_quant, db_info, metadata = {}, ppp = 2):
        self = cls.__new__(cls)
        self.sample_name = sample_name
        self.db_info = db_info
        self.set_metadata(metadata)
        self.dta_select = []
        self.pep_quant = pep_quant
        self.peptides = set(pep_quant.keys())
        self.loci_pep_prot = build_loci_from_all_peptides(self.peptides, group_subsets=False, ppp=ppp, seqDB=db_info.seqDB)
        self.prot_ids = set(chain(*self.loci_pep_prot.values()))
        self.prot_to_pep = {protid: set(peptides) for (peptides, protid) in chain(
            *[zip(repeat(peptides), prot_ids) for (peptides, prot_ids) in self.loci_pep_prot.items()])}
        return self
    
    def set_metadata(self, metadata):
        if isinstance(metadata,pd.core.frame.DataFrame):
            if self.sample_name in metadata:
                metadata = metadata[self.sample_name]
            else:
                raise ValueError("{} not found in metadata DataFrame".format(self.sample_name))
        if isinstance(metadata,pd.core.frame.Series):
            metadata = dict(metadata)
        self.update(metadata)
    
    def summary(self):
        summary = pd.DataFrame(columns=['Proteins','Peptides','Spectra'], index=['Filtered','Forward','Decoy','FDR'])
        decoy_peptides = set(chain(*[x['unmod_peptide_seq'] for x in self.dta_select if x['reverse']]))
        summary.Peptides.Decoy = len(decoy_peptides)
        forward_peptides = set(chain(*[x['unmod_peptide_seq'] for x in self.dta_select if not x['reverse']]))
        summary.Peptides.Forward = len(forward_peptides)
        summary.Peptides.Filtered = summary.Peptides.Decoy + summary.Peptides.Forward
        summary.Peptides.FDR = 100 * summary.Peptides.Decoy / summary.Peptides.Forward
        
        summary.Proteins.Decoy = len([x for x in self.dta_select if x['reverse']])
        summary.Proteins.Forward = len([x for x in self.dta_select if not x['reverse']])
        summary.Proteins.Filtered = summary.Proteins.Decoy + summary.Proteins.Forward
        summary.Proteins.FDR = 100 * summary.Proteins.Decoy / summary.Proteins.Forward
        
        summary.Spectra.Decoy = sum(self.pep_quant[peptide] for peptide in decoy_peptides)
        summary.Spectra.Forward = sum(self.pep_quant[peptide] for peptide in forward_peptides)
        summary.Spectra.Filtered = summary.Spectra.Decoy + summary.Spectra.Forward
        summary.Spectra.FDR = 100 * summary.Spectra.Decoy / summary.Spectra.Forward
        
        return summary
    
    def build_protein_clusters(self):
        """ return a list of ProteinClusters from a Sample """
        clusterDB = self.db_info.clusterDB
        return [ProteinCluster(cluster_doc, self) for cluster_doc in clusterDB.find({'pID': {'$in': list(self.prot_ids)}})]
        
        
