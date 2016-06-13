"""

From home run:
$ ssh -fNL 27028:localhost:27018 gbwl
$ ssh -fNL 27027:localhost:27017 gbwl


To Do:
when you have pfam and go sets, merge redundant sets


IPR descriptions:
from analysis import parse_IPR_tree
IPR_terms = parse_IPR_tree.IPR_terms


pfam descriptions:
from interproscan import pfam_sets
pfam_dict = pfam_sets.parse_pfam_clans()
for locus in loci:
    locus['descr'] = {x:pfam_dict[x]['name'] for x in locus['Pfam']}


go descriptions:
from goatools import obo_parser
ontology = obo_parser.GODag()

"""

from itertools import chain
from collections import Counter
import hashlib
import warnings
from pymongo import MongoClient
from collections import defaultdict

class Functionizer:
    
    # So we can pickle a Functionizer instance
    def __getstate__(self):
        state = dict(self.__dict__)
        del state['protDB']
        del state['domainDB']
        del state['hashDB']
        state['protDB_tuple'] = self.protDB.database.client.address + (self.protDB.database.name,self.protDB.name)
        state['domainDB_tuple'] = self.domainDB.database.client.address + (self.domainDB.database.name,self.domainDB.name)
        state['hashDB_tuple'] = self.hashDB.database.client.address + (self.hashDB.database.name,self.hashDB.name)
        return state
    
    def __setstate__(self, d):
        self.__dict__.update(d)
        self.protDB = Functionizer.mongoclient_builder(d['protDB_tuple'])
        self.domainDB = Functionizer.mongoclient_builder(d['domainDB_tuple'])
        self.hashDB = Functionizer.mongoclient_builder(d['hashDB_tuple'])
        
    def __init__(self, protDB=None, domainDB=None, hashDB=None):
        """ Accepts mongo collection objects 
        protDB is only needed to look up protein sequences to hash them, if hashDB is not provided
        domainDB contains interproscan info indexed by protein sequence md5sum
        hashDB: protDB ids <-> md5sum
        """
        default_client = MongoClient('wl-cmadmin', 27018)
        if not protDB:
            protDB = default_client.ProtDB_072114.ProtDB_072114
        if not domainDB:
            domainDB = default_client.DomainDB_072114.DomainDB_072114
        if not hashDB:
            hashDB = default_client.HashDB_072114.HashDB_072114
        assert protDB.count() > 0
        assert domainDB.count() > 0
        if hashDB.count() == 0:
            warnings.warn("HashDB is empty. Please check")
        self.protDB = protDB
        self.domainDB = domainDB
        self.hashDB = hashDB

    
    @staticmethod
    def mongoclient_builder(db_info_tuple):
        return MongoClient(host=db_info_tuple[0], port=db_info_tuple[1])[db_info_tuple[2]][db_info_tuple[3]]

    def compute_hashes(self, protIDs):
        # If you don't have a hashDB, for getting hashes of proteins by computing md5sum
        return set(hashlib.md5(d['s'].encode('utf-8')).hexdigest() for d in self.protDB.find({'_id': {'$in': list(protIDs)}}))
    
    def get_hashes(self, protIDs):
        return set(x['_id'] for x in self.hashDB.find({'pID': {'$in': list(protIDs)}}, {'_id': True}))

    def get_all_go_from_hashes(self, hashes):
        # return set of all go_terms for all domains in list of hashes
        domain_list = self.get_domains_from_hashes(list(hashes))
        return set(chain(*[set(chain(*[x['g'] for x in chain(*[x for x in values.values()]) if 'g' in x])) for values in domain_list]))
    
    def get_domains_from_hashes(self, hashes):
        return [x['d'] for x in self.domainDB.find({'_id': {'$in': list(hashes)}})]
    """ test
    from metaproteomics.analysis import functional_analysis
    f = functional_analysis.Functionizer()
    f.get_annotations_from_protIDs(list(range(10)))

    """
    
    def get_annotations_from_protIDs(self, protIDs, from_db=[], return_ipa=True, return_go=True):
        """
        If `from_db` is empty, use all member databases.
        `from_db` is a list of strings. for example `['pfam','superfamily']`
        
        If `from_db` is given, the GO and/or IPA accessions returned will be only from the databases in `from_db`
        
        Annotations are returned as a dictionary
        """
        hashes = self.get_hashes(protIDs)
        domain_result = self.get_domains_from_hashes(hashes)
    
        if from_db:
            from_db = [x.lower() for x in from_db]
            # Get only the domains from dbs in from_db
            domain_result = [{key: value for key, value in x.items() if key.lower() in from_db} for x in domain_result]
        else:
            # Get all annotations from any dbs
            pass
    
        # Collapse db : domains dict into list of domains
        domain_list = list(chain(*[list(chain(*x.values())) for x in domain_result]))
        if not domain_list:
            return {}
        
        annotations = defaultdict(set)
        if return_ipa:
            _ipa = set([domain['ipa'] for domain in domain_list if 'ipa' in domain])
            if _ipa: annotations['ipa'] = _ipa
        if return_go:
            _go = set(chain(*[domain['g'] for domain in domain_list if 'g' in domain]))
            if _go: annotations['go'] = _go
        for domain in domain_list:
            if 'sa' in domain:
                annotations[domain['a']].add(domain['sa'])
        return dict(annotations)
        
    def pfam_info(self, pfams):
        from ..interproscan import pfam_sets
        if "pfam_dict" not in self.__dict__:
            self.pfam_dict = pfam_sets.parse_pfam_clans()
        return {x:self.pfam_dict[x]['name'] for x in pfams}
        
    
def get_go_set_parent(p):
    # Get the most common set of GO terms for the parent loci only
    parent_hashes = [x['_id'] for x in hashDB.find({'pID': {'$in': p['parent_forward_loci']}}, {'_id': True})]
    parent_domain_result = domainDB.find({'_id': {'$in': parent_hashes}})
    parent_domain_dict = dict([(x['_id'], x['d']) for x in parent_domain_result])

    hash_go_dict = {hash: set(chain(*[x['g'] for x in chain(*[x for x in values.values()]) if 'g' in x])) for
                    hash, values in parent_domain_dict.items()}
    hash_go_dict = {key: value for key, value in hash_go_dict.items() if value != set()}
    # c = Counter(chain(*[list(x) for x in hash_go_dict.values()])) # Counter for all Go terms (Singly)
    hash_go_dict_counter = Counter(
        [frozenset(x) for x in hash_go_dict.values()])  # What are the unique sets of go terms?

    if len(hash_go_dict_counter.most_common(1)) > 1:
        print('Theres a tie')
    return list(hash_go_dict_counter.most_common(1)[0][0]) if hash_go_dict_counter else None


def get_go_set(p):
    # Get the most common set of GO terms for the parent loci
    # If the parent loci has no GO terms, use all the subsets
    parent_go_terms = get_go_set_parent(p)
    if parent_go_terms:
        return parent_go_terms
    hashes = [x['_id'] for x in hashDB.find({'pID': {'$in': p['forward_loci']}}, {'_id': True})]
    domain_dict = dict([(x['_id'], x['d']) for x in domainDB.find({'_id': {'$in': hashes}})])
    hash_go_dict = {hash: set(chain(*[x['g'] for x in chain(*[x for x in values.values()]) if 'g' in x])) for
                    hash, values in domain_dict.items()}
    hash_go_dict = {key: value for key, value in hash_go_dict.items() if value != set()}
    hash_go_dict_counter = Counter([frozenset(x) for x in hash_go_dict.values()])
    if len(hash_go_dict_counter.most_common(1)) > 1:
        print('Theres a tie')
    return list(hash_go_dict_counter.most_common(1)[0][0]) if hash_go_dict_counter else None


def get_all_go_parent(p):
    # Using the parent loci only, return the set of all go_terms for all domains in all hashes
    parent_hashes = [x['_id'] for x in hashDB.find({'pID': {'$in': p['parent_forward_loci']}}, {'_id': True})]
    return get_all_go_from_hashes(parent_hashes)

def get_all_go_parent_sub(p):
    # return the set of all go_terms for all domains in all hashes in all loci, parents and sub
    hashes = [x['_id'] for x in hashDB.find({'pID': {'$in': p['forward_loci']}}, {'_id': True})]
    return get_all_go_from_hashes(hashes)




def get_all_go(p, verbose=False):
    # Using the parent loci only, get the set of all go_terms for all domains in all hashes
    # If the parent loci has no GO terms, use the next largest subsets (by psm count), and so forth
    parent_go_terms = get_all_go_parent(p)
    if parent_go_terms:
        if verbose: print('Using parent only')
        return parent_go_terms
    # Nothing from just using parent. What is the next largest subsets (or more, if tie)
    if len(p['subset']) == 0:
        if verbose: print('No parent match, no subsets to try')
        return None
    subset_lens = [len(x['peptide_seq']) for x in p['subset']]
    next_biggest = max(subset_lens)
    while next_biggest > 0:
        if verbose: print('Trying subsets with ' + str(next_biggest) + ' peptides')
        loci = list(chain(*[x['forward_loci'] for x in p['subset'] if len(x['peptide_seq']) == next_biggest]))
        if len(loci) == 0:
            next_biggest -= 1
            continue
        hashes = [x['_id'] for x in hashDB.find({'pID': {'$in': loci}}, {'_id': True})]
        result = get_all_go_from_hashes(hashes)
        if not result:
            next_biggest -= 1
        else:
            if verbose: print('Found using subsets with ' + str(next_biggest) + ' peptides')
            return result
    if verbose: print('No parent match, no subsets matches')
    return None


def get_annotations(p, from_db=None, return_ipa=True, return_go=False, parent_only=False,
                    root_protID_field='root_protID', protID_field='protID',
                    verbose=False):
    '''
    # Get the set of all annotations for all domains in all hashes from the dbs in from_db
    # If from_db == None, use all member databases.
    # from_db is a list of strings. for example ['pfam','superfamily']
    # If return_ipa: Return set of interproscan acessions,supercedes anything set in from_db
    # else, return member db annotations (collectively)
    # If parent_only == True, query ONLY using parent forward loci.

    # First query using root_loci. If no results, of the type requested. Successively
    # use the next largest subsets (by psm count), and so forth until results or return None.
    # Looks in 'p' for field root_protID_field as the root. Then under 'subset' should be a list of dicts
    # where each contains fields 'peptide_seq' and protID_field

    example: p = {'root_protID': [1000000000], 'subset': [
            {'protID': [1236666666, 1212666666666666], 'peptide_seq': {'LASYLDKVR', 'VTMQNLNDR', 'AAQQW'}},
            {'protID': (61518136, 61518138, 61518135), 'peptide_seq': {'DAEAWFNEK', 'VTMQNLNDR'}},
            {'protID': (10000000), 'peptide_seq': {'DAEAWFNEK', 'VTMQNLNDE'}}]}
    hashes = ['1366487868c1d8201018fdef80e96429', '12ecfb292578bd34964872c51efcec7a']
    '''
    protID = p[root_protID_field]
    annotations = get_annotations_from_protIDs(protID, from_db=from_db, return_ipa=return_ipa, return_go=return_go,
                                               verbose=verbose)

    if annotations:
        if verbose: print('Using parent only')
        return annotations

    if parent_only:
        return annotations

    # Nothing from just using parent. What is the next largest subsets (or more, if tie)
    if len(p['subset']) == 0:
        if verbose: print('No parent match, no subsets to try')
        return None
    subset_lens = [len(x['peptide_seq']) for x in p['subset']]
    next_biggest = max(subset_lens)
    while next_biggest > 0:
        if verbose: print('Trying subsets with ' + str(next_biggest) + ' peptides')
        protID = list(chain(*[x[protID_field] for x in p['subset'] if len(x['peptide_seq']) == next_biggest]))
        if len(protID) == 0:
            next_biggest -= 1
            continue
        annotations = get_annotations_from_protIDs(protID, from_db=from_db, return_ipa=return_ipa, verbose=verbose)
        if not annotations:
            next_biggest -= 1
        else:
            if verbose: print('Found using subsets with ' + str(next_biggest) + ' peptides')
            return annotations
    if verbose: print('No parent match, no subsets matches')
    return None

