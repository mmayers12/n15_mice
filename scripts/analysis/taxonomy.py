# -*- coding: utf-8 -*-
"""
Created on Tue Oct  7 09:18:00 2014

@author: Greg Stupp

# Example: Parse refseq defline

# from Sandip's db (with db and go info at the end)
taxonomy.fasta_to_organism('>gi|526245011|ref|YP_008320337.1| terminase small [Euphorbia mosaic virus - A [Mexico:Yucatan:2004]]|refseq|[]')

# Or, as a dictionary with source database specified
taxonomy.fasta_to_organism({'d':'>gi|526245011|ref|YP_008320337.1| terminase small [Euphorbia mosaic virus - A [Mexico:Yucatan:2004]]', 'r': 'refseq'})

# or, call the right function
taxonomy.fasta_to_organism_refseq('>gi|526245011|ref|YP_008320337.1| terminase small [Euphorbia mosaic virus - A [Mexico:Yucatan:2004]]')

# Get taxID:
t = taxonomy.Taxonomy(preset = 'wl') # connects to wolan lab's mongoDB
t.organism_to_taxonomy('Euphorbia mosaic virus - A [Mexico:Yucatan:2004]')

"""

import re
import sys
import logging

from pymongo import MongoClient

logger = logging.getLogger('my_logger')
if not logger.handlers:
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARN) # show warning and higher priority in console. Don't show INFO
    ch.setFormatter(logging.Formatter('%(levelname)s - %(funcName)s - %(message)s'))
    logger.addHandler(ch)

class Taxonomy(object):
    def __init__(self, host='localhost', port=27017, mongo_coll = None):
        # setting mongo_coll overrides `host`,`port`, and the default db & coll "taxonomy"
        if mongo_coll:
            self.taxonomy_coll = mongo_coll
        else:
            self.host = host
            self.port = port
            client = MongoClient(self.host, port=self.port)
            self.taxonomy_coll = client.taxonomy.taxonomy

        self.all_scientific_names = []
        if self.taxonomy_coll.find_one() is None:
            print("taxonomy database is empty. You should check this.")
        self.ranks = ['superkingdom', 'kingdom', 'subkingdom', 'superphylum', 'phylum', 'subphylum',
                        'superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order', 'suborder',
                        'infraorder', 'parvorder', 'superfamily', 'family', 'subfamily', 'tribe', 'subtribe', 'genus',
                        'subgenus', 'species group', 'species subgroup', 'species', 'subspecies', 'varietas', 'forma', 'no rank']
        self.ranks = [x for x in self.ranks[::-1]]


    def get_organism_parent_taxid(self, taxid):
        if type(taxid) != int:
            return None
        tax_dict = self.taxonomy_coll.find_one({'taxid': taxid})
        if tax_dict:
            return tax_dict['parent_taxid']
        else:
            sys.stderr.write("Warning: taxid {} not found\n".format(taxid))
            return None

    def get_organism_parent(self, taxid):
        p_taxid = self.get_organism_parent_taxid(taxid)
        if p_taxid == None:
            return None
        return self.taxonomy_coll.find_one({'taxid':p_taxid})

    def get_family(self, taxid):
        return self.get_rank(taxid, 'family')

    def get_rank(self, taxid, rank, field_to_return = None):
        parent = self.get_organism_parent(taxid)
        if parent == None:
            return None
        while (parent['rank'] != rank and parent['parent_taxid'] != None
                and parent['parent_taxid'] != parent['taxid']):
            parent = self.get_organism_parent(parent['taxid'])
        if parent['rank'] == rank:
            if field_to_return:
                return parent[field_to_return]
            return parent
        else:
            return None

    def get_rank_lineage(self, taxid):
        """Get lineage only including ranks that are not 'no rank'
        Return as a dictionary containing {'rank': taxid}
        """
        rank = self.taxid_to_taxonomy(taxid)['rank']
        lineage = {}
        if rank != 'no rank':
            lineage = {rank: taxid}
        parent = self.get_organism_parent(taxid)
        while parent['taxid'] != 1:
            rank = parent['rank']
            if rank != 'no rank':
                lineage[rank] = parent['taxid']
            parent = self.get_organism_parent(parent['taxid'])
        if rank != 'no rank':
            lineage[rank] = parent['taxid']
        return lineage

    def get_children_taxid(self, taxid):
        if type(taxid) is int:
            taxid = [taxid]
        return [x['taxid'] for x in self.taxonomy_coll.find({'parent_taxid': {'$in': taxid}}, {'taxid': 1})]

    def get_all_children_taxid(self, taxid, has_gene = False, include_self = True):
        if has_gene:
            children = self.taxonomy_coll.find({'lineage': taxid, 'has_gene': has_gene})
        else:
            children = self.taxonomy_coll.find({'lineage': taxid})
        if include_self:
            return [taxid] + [c['taxid'] for c in children]
        else:
            return [c['taxid'] for c in children]

    def organism_to_taxonomy(self, organism, _search_sub = True):
        # example: 'Krokinobacter sp. (strain 4H-3-7-5)'
        # example: 'Clostera anastomosis granulovirus M3246'
        if organism == None:
            return None
        organism = organism.lower().strip()
        match = list(self.taxonomy_coll.find({'scientific_name': organism}))
        if len(match) > 1:
            print('More than one match for scientific name: ' + organism + '. Matches: ' + str(match), file=sys.stderr)
        if len(match) >= 1:
            return match[0]
        match = list(self.taxonomy_coll.find({'uniprot_name': organism}))
        if len(match) > 1:
            print('More than one match for uniprot name: ' + organism + '. Matches: ' + str(match), file=sys.stderr)
        if len(match) >= 1:
            return match[0]
        match = list(self.taxonomy_coll.find({'synonyms': organism}))
        if len(match) > 1:
            print('More than one match for synonyms:' + organism + '. Matches: ' + str(match), file=sys.stderr)
        if len(match) >= 1:
            return match[0]

        # Remove words from the end one by one and look up resulting string
        # Example: 'Clostera anastomosis granulovirus M3246' -> 'Clostera anastomosis granulovirus' -> 'Clostera anastomosis'
        if _search_sub:
            print('Search for substrings of: ' + organism, file=sys.stderr)
            organism_cut = ' '.join(organism.split(' ')[:-1])
            while organism_cut:
                tax = self.organism_to_taxonomy(organism_cut, _search_sub = False)
                if tax:
                    print('Matched: ' + organism + ' to: ' + organism_cut, file=sys.stderr)
                    return tax
                organism_cut = ' '.join(organism_cut.split(' ')[:-1])

        print('No taxonomy match: ' + organism, file=sys.stderr)
        return None

    def taxid_to_taxonomy(self, taxid):
        # accepts a tax_id or list of tax_ids
        # returns taxonomy dicts. Not neccesarily in the same order as input
        if type(taxid) == int:
            return self.taxonomy_coll.find_one({'taxid': taxid})
        if type(taxid) == str:
            return None
        taxid = [t for t in taxid if t]
        return list(self.taxonomy_coll.find({'taxid': {'$in': taxid}}))

    def get_lineages(self, tax_ids):
        # Skip taxids that are None, remove lineages that are None
        # Therefore, list of lineages may not be the same length as list of tax_ids
        taxonomy = self.taxonomy_coll
        lineages = []
        tax_ids = [x for x in tax_ids if x is not None]

        found_taxids, lineages = zip(*[(x['taxid'],x['lineage']) for x in taxonomy.find({'taxid': {'$in': tax_ids}}, {'lineage':True, 'taxid': True})])
        missed_taxids = set(tax_ids) - set(found_taxids)
        for i in missed_taxids:
            print('Failed to find taxonomy entry for: ' + str(i))

        return lineages

    def LCA(self, tax_ids):
        # Get LCA using all nodes
        # Accepts list of taxIDs
        tax_ids = [x for x in tax_ids if x]
        if not tax_ids:
            return None
        if len(tax_ids) == 1:
            return tax_ids[0]
        lineages = self.get_lineages(tax_ids)
        if lineages == [None]:
            return None
        common_ancestors = set(lineages[0]).intersection(*lineages[1:])
        for id in lineages[0]:
            if id in common_ancestors:
                return id

    def LCA_rank(self, tax_ids):
        # Get LCA using only nodes with a rank
        #tax_ids = [457421, 457412, 123]
        lineages=[]
        for taxid in tax_ids:
            lineages.append(self.get_rank_lineage(taxid))

        for rank in self.ranks:
            matches = set([x.get(rank, None) for x in lineages])
            if len(matches) == 1 and matches != {None}:
                return (rank, lineages[0][rank])


    def demo(self):
        taxid = 639282
        print('Get entry for the genus of organism: ' + str(taxid))
        print(self.taxid_to_taxonomy(taxid)[0]['scientific name'])
        genus = self.get_rank(taxid, 'genus')
        print('Get all species and strains within this genus')
        child_ids = self.get_all_children_taxid(genus['taxid'])
        print([x['scientific name'] for x in self.taxid_to_taxonomy(child_ids)])

        # takes too long
        tax = self.get_all_children_taxid(7776)
        print('done')


def fasta_to_organism_uniprot(fasta_defline):
    g = re.search('OS=(.*?)( [A-Z]{2}=|\|UniProt)', fasta_defline)
    # Matches everything after 'OS=' until either : (any two uppercase letters, = ) or (|UniProt)
    if g:
        return g.groups()[0]
    else:
        #logger.warn('Unmatched organism' + fasta_defline)
        return None


def fasta_to_organism_refseq(fasta_defline):
    '''
    refseq defline are retarded and have no standardized way of noting organism with square brackets
    testers:
    string = "[GSEE] tandem repeats [Invertebrate iridescent virus 30]"
    string = "coat protein [Euphorbia mosaic virus - A [Mexico:Yucatan:2004]]" #pID: 745
    string = "[citrate [pro-3S]-lyase] ligase [Vibrio cholerae]"
    '''
    # If there are 8 pipes (coming from fasta file)
    if fasta_defline.count('|') == 8:
        # Split defline by '|'. Take the 6th
        txt = fasta_defline.split('|')[6].strip()
    elif fasta_defline.count('|') == 4:
        txt = fasta_defline.split('|')[-1].strip()
    else:
        txt = fasta_defline

    brackets = list(parse_brackets(txt))
    if not brackets:
        #logger.warn('Malformed defline: ' + fasta_defline)
        try:     # Just return whatever is between the last '[' and the last ']'
            txt_flip = txt[::-1]
            organism = txt_flip[txt_flip.find(']')+1:txt_flip.find('[')][::-1]
            if len(organism) > 4:
                return organism
            else:
                return None
        except Exception:
            return None

    # if there are nested brackets
    if max(list(zip(*brackets))[0]) > 0:
        # take the string at level 0
        organism = [x[1] for x in brackets if x[0] == 0]
        if len(organism) > 1:
            # take the last bracket at level 0
            #logger.warn('Confusing defline: ' + fasta_defline)
            return organism[-1]
        else:
            return organism[0]
    else:
        # there are no nested brackets
        # take the last brackets
        return brackets[-1][1]

def parse_brackets(string):
    """Generate parenthesized contents in string as pairs (level, contents).
    http://stackoverflow.com/questions/4284991/parsing-nested-parentheses-in-python-grab-content-by-level
    """
    if string.count('[') != string.count(']'):
        return None

    stack = []
    for i, c in enumerate(string):
        if c == '[':
            stack.append(i)
        elif c == ']' and stack:
            start = stack.pop()
            yield (len(stack), string[start + 1: i])

def fasta_to_organism(fasta_defline):
    # if passing a dictionary from ProtDB
    if type(fasta_defline) is dict:
        if 'r' in fasta_defline:
            if fasta_defline['r'].lower().count('refseq') or fasta_defline['r'].lower().count('hmp_reference_genomes'):
                return fasta_to_organism_refseq(fasta_defline['d'])
            elif fasta_defline['r'].lower().count('uniprot'):
                return fasta_to_organism_uniprot(fasta_defline['d'])
        print('No db info: ' + str(fasta_defline), file = sys.stderr)
    # If passing a string from a fasta file
    else:
        linesplit = fasta_defline.split('|')
        repository = linesplit[-2].lower()
        if repository == 'refseq' or repository == 'hmp_reference_genomes':
            return fasta_to_organism_refseq(linesplit[-3])
        elif repository.count('uniprot'):
            return fasta_to_organism_uniprot(fasta_defline)
        print('No db info: ' + fasta_defline, file = sys.stderr)
    return None
