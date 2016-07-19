from collections import defaultdict

from .. import obo_parser
from . import interproscan


class GSEA:
    def __init__(self):
        self.go_ontology = obo_parser.GODag()
        self.ipr_terms = interproscan.Interproscan()

    def generate_go_sets(self, locus_go):
        ## Generate go-term locus-sets.
        # Input: Dict of locus IDs (keys), where a locus ID is a unique ID for that set of PSMs
        #        Where each locus contains its set of go terms (values)
        # Output: Dict of go-terms (keys), and associated loci (values)
        #
        # Go terms are expanded up through their ancestors through "Is a" and "part of" relationships
        # For example a locus annotated at level 3 is part of sets for level 2 and 1 also
        # Trying to emulate msigdb
        # http://www.broadinstitute.org/gsea/msigdb/collection_details.jsp#C5
        # GO gene sets for very broad categories, such as Biological Process, have been omitted from MSigDB.
        # GO gene sets with fewer than 10 genes have also been omitted.
        # Gene sets with the same members have been resolved based on the GO tree structure:
        # if a parent term has only one child term and their gene sets have the same members, the child gene set is omitted;
        # if the gene sets of sibling terms have the same members, the sibling gene sets are omitted.

        """
        #Example Input:
        from random import shuffle
        from itertools import chain
        go_terms = [x.id for x in ontology.values()]
        shuffle(go_terms)
        locus_go = {x:set(go_terms[x:x+3]) for x in range(100)}
        """

        # Make dict of go_term (key) : IDs (value)
        go_locus = defaultdict(set)
        for key, value in locus_go.items():
            for go in value:
                go_locus[go].add(key)

        # Expand up ancestors
        for go_term in list(go_locus.keys())[:]:
            for parent in self.go_ontology[go_term].get_all_parents():
                go_locus[parent] = go_locus[parent].union(go_locus[go_term])

        # Remove "very broad" gene sets. Arbitrary definition: gene sets that emcompass >50% of all IDs
        all_ids = set(locus_go.keys())
        go_locus = {key: value for (key, value) in go_locus.items() if len(value) / len(all_ids) <= 0.5}

        # Remove gene sets with less than 10 members
        go_locus = {key: value for (key, value) in go_locus.items() if len(value) >= 10}

        # Remove child terms with identical gene sets as their parents
        to_remove = set()
        for parent in go_locus.keys():
            # If one child term in go_locus
            cil = [x.id for x in self.go_ontology[parent].children if x.id in go_locus.keys()]
            if len(cil) == 1 and go_locus[cil[0]] == go_locus[parent]:
                to_remove.add(cil[0])
        go_locus = {key: value for (key, value) in go_locus.items() if key not in to_remove}

        # Remove sibling terms with identical gene sets
        to_remove = set()
        for brother in go_locus.keys():
            for parent in self.go_ontology[brother].parents:
                siblings = set([y.id for y in parent.children])
                siblings.remove(brother)
                for sibling in siblings:
                    if sibling in go_locus.keys() and go_locus[brother] == go_locus[sibling]:
                        to_remove.add(sibling)
        go_locus = {key: value for (key, value) in go_locus.items() if key not in to_remove}

        return go_locus

    def generate_ipa_sets(self, locus_ipa):
        ## Generate IPA-term locus-sets.
        # Input: Dict of locus IDs (keys), where a locus ID is a unique ID for that set of PSMs
        #        Where each locus contains its set of IPA terms (values)
        # Output: Dict of IPAs (keys), and associated loci (values)
        #
        # IPA terms are expanded up through their parents
        # For example a locus annotated at level 3 is part of sets for level 2 and 1 also
        # Sets for very broad categories, (making up more than 50% of all loci) are omitted
        # Sets with fewer than 10 proteins have also been omitted.


        """
        #Example Input:
        from random import shuffle
        from itertools import chain
        terms = [x.ipa for x in IPR_terms.values()]
        shuffle(terms)
        locus_ipa = {x:set(terms[x:x+3]) for x in range(100)}
        locus_ipa[0].add('IPR028515')
        locus_ipa[1].add('IPR001972')
        """

        def get_all_parents(term):
            parents = [term.parent]
            for i in range(term.level - 1):
                parents.append(parents[-1].parent)
            return parents if term.parent else []

        def get_children(term, IPR_terms):
            children = [x for x in IPR_terms.values() if term == x.parent]
            return children if children else None

        # Make dict of go_term (key) : IDs (value)
        ipa_locus = defaultdict(set)
        for key, value in locus_ipa.items():
            for ipa in value:
                ipa_locus[ipa].add(key)

        # Expand up ancestors
        for ipa_term in list(ipa_locus.keys())[:]:
            if ipa_term in self.ipr_terms:  # IPR_terms only contains one that have parents. Not all IPR terms are in that file
                for parent in get_all_parents(self.ipr_terms[ipa_term]):
                    ipa_locus[parent.ipa] = ipa_locus[parent.ipa].union(ipa_locus[ipa_term])
        ipa_locus = dict(ipa_locus)

        # Remove "very broad" gene sets. Arbitrary definition: gene sets that emcompass >50% of all IDs
        all_ids = set(locus_ipa.keys())
        ipa_locus = {key: value for (key, value) in ipa_locus.items() if len(value) / len(all_ids) <= 0.5}

        # Remove gene sets with less than 10 members
        ipa_locus = {key: value for (key, value) in ipa_locus.items() if len(value) >= 10}

        # TODO: Remove sets that contains the same sets as their children

        return ipa_locus
