'''
See https://github.com/tanghaibao/goatools
Modified from commit: https://github.com/tanghaibao/goatools/commit/540935f8bfe1c1278c76a1b3032c099c1cd09bae
Parses "relationship: part_of" and adds to parents
'''

#!/usr/bin/env python
# -*- coding: UTF-8 -*-
from __future__ import print_function
import sys
import os
try:
    from exceptions import EOFError  # python3
except ImportError:
    pass # in default namespace in python2 

typedef_tag, term_tag = "[Typedef]", "[Term]"
short_namespace = {'biological_process':'BP',
                   'molecular_function': 'MF',
                   'cellular_component': 'CC'}


def after_colon(line):
    # macro for getting anything after the :
    return line.split(":", 1)[1].strip()


def read_until(handle, start):
    # read each line until it has a certain start, and then puts
    # the start tag back
    while 1:
        pos = handle.tell()
        line = handle.readline()
        if not line:
            break
        if line.startswith(start):
            handle.seek(pos)
            return
    raise EOFError("%s tag cannot be found" % start)


class OBOReader:
    """
    parse obo file, usually the most updated can be downloaded from
    http://purl.obolibrary.org/obo/go/go-basic.obo

    >>> reader = OBOReader()
    >>> for rec in reader:
            print rec

    """

    def __init__(self, obo_file="go-basic.obo"):

        try:
            self._handle = open(obo_file)
        except:
            print(("download obo file first\n "
                                 "[http://purl.obolibrary.org/obo/"
                                 "go/go-basic.obo]"), file=sys.stderr)
            sys.exit(1)

    def __iter__(self):

        line = self._handle.readline()
        if not line.startswith(term_tag):
            read_until(self._handle, term_tag)
        while 1:
            yield self.__next__()

    def __next__(self):

        lines = []
        line = self._handle.readline()
        if not line or line.startswith(typedef_tag):
            raise StopIteration

        # read until the next tag and save everything in between
        while 1:
            pos = self._handle.tell()   # save current postion for roll-back
            line = self._handle.readline()
            if not line or (line.startswith(typedef_tag)
                            or line.startswith(term_tag)):
                self._handle.seek(pos)  # roll-back
                break
            lines.append(line)

        rec = GOTerm()
        for line in lines:
            if line.startswith("id:"):
                rec.id = after_colon(line)
            if line.startswith("alt_id:"):
                rec.alt_ids.append(after_colon(line))
            elif line.startswith("name:"):
                rec.name = after_colon(line)
            elif line.startswith("namespace:"):
                rec.namespace = after_colon(line)
            elif line.startswith("is_a:"):
                rec._parents.append(after_colon(line).split()[0])
            elif (line.startswith("is_obsolete:") and
                  after_colon(line) == "true"):
                rec.is_obsolete = True
            elif line.startswith("relationship: part_of"):
                rec._parents.append(after_colon(line).split()[1].strip())

        return rec


class GOTerm:
    """
    GO term, actually contain a lot more properties than interfaced here
    """

    def __init__(self):
        self.id = ""                # GO:xxxxxx
        self.name = ""              # description
        self.namespace = ""         # BP, CC, MF
        self._parents = []          # is_a basestring of parents
        self.parents = []           # parent records
        self.children = []          # children records
        self.level = -1             # distance from root node
        self.is_obsolete = False    # is_obsolete
        self.alt_ids = []           # alternative identifiers

    def __str__(self):
        return "{id}: {name} [{namespace}]".format(id=self.id, name=self.name, namespace=short_namespace[self.namespace])

    def __repr__(self):
        return "GOTerm('%s')" % (self.id)

    def has_parent(self, term):
        for p in self.parents:
            if p.id == term or p.has_parent(term):
                return True
        return False

    def has_child(self, term):
        for p in self.children:
            if p.id == term or p.has_child(term):
                return True
        return False

    def get_all_parents(self):
        all_parents = set()
        for p in self.parents:
            all_parents.add(p.id)
            all_parents |= p.get_all_parents()
        return all_parents

    def get_all_children(self):
        all_children = set()
        for p in self.children:
            all_children.add(p.id)
            all_children |= p.get_all_children()
        return all_children

    def get_all_parent_edges(self):
        all_parent_edges = set()
        for p in self.parents:
            all_parent_edges.add((self.id, p.id))
            all_parent_edges |= p.get_all_parent_edges()
        return all_parent_edges

    def get_all_child_edges(self):
        all_child_edges = set()
        for p in self.children:
            all_child_edges.add((p.id, self.id))
            all_child_edges |= p.get_all_child_edges()
        return all_child_edges


class GODag(dict):

    def __init__(self, obo_file = None):
        if not obo_file:
            obo_file = os.path.join(os.path.dirname(__file__), 'go-basic.obo')
        self.load_obo_file(obo_file)

    def load_obo_file(self, obo_file):

        print("load obo file %s" % obo_file, file=sys.stderr)
        obo_reader = OBOReader(obo_file)
        for rec in obo_reader:
            self[rec.id] = rec
            for alt in rec.alt_ids:
                self[alt] = rec

        self.populate_terms()
        print(len(self), "nodes imported", file=sys.stderr)

    def populate_terms(self):

        def depth(rec):
            if rec.level < 0:
                if not rec.parents:
                    rec.level = 0
                else:
                    rec.level = min(depth(rec) for rec in rec.parents) + 1
            return rec.level

        # make the parents references to the GO terms
        for rec in self.values():
            rec.parents = [self[x] for x in rec._parents]

        # populate children and levels
        for rec in self.values():
            for p in rec.parents:
                p.children.append(rec)

            if rec.level < 0:
                depth(rec)

    def write_dag(self, out=sys.stdout):
        for rec_id, rec in sorted(self.items()):
            print(rec, file=out)

    def query_term(self, term, verbose=False):
        if term not in self:
            print("Term %s not found!" % term, file=sys.stderr)
            return

        rec = self[term]
        print(rec, file=sys.stderr)
        if verbose:
            print("all parents:", rec.get_all_parents(), file=sys.stderr)
            print("all children:", rec.get_all_children(), file=sys.stderr)

        return rec

    def paths_to_top(self, term, verbose=False):
        """ Returns all possible paths to the root node

            Each path includes the term given. The order of the path is
            top -> bottom, i.e. it starts with the root and ends with the
            given term (inclusively).

            Parameters:
            -----------
            - term:
                the id of the GO term, where the paths begin (i.e. the
                accession 'GO:0003682')

            Returns:
            --------
            - a list of lists of GO Terms
        """
        # error handling consistent with original authors
        if term not in self:
            print("Term %s not found!" % term, file=sys.stderr)
            return

        def _paths_to_top_recursive(rec):
            if rec.level == 0:
                return [[rec]]
            paths = []
            for parent in rec.parents:
                top_paths = _paths_to_top_recursive(parent)
                for top_path in top_paths:
                    top_path.append(rec)
                    paths.append(top_path)
            return paths

        go_term = self[term]
        return _paths_to_top_recursive(go_term)

    def update_association(self, association):
        bad_terms = set()
        for key, terms in list(association.items()):
            parents = set()
            for term in terms:
                try:
                    parents.update(self[term].get_all_parents())
                except:
                    bad_terms.add(term.strip())
            terms.update(parents)
        if bad_terms:
            print("terms not found: %s" % (bad_terms,), file=sys.stderr)

