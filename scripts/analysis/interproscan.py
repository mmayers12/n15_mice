"""
For representing interproscan terms
attribute: name, accesion ID, reference to parent and children terms, level #, and type

Flat files:
ftp://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt
ftp://ftp.ebi.ac.uk/pub/databases/interpro/entry.list

@author: gstupp
interproscan.py

"""
import os
import sys

class IPRTerm():
    # Represents one interproscan term
    def __init__(self, name='', id='', parent=None, level=0, type=None):
        self.name = name
        self.ipa = id # maintain for compatability, equals id
        self.id = id
        self.parent = parent
        self.children = []
        self.level = level
        self.type = type

    def __repr__(self):
        return self.ipa + ': ' + self.name

    def __str__(self):
        return self.ipa + ': ' + self.name

class InterproscanParser():
    def __init__(self):
        try:
            self._handleTree = open(os.path.join(os.path.dirname(__file__), 'ParentChildTreeFile.txt'))
            self._handleEntry = open(os.path.join(os.path.dirname(__file__), 'entry.list'))
        except:
            print("download flat files first\n" +
            "ftp://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt\n" +
            "ftp://ftp.ebi.ac.uk/pub/databases/interpro/entry.list", file=sys.stderr)
            sys.exit(1)
        self.d = {}
        self.parse()
        self.parse_ipa_name()

    def parse_ipa_name(self):
        """ Looks like:
        Active_site
        IPR000126 Serine proteases, V8 family, serine active site
        IPR000138 Hydroxymethylglutaryl-CoA lyase, active site
        Binding_site
        IPR000048 IQ motif, EF-hand binding site
        IPR000222 PPM-type phosphatase, divalent cation binding
        """
        ipr_type = None
        for line in self._handleEntry:
            if not line.strip(): continue
            line = line.strip().split(' ',1)
            if len(line) == 1:
                # this is a `type`
                ipr_type = line[0]
                continue
            id = line[0]
            term = IPRTerm(id=id, name=line[1], type=ipr_type, level=0)
            if id in self.d:
                self.d[id].name = term.name
                self.d[id].type = term.type
            else:
                self.d[id] = term

    def parse(self):
        prev_level = 0
        prev_parents = []

        for line in self._handleTree:
            line_split = line.strip().split('::')
            id = line_split[0][line_split[0].index('I'):]
            curr_level = int(line_split[0].count('-') / 2)

            if curr_level == 0:
                term = IPRTerm(id=id, name=line_split[1], parent=None, level=0)
                prev_parents = [term]

            elif curr_level > prev_level:
                term = IPRTerm(id=id, name=line_split[1], parent=prev_parents[-1], level=curr_level)
                prev_parents.append(term)

            elif curr_level < prev_level:
                for i in range(prev_level + 1 - curr_level): prev_parents.pop()
                term = IPRTerm(id=id, name=line_split[1], parent=prev_parents[-1], level=curr_level)
                prev_parents.append(term)

            elif curr_level == prev_level:
                prev_parents.pop()
                term = IPRTerm(id=id, name=line_split[1], parent=prev_parents[-1], level=curr_level)
                prev_parents.append(term)

            prev_level = curr_level
            self.d[term.id] = term

        # populate children
        for rec in self.d.values():
            if rec.parent:
                # this works because a child can only have one parent
                # and IPRTerm.parent is a reference to the IPRTerm object
                rec.parent.children.append(rec)

class Interproscan(dict):
    # object that handles parsing flat file and doing lookups
    def __init__(self):
        # load from flat files
        parser = InterproscanParser()
        for key, value in parser.d.items():
            self[key] = value


