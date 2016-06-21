"""
SQT File Parser
DTASelect-filter.txt File Parser

"""
import re
from itertools import chain
from collections import defaultdict
from .utils import get_lcstep, get_unmod_seq


def get_peptides(dta_select, **kwargs):
    """ Return set of peptides in this dtaselect-filter file
    :param dta_select: list or iterator of Loci
    :return: set.
    """
    return set(chain(*[x['peptide_seq'] for x in dta_select]))
    
def get_unmod_peptides(dta_select, **kwargs):
    """ Return set of peptides in this dtaselect-filter file
    :param dta_select: list or iterator of Loci
    :return: set.
    """
    return set(chain(*[x['unmod_peptide_seq'] for x in dta_select]))

def dta_select_header(in_file):
    """
    Just grab the header from a DTASelect-filter.txt file file.
    --- Example:
    DTASelect v2.1.3
    /mongoc/gstupp/DTASelect/2015_10_23_CM1E2_d56col_unenr/RagcoH
    /mongoc/gstupp/DTASelect/2015_10_23_CM1E2_d56col_unenr/RagcoH/output.fasta
    Blazmass ? in SQT format.
    --quiet --brief --trypstat --modstat -y 2 -DM 10 --extra --dm --pfp 0.01 -p 2
    ---
    """
    header = {}
    with open(in_file) as f:
        header['dta_select_version'] = next(f).strip()
        header['folder_path'] = next(f).strip()
        header['fasta_db'] = next(f).strip()
        header['search_type'] = next(f).strip()
        header['dta_select_params'] = next(f).strip()
    header['file_name'] = in_file
    return header

def dta_select_parser(in_file, small=False, get_tax=False, check_peptides=False, get_hashes=False, 
                      return_reverse=True, taxDB=None, protDB=None):
    """ steps through and parses a DTASelect-filter.txt file (generator function)
    :param in_file: path to DTASelect-filter.txt file
    :param small: get rid of keys: 'loci' and 'peptides'
    :param get_tax: Look up taxonomy information for each protDB ID
    :param check_peptides: Check if all peptide sequences are in all proteins sequences in 'forward_loci'
    :param get_hashes:
    :param return_reverse: include reverse loci
    :pararm taxDB: mongo collection to the taxDB (containing protdbID -> taxid mapping) or None
    :param protDB: mongo collection to protDB (protdbID information). If protDB is given, a better protein['name'] will be attempted
    :type in_file: str

        get_forward_loci is removed
        noProtDB is removed. If the locus contains a "||", the part before the pipes is determined to be the protID

    protein is a dict with (possible) fields:
        loci: list of locus dicts. described below
        peptides: list of peptide dicts. described below
        reverse: boolean. True if all loci for a protein are reverse
        peptide_seq: set. All peptide amino acid sequences
        forward_loci: list. 'Locus' fields in loci for forward loci
        all_loci: list. 'Locus' fields in loci
        name: string. The "representative" locus (the largest one by AA length)
        tax_id: list. Unique list of taxonomy IDs for all forward loci
        lca: Int or None. Lowest common ancestor of tax_ids
        hashes: list. MD5sums of protein sequences for forward_loci. len(hashes) gives the number of unique proteins matching

    loci: dict parsed from Loci lines in file. Mostly unchanged
        fields from file: Locus, Sequence Count, Spectrum Count, Sequence Coverage, Length, MolWt, pI, Validation Status, NSAF, EMPAI, Descriptive Name
        fields added:
            reverse: boolean. True if `locus` starts with "Reverse_"
            description: part after'||', if exists
    peptides: dict parsed from peptide lines in file. Mostly unchanged.
        fields from file: Unique, FileName, XCorr, DeltCN, Conf%, M+H+m CalcM+H+, TotalIntensity, SpR, SpScore, IonProportion, Redundancy, Sequence
        fields added:
            aa_sequence: `Sequence` with the left and right sequence stripped
            is_modified: boolean. True if the peptide has PTMs
            unmod_peptide: peptide sequence without PTMs
            diff_mass: mass difference of PTMs
            mods: list of tuples: (AA (amino acid that is modified), pos (1-based position within peptide), mass (mass of this PTM))
            lc_step, scan, charge_state:  parsed from `FileName`

    """

    def per_to_float(x):
        # '50%' -> 0.50
        return float(x[:-1]) / 100

    if get_tax:
        from ..analysis import taxonomy
        from pymongo import MongoClient
        if not taxDB:
            taxDB = MongoClient('wl-cmadmin.scripps.edu', 27017).TaxDB_072114.TaxDB_072114
        t = taxonomy.Taxonomy(host='wl-cmadmin.scripps.edu', port=27017)

    if check_peptides:
        protDB = MongoClient('wl-cmadmin.scripps.edu', 27018).ProtDB_072114.ProtDB_072114

    if get_hashes:
        client = MongoClient('wl-cmadmin.scripps.edu', 27017)
        redunDB = client.redunDB.redunDB
        
        
    locus_types = [str, int, int, per_to_float, int, int, float, str, float, float, str]
    peptide_types = [str, str, float, float, float, float, float, float, int, float, float, int, str]
    peptide_types_ppm = [str, str, float, float, float, float, float, float, float, int, float, float, int, str]
    with open(in_file) as f:
        line = next(f)

        # Skip header.
        # TODO: parse it
        while not line.startswith('Locus'):
            line = next(f)
            if "--dm" in line:
                peptide_types = peptide_types_ppm

        # Read locus header
        locus_columns = line.strip().split('\t')

        # Read peptide header
        line = next(f)
        peptide_columns = line.rstrip().split('\t')
        line = next(f)

        # Read the rest of the file
        while line != '\tProteins\tPeptide IDs\tSpectra\n':
            protein = dict()
            # Read loci for a protein
            loci = []
            while not line.startswith('\t') and not line.startswith('*'):  # While it starts with a number or "Rev'
                try:
                    loci.append(dict(zip(locus_columns, [x(y) for x, y in zip(locus_types, line.strip().split('\t'))])))
                except:
                    print("parsing line failed: " + line)
                line = next(f)

            # Read peptides
            peptides = []
            while line[0] in ['\t', '*'] and line != '\tProteins\tPeptide IDs\tSpectra\n':
                peptides.append(
                    dict(zip(peptide_columns, [x(y) for x, y in zip(peptide_types, line.rstrip().split('\t'))])))
                line = next(f)

            protein['peptides'] = peptides
            protein['loci'] = loci

            # Parse out reverse loci
            for l in protein['loci']:
                if l['Locus'].startswith('Reverse_'):
                    l['reverse'] = True
                    l['Locus'] = l['Locus'][8:]
                else:
                    l['reverse'] = False
                if l['Locus'].count("||"):
                    locus_split = l['Locus'].split('||')
                    l['description'] = '||'.join(locus_split[1:])
                    l['Locus'] = int(locus_split[0])
                else:
                    l['description'] = l['Locus']
                    l['Locus'] = int(l['Locus'])

            # Are all loci for a protein reverse?
            ## - shouldn't this logic be changed to: 'are any loci reverse?' and if so, set 'Reverse' to True?
            ## (because then we can't distinguish this peptide match from a fictional protein from a real protein)
            ## - There aren't typically many overlapping peptides between forward and reverse proteins anyway (~1%)
            if all([l['reverse'] for x in protein['loci']]):
                protein['reverse'] = True
            else:
                protein['reverse'] = False

            if not return_reverse and protein['reverse']:  # skip this one if we are skipping reverse loci
                continue

            # Pull out peptide sequences
            for p in protein['peptides']:
                p['aa_sequence'] = re.findall('\.(.*)\.', p['Sequence'])[0]
                p['is_modified'] = True if ')' in p['aa_sequence'] else False
                if p['is_modified']:
                    seq, mods = get_unmod_seq(p['aa_sequence'])                 
                    p.update(mods)
                else:
                    p['unmod_peptide'] = p['aa_sequence']

                # MudPIT salt step (chromatography method from Xcalibur)
                p['lc_step'] = get_lcstep(p['FileName'])

                # Scan number from instrument (unique per salt step) - from MS2 / SQT file
                p['scan'] = int(p['FileName'].split('.')[1])

                # predicted ion charge from instrument - from MS2 / SQT file
                p['charge_state'] = int(p['FileName'].split('.')[3])

                # To try to not break things
                p['LCStep'] = p['lc_step']
                p['Scan'] = p['scan']
                p['ChargeState'] = p['charge_state']
                p['AA_Sequence'] = p['aa_sequence']
                p['isModified'] = p['is_modified']

            protein['peptide_seq'] = list(set((x['aa_sequence'] for x in protein['peptides'])))
            protein['unmod_peptide_seq'] = list(set((x['unmod_peptide'] for x in protein['peptides'])))
            protein['quantification'] = sum([x['Redundancy'] for x in protein['peptides']])
            protein['forward_loci'] = [l['Locus'] for l in protein['loci'] if not l['reverse']]
            protein['all_loci'] = [l['Locus'] for l in protein['loci']]
            
            def is_good_db(s):
                s = s.lower()
                if "refseq" in s or "uniprot" in s or s=="hmp_reference_genomes":
                    return True
                else:
                    return False

            # get a "representative" locus
            if protDB:
                # pick the largest protein in one of dbs: ['RefSeq','UniProt*', 'HMP_Reference_Genomes']
                p_result = [x for x in protDB.find({'_id':{'$in':protein['forward_loci']}}) if is_good_db(x['r'])]
                if p_result:
                    protein['name'] = max([(len(p['s']),p['d']) for p in p_result], key=lambda x:x[0])[1]
            if 'name' not in protein:                    
                # the largest one in any db
                max_length = max(l['Length'] for l in protein['loci'])
                protein['name'] = [l['description'] for l in protein['loci'] if l['Length'] == max_length][0]
            
            if get_tax:
                # get all possible taxIDs
                protDB_ids = protein['forward_loci']
                assert all(isinstance(x, int) for x in protDB_ids)
                taxIDs_doc = list(taxDB.aggregate(
                    [{'$match': {'_id': {'$in': protDB_ids}}},
                     {'$group': {'_id': None, 'taxid': {'$addToSet': '$taxid'}}}]))
                if taxIDs_doc:
                    protein['tax_id'] = [x for x in taxIDs_doc[0]['taxid'] if x]
                    protein['lca'] = t.LCA(taxIDs_doc[0]['taxid'])
                else:
                    protein['tax_id'] = []
                    protein['lca'] = None
                # To try to not break things
                protein['LCA'] = protein['lca']

            if check_peptides:
                # Are all peptides found within the fasta sequences for all possible forward_loci ?
                # Skip reverse loci. May want to change this to use all loci, regardless of forward or reverse
                # to avoid some proteins not having these entries.
                # Keeping like this for now to keep compatibility with get_forward_loci lookup
                if protein['forward_loci']:
                    protein['protDB'] = list(protDB.find({'_id': {'$in': protein['forward_loci']}}))
                    defline, seq = zip(*[(x['d'], x['s']) for x in protein['protDB']])
                    protein['all_peptides_in_proteins'] = all(
                        [all([p in s for p in protein['peptide_seq']]) for s in seq])
                    if not protein['all_peptides_in_proteins']:
                        print('not all peptides in proteins' + str(protein['forward_loci'][0]))

            if get_hashes:
                protein['hashes'] = [x['_id'] for x in redunDB.find({'pID': {'$in': protein['forward_loci']}})]

            if small:
                protein['descriptions'] = [l['description'] for l in protein['loci']]
                del protein['loci']
                del protein['peptides']

            yield protein


def sqt_parser_demo():
    import matplotlib.pyplot as plt
    filepath = "/mongoa/DTASelect/10_2014_mass_spec/102414_Bfragilis_HCD35/indexDB_search_noProtDB/102414_SC_Bfragilis_25ug_HCD_FTMS_MS2_05.sqt"
    Xcorr, deltCN, Xcorr_r, deltCN_r = ([] for i in range(4))
    for chunk in sqt_chunker(filepath):
        if chunk['matches'][0]['reverse']:
            Xcorr_r.append(chunk['Xcorr'])
            deltCN_r.append(chunk['deltCN'])
        else:
            Xcorr.appen(chunk['Xcorr'])
            deltCN.appen(chunk['deltCN'])
    fig = plt.figure(figsize=(12, 8))
    fig.gca()
    ax = fig.add_subplot(111)
    ax.scatter(Xcorr, deltCN, alpha=.5, marker='.', color='blue')
    ax.scatter(Xcorr_r, deltCN_r, alpha=0.5, marker='.', color='red')


def sqt_chunker(sqt_filename):
    """ yield "chunks" of the MS2 file, with one chunk per scan
    (generator function) """

    current_chunk = []
    with open(sqt_filename) as f:
        for line in f:
            if line[0] != 'H':
                if line[:2] == 'S\t' and current_chunk:
                    # Some spectra have no matches. Skip them
                    if len(current_chunk) > 1:
                        yield parse_sqt_chunk(current_chunk)
                    current_chunk = [line]  # reset current block
                else:
                    current_chunk.append(line)
        else:
            yield parse_sqt_chunk(current_chunk)


def parse_sqt_chunk(chunk):
    # Parse spectra and matches
    # Define column names and types
    spectrum_keys = ['low_scan', 'high_scan', 'charge', 'process_time', 'server', 'obs_mass', 'total_intensity',
                     'lowest_SP', '#seq_matching', 'top_match']
    spectrum_type = [int, int, int, int, str, float, float, float, int]
    match_keys = ['rank_Xcorr', 'rank_Sp', 'calc_mass', 'deltCN', 'Xcorr', 'Zscore', 'matched_ions', 'expected_ions',
                  'seq_matched', 'valid_status']
    match_type = [int, int, float, float, float, float, int, int, str, str]
    # Accept a list of raw strings. Split them and turn into an iterator that yields lines
    chunk = [l.split() for l in chunk]
    chunk = iter(chunk)
    line = next(chunk)
    # First line is Spectrum line
    s = [x(y) for x, y in zip(spectrum_type, line[1:])]
    spect_dict = dict(zip(spectrum_keys, s))
    spect_dict['matches'] = []
    while True:
        if line[0] == 'M':
            m = [x(y) for x, y in zip(match_type, line[1:])]
            match_dict = dict(zip(match_keys, m))
            match_dict['L'] = []
            match_dict['reverse'] = False  # if any L within a M is Reverse, the whole M is reverse
            while True:
                try:
                    line = next(chunk)
                except StopIteration:
                    break
                if line[0] == 'L':
                    # Sometimes, the L line is blank for no apparent reason
                    if len(line) > 1:
                        if line[1].startswith('Reverse'):
                            match_dict['reverse'] = True
                            L = int(line[1][line[1].index('_') + 1:])
                        else:
                            L = int(line[1])
                        match_dict['L'].append(L)
                    else:
                        match_dict['L'].append(None)
                else:
                    spect_dict['matches'].append(match_dict)
                    break
        else:
            try:
                line = next(chunk)
            except StopIteration:
                spect_dict['matches'].append(match_dict)
                break
    # Give delt_CN and Xcorr for the spectrum
    spect_dict['Xcorr'] = spect_dict['matches'][0]['Xcorr']
    if len(spect_dict['matches']) > 1:
        spect_dict['deltCN'] = spect_dict['matches'][1]['deltCN']
    else:
        spect_dict['deltCN'] = 0

    return spect_dict


def build_pep_quant_dict(dtaselect_filter, field='Redundancy', **kwargs):
    """
    From a dtaselect-filter file, return a peptide quantification dict
    return {peptide_seq:field}
    output example: {'AAA': 43, 'ADB': 5}

    We can't use the unique list of peptides because dtaselect reports different isotopes of the same peptide
    (in different scans) with different quantifications. Also need to handle PTMs
    Example:
        121514_SC_sampleH1sol_25ug_HCD_FTMS_MS2_04.8403.8403.2  5.722176        0.36840343      100.0   1763.8446       1763.849        88.0    1       22.869139       0.0     3       K.AANLDAIHDTVHEMAR.D
        121514_SC_sampleH1sol_25ug_HCD_FTMS_MS2_04.8402.8402.3  5.1739216       0.17002201      99.9    1763.846        1763.849        88.0    1       15.8237915      0.0     3       K.AANLDAIHDTVHEMAR.D
        121514_SC_sampleH1sol_25ug_HCD_FTMS_MS2_05.6199.6199.2  4.840225        0.42235494      100.0   1779.8467       1779.8439       88.0    1       16.761187       0.0     8       K.AANLDAIHDTVHEM(15.99489974975586)AR.D
        121514_SC_sampleH1sol_25ug_HCD_FTMS_MS2_04.6378.6378.3  4.192363        0.21340305      100.0   1779.8423       1779.8439       88.0    1       10.633918       0.0     5       K.AANLDAIHDTVHEM(15.99489974975586)AR.D
    becomes:
        pep['AANLDAIHDTVHEMAR'] ->  {'121514_SC_sampleH1sol_25ug_HCD_FTMS_MS2_04.8402.8402.3': 3, '121514_SC_sampleH1sol_25ug_HCD_FTMS_MS2_04.8403.8403.2': 3,
                                     '121514_SC_sampleH1sol_25ug_HCD_FTMS_MS2_05.6199.6199.2': 8, '121514_SC_sampleH1sol_25ug_HCD_FTMS_MS2_04.6378.6378.3': 5}
    sum them together
    becomes: pep['AANLDAIHDTVHEMAR'] -> 19

    :param dtaselect_filter: if type=str, path to dtaselect_filter file.
        otherwise, assumes is a list or iterator of loci
    :param field: the field to use as the quantification field
    :param kwargs: passed to dta_select_parser if `dtaselect_filter` is a path
    :return: dict.

    """
    if isinstance(dtaselect_filter, str):
        path_to_dtaselect_filter = dtaselect_filter
        parser = dta_select_parser(path_to_dtaselect_filter, **kwargs)
    else:
        parser = dtaselect_filter
    pep = defaultdict(dict)
    for x in chain(*[locus['peptides'] for locus in parser]):
        pep[x.get('unmod_peptide', x['aa_sequence'])][x['FileName']] = x[field]
    pep_quant = dict()
    for peptide,fname_dict in pep.items():
        pep_quant[peptide] = sum(fname_dict.values())
    return pep_quant


