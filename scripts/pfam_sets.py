"""
Build domain architecture db from compil
For each protID. Generate the pfam domain architecture.

Domain architecture is the list of pfam domains in order N-terminus to C-terminus.
Domains entireley contained within another are delimited by a tilde.
Consecutive repeats of a domain are NOT collapsed to a single repeat (at this time)

ArchDB was made in late April eithout removing repeats.
Repeat removal added in May 5th. Need to remake ArchDB


https://www.biostars.org/p/136670/

"""

#%% 
import os
import re
import gzip
import pickle
import subprocess
from pymongo import MongoClient
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

HOST = 'wl-cmadmin'
PORT = 27018
client = MongoClient(HOST, PORT)
hashDB = client.HashDB_072114.HashDB_072114
domainDB = client.DomainDB_072114.DomainDB_072114
protDB = client.ProtDB_072114.ProtDB_072114

def gsave(obj, filename, protocol=3):
    # should check if file exists first...
    with gzip.GzipFile(filename, 'wb') as f:
        pickle.dump(obj, f, protocol)

def gload(filename):
    f = gzip.open(filename,'rb')
    myNewObject = pickle.load(f)
    f.close()
    return myNewObject

def unique_dicts_by_key(dicts, key):
    #Remove duplicate domains with exact same start, stop, sa.
    seen = set()
    return [d for d in dicts if d[key] not in seen and not seen.add(d[key])]

def remove_repeats(domains):
    domains = sorted(domains, key = lambda x: x['st'])
    remove = []
    for idx in range(1, len(domains)):
        if domains[idx]['sa'] == domains[idx-1]['sa']:
            remove.append(idx)
    domains = [x for (idx,x) in enumerate(domains) if idx not in remove]
    return domains
    
def domains_to_str(domains):
    """
    Take list of domain dicts and convert to domain string
    domains = [{'sa': 'A', 'st': 14, 'sp': 300},{'sa': 'A', 'st': 14, 'sp': 300},
     {'sa': 'B', 'st': 100, 'sp': 150},
     {'sa': 'C', 'st': 350, 'sp': 360},
     {'sa': 'D', 'st': 340, 'sp': 400}]
    Returns: 'A~B-D~C'
    """
    # Make id field so we can uniquify the list of domains
    for d in domains:
        d['id'] = ';'.join(map(str,[d['sa'], d['st'], d['sp']]))
    #Remove duplicate domains with exact same start, stop, sa.
    domains = unique_dicts_by_key(domains, 'id')
    # Remove repeats
    domains = remove_repeats(domains)
    # unique list of domains. Sort by start position
    domains = sorted(domains, key = lambda x: x['st'])
    domain_str = ''
    last_stop = 0
    for domain in domains:
        if last_stop == 0:
            delim = ''
        else:
            delim = '-'
        if domain['sp'] <= last_stop:
            domain_str += '~' + domain['sa']
        else:
            domain_str += delim + domain['sa']
            last_stop = domain['sp']
    return domain_str

def process_range(r):
    # process a range of protIDs, so we can parallelize this
    pfam_sets = defaultdict(set)
    ids = {doc['_id']:doc['pID'] for doc in hashDB.find({'pID':{'$in': list(range(r[0],r[1]))}})}
    for doc in domainDB.find({'_id':{'$in': list(ids.keys())}}):
        domain_list = doc['d'].get('Pfam',None)
        if not domain_list:
            continue
        domain_str = domains_to_str(domain_list)
        
        pfam_sets[domain_str].update(set(ids[doc['_id']]))
    return pfam_sets

def merge(dicts):
    # merge a list of dicts by taking the union of the values of each key
    merged = defaultdict(set)
    for d in dicts:
        for k in d:
            merged[k] |= d[k]
    return merged
        
def make_arch_db():
    executor = ProcessPoolExecutor(max_workers=8)
    by = 10000
    m = 60000000
    #by = 2000
    #m = 10000
    e = executor.map(process_range, zip(range(0, m, by),range(by, m+by, by)))
    executor.shutdown()
    print('done calculating architectures')
    pfam_sets = merge(e)
    print(len(pfam_sets))
    gsave(pfam_sets,'pfam_sets.pkl.gz')
    
    # mongodb
    db = MongoClient('wl-cmadmin', 27017).ArchDB_Pfam_071414.ArchDB_Pfam_071414
    db.insert(map(lambda item: {'_id': min(item[1]), 'pID': list(item[1]), 'Pfam': item[0]}, pfam_sets.items()))
    db.ensure_index('pID')
    db.ensure_index('Pfam')

def parse_pfam_clans():
    """
    ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
    Pfam-A.clans.tsv
      This file contains a list of all Pfam-A families that are in clans.
      The columns are: Pfam-A accession, clan accession, clan ID, Pfam-A
      ID, Pfam-A description.
    
    {'accession': Pfam-A accession, 'clan': clan accession, 'clan_name': clan ID, 'name': Pfam-A description}
    """

    pfam = []
    f_path = os.path.join(os.path.dirname(__file__), "Pfam-A.clans.tsv")
    with open(f_path) as f:
        for line in f:
            line = line.strip().split('\t')
            line = [x if not x == '\\N' else None for x in line ]
            pfam.append({'accession': line[0], 'clan': line[1], 'clan_name': line[2], 'name': line[4]})
    
    pfam_dict = {x['accession']:x for x in pfam}
    return pfam_dict

def domain_str_to_descr(domain_str):
    # 'PF02861-PF02861-PF00004-PF07724' --> 'Clp amino terminal domain; Clp amino terminal domain; ATPase family associated with various cellular activities (AAA); AAA domain (Cdc48 subfamily)'
    domain_str = re.split('-|~', domain_str)
    return '; '.join([pfam_dict[x]['name'] for x in domain_str])


def make_fasta_from_group(group):
    """
    From a domain architecture group, generate a fasta file for that group so it can be clustered with an external program
    group is a dict with keys: 'Pfam', '_id', 'pID'    
    """
    seq_dict = {'>' + str(x['_id']) + '||' + x['d']:x['s'] for x in protDB.find({'_id':{'$in':group['pID']}})}
    fasta_out = str(group['_id']) + ".fasta"
    with open(fasta_out,"w") as f:
        for k,v in seq_dict.items():
            print(k, file = f)
            print(v, file = f)
            
def mmseq_group(fasta_out):
    """
    cluster a fasta file using mmseq
    """
    group_id = fasta_out.split('.')[0]
    fasta_db = fasta_out.replace(".fasta", ".db")
    out_db = fasta_db.replace(".db", ".out")
    out_fasta = out_db + '.fasta'
    tmp_dir = "tmp" + group_id
    
    #setup mmseq
    os.environ["MMDIR"] = "/home/gstupp/bin/mmseqs/"
    os.environ["LD_LIBRARY_PATH"] = os.environ["MMDIR"] + "lib/ffindex/src"
    mmseq_path = os.environ["MMDIR"] + "bin/"
    
    os.system(mmseq_path + "fasta2ffindex " + fasta_out + " " + fasta_db)
    os.system("mkdir " + tmp_dir)
    #os.system(mmseq_path + "mmseqs_cluster " + fasta_db + " " + out_db + " " + tmp_dir)
    command = [mmseq_path + "mmseqs_cluster", fasta_db, out_db, tmp_dir]
    p = subprocess.Popen(command, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        print("something happened")
        return None
    
    os.system(mmseq_path + "ffindex2fasta " + out_db + " " + out_fasta)
    # puts weird shit at the end
    # $ tr -cd $'[[:alnum:]]._-|\\n>' < 9529635.out.fasta > 9529635.out.strip.fasta


def parser(fasta_file_handle):
    defline, members = '', []
    for line in fasta_file_handle:
        line = line.decode('utf-8')
        if line[0] == '>':
            if defline:
                yield {'id': defline, 'members': members}
            defline, members = line[1:].split('||')[0], []
        else:
            members.append(line.split('||')[0])
    if defline:
        yield {'id': defline, 'members': members}
  
def parse_mmseq(out_fasta):
    """
    Parse resulting cluster file from mmseq
    """
    f = open(out_fasta, 'rb')
    next(parser(f))
    
    

pfam_dict = parse_pfam_clans()


