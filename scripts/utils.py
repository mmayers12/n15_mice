""" 
Bunch of utility functions

"""

import dill
import gzip
import os
import shelve
import re


def get_lcstep(filename):
    """Grab the LC step number from the step file name
    Examples of naming conventions accepted:
    20151113-unenriched-4a-s2.ms2
    20151113-unenriched-4a-s1_151113145120.ms2    
    05282015_lysed_AW_0518_Phe3_4.1364.1364.2
    05282015_lysed_AW_0518_Phe3_4_150529142251.1364.1364.2
    UnerFecal_C3_1016_S3.8590.8590.2
    UnerFecal_C3_1016_S3_151113145120.8590.8590.2
    :param filename: Name of the file for the given LC step        
    """
    try:
        lcstep = re.split('_|-', filename.split('.')[0])[-1]
        if len(lcstep) > 3:
            lcstep = re.split('_|-', filename.split('.')[0])[-2]
        if lcstep[0] == 's' or lcstep[0] == 'S':
            lcstep = lcstep[1:]
        return int(lcstep)
    except ValueError:
        return None

def get_unmod_seq(peptide_seq):	

    unmod_peptide = peptide_seq
    diff_mass = 0
    mod_list = []
    mods = dict()
    for match in re.finditer('\((.*?)\)', peptide_seq):
        mass = float(match.groups()[0])
        AA = peptide_seq[match.span()[0] - 1]
        diff_mass += mass
        pos = unmod_peptide.index(match.group())  # 1-base
        mod_list.append((AA, pos, mass))
        unmod_peptide = unmod_peptide.replace(match.group(), '', 1)

        mods['unmod_peptide'] = unmod_peptide
        mods['diff_mass'] = diff_mass
        mods['mods'] = mod_list

    return unmod_peptide, mods


def save(obj, filename, protocol=None, force=False):
    """Save a pickleable object to disk, piped through gzip
    :rtype : object
    :param obj:
    :param filename: string
    :param protocol: int
    :param force: boolean
    """
    if os.path.isfile(filename) and not force:
        raise ValueError(filename + " exists")
    with gzip.GzipFile(filename, 'wb') as f:
        dill.dump(obj, f, protocol)


def load(filename):
    """Load a pickleable object to disk, piped through gzip"""
    with gzip.open(filename, 'rb') as f:
        myNewObject = dill.load(f)
    return myNewObject


def shelve_workspace(filename, save=[], skip=[], overwrite=False, skip_prefix='_'):
    """
    Save local variables in the caller's frame
    
    save: list of strings. If given, save only these
    skip: list of strings. If given, save everything but these
    overwrite: If true, replace existing file. If false, and file exists, it will be updated
    skip_prefix: skip vars that start with this
    """

    if os.path.isfile(filename) and not overwrite:
        print("warning: " + filename + " exists")
    if os.path.isfile(filename) and overwrite:
        os.remove(filename)

    if save and skip:
        raise ValueError("Can't do both save and no_save")

    import inspect
    d = shelve.open(filename)
    frame = inspect.currentframe()
    try:
        frame_vars = frame.f_back.f_locals
        for key, value in frame_vars.items():
            if (save and key not in save) or (skip and key in skip) or key.startswith(skip_prefix):
                continue
            try:
                d[key] = value
            except:
                # PicklingError, TypeError, NotImplementedError
                print('ERROR shelving: {0}'.format(key))
    finally:
        del frame
        print(list(d.keys()))
        d.close()


def restore_workspace(filename):
    """
    Restore variables to the caller's frame
    """
    import inspect
    d = shelve.open(filename)
    frame = inspect.currentframe()
    try:
        frame_vars = frame.f_back.f_locals
        for key in d.keys():
            if key in frame_vars:
                print("Clobbering: " + key)
            try:
                frame_vars[key] = d[key]
            except:
                # PicklingError, TypeError, NotImplementedError
                print('ERROR restoring: {0}'.format(key))
    except:
        print(key)
    finally:
        del frame
        d.close()


def uniquify_by_key(d, key):
    # uniquify a list of dicts by key `key`
    seen = set()
    return [x for x in d if not (x[key] in seen or seen.add(x[key]))]
