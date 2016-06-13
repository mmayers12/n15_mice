""" 
Bunch of utility functions

"""

import dill
import gzip
import os
import shelve


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
