#! /usr/bin/env python3

"""
Generate qc plots using dtaselect-filter.txt

"""
import os
import matplotlib
if __name__ == '__main__':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import chain
from metaproteomics.file_processing import blazmass_tools
import numpy as np
from collections import Counter
import sys


def num_PSM_per_scan(dta_parser):
    num_PSM_per_scan = Counter({p["AA_Sequence"]:p['LCStep'] for p in chain(*[x['peptides'] for x in dta_parser])}.values())
    x,y = zip(*sorted(num_PSM_per_scan.items(), key = lambda x:x[0]))
    return x,y

def get_lcstep(filename):
    # Parse lc step out of filename
    # failure to parse -> returns -1
    split_on = '_'
    if '-' in filename:
        split_on = '-'
    try:
        lcstep = int(filename.split('.')[0].split(split_on)[-1])
        if lcstep > 100:
            lcstep = int(filename.split('.')[0].split(split_on)[-2])
        return lcstep
    except ValueError:
        return -1

def plot_hist(data, title, xlabel, ylabel, save=None):
    # save: folder to save plot in
    plt.figure()     
    n,bins,patches = plt.hist(data, bins=50)
    plt.title(title)
    data_median = np.median(data)
    plt.plot([data_median, data_median], [min(n),max(n)], color="r")
    plt.text(x=data_median, y=max(n)*.9, s="<- median: {0:.2f}".format(data_median))
    ax = plt.gca()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel) 
    if save:
        plt.savefig(os.path.join(save, title.replace("\n","_") + ".png"), orientation='landscape', dpi=300)
        plt.savefig(os.path.join(save, title.replace("\n","_") + ".pdf"), orientation='landscape')
        
def plot_bar(x, y, title, xlabel, ylabel, save=None):
    # save: folder to save plot in
    plt.figure()
    plt.bar(x,y)
    plt.title(title)
    ax = plt.gca()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if save:
        plt.savefig(os.path.join(save, title.replace("\n","_") + ".png"), orientation='landscape', dpi=300)
        plt.savefig(os.path.join(save, title.replace("\n","_") + ".pdf"), orientation='landscape')

def get_name_parser(dta_path):
    dta_path = os.path.abspath(os.path.expanduser(dta_path))
    name = os.path.basename(os.path.dirname(dta_path))
    dta_parser = list(blazmass_tools.dta_select_parser(dta_path))
    return name, dta_parser

def plot_ppm(dta_path, name = '', dta_parser = None, save = None):
    # method for plotting ppm, can be imported and used in notebooks
    if not name or not dta_parser:
        name, dta_parser = get_name_parser(dta_path)

    ppm_error = list(chain(*[[x['PPM'] for x in locus['peptides']] for locus in dta_parser]))
    plot_hist(ppm_error, 'ppm error\n{}'.format(name), 'ppm', 'count', save=save)

def plot_psm(dta_path, name = '', dta_parser = None, save = None):
    # method for plotting psms, can be imported and used in notebooks
    if not name or not dta_parser:    
        name, dta_parser = get_name_parser(dta_path)
        
    step,num_PSM = num_PSM_per_scan(dta_parser)
    plot_bar(step, num_PSM, "Number of PSMs per LC-Step\n{}".format(name), "LC-Step", "Number of PSMs", save=save)

if len(sys.argv) == 2:
    dta_path = sys.argv[1]
else:
    dta_path = os.getcwd()

if __name__ == "__main__":
    name, dta_parser = get_name_parser(os.path.join(dta_path, 'DTASelect-filter.txt'))  
    plot_ppm(dta_path, name, dta_parser, save = dta_path)    
    plot_psm(dta_path, name, dta_parser, save = dta_path)

