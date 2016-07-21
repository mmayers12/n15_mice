import os
import pandas as pd
import numpy as np

BASE = os.getcwd()

def strip_date(name):
    return '_'.join(name.split('-')[1:])

"""
Names
20151112-UnerFecal-Lys-N15-1016     20160122-UnerFecal-Lys-N14N15-1121  20160223-BioGlyAOMK-Lys-IL10-N14
20151113-UnerFecal-Lys-N15-1019     20160124-AOMKFecal-Lys-N14N15-1120  20160312-BioGlyCMK-Lys-UC13
20151114-UnerFecal-Lys-N15-1021     20160125-CMKFecal-Lys-N14N15-1120   20160319-DMSO-Lys-UC13
20151207-UnerFecal-Lys-N14N15-1021  20160126-UnerFecal-Lys-N14N15-1121  20160320-DMSO-Lys-UC13
20151208-UnerFecal-Lys-N14N15-1019  20160127-UnerFecal-Lys-N14N15-1121  20160331-UnerFecal-Lys-N14N15-Pool
20151209-UnerFecal-Lys-N14N15-1016  20160129-UnerFecal-Lys-N14N15-1021  20160402-UnerFecal-Lys-N14-Pool
20151210-UnerFecal-Sup-N14N15-1021  20160216-UnerFecal-Lys-N14N15-1111  20160403-UnerFecal-Lys-N15-Pool
20151211-UnerFecal-Sup-N14N15-1019  20160218-UnerFecal-Lys-N15-1111     20160409-BioGlyCMK-Lys-N14N15-Pool
20151212-UnerFecal-Sup-N14N15-1016  20160219-UnerFecal-Lys-N14N15-1111  20160410-DMSO-Lys-N14N15-Pool
20151213-UnerFecal-Lys-N14N15-1121  20160220-UnerFecal-Lys-N14N15-1111  20160411-AOMK-Lys-N14N15-Pool
"""
def make_name(name, N14 = True):
    out_name = ''
    #Enrichment type
    if 'Uner' in name:
        out_name += 'U'
    elif 'AOMK' in name:
        out_name += 'A'
    elif 'CMK' in name:
        out_name += 'C'
    elif 'DMSO' in name:
        out_name += 'D'

    #Sample_type
    if 'Lys' in name:
        out_name += 'L'
    elif 'Sup' in name:
        out_name += 'S'

    #Sample Type (N14, N15 or Mix)
    if 'N14N15' in name:
        out_name += '_Mix'
    elif 'N14' in name:
        out_name += '_Tc'
    elif 'N15' in name:
        out_name += '_Rg'

    #Collection Date and Rep
    if N14:
        out_name += '_' + '_'.join(name.split('_')[-2:])
    else:
        out_name += '_' + '_'.join(name.split('_')[-3:]).replace('n15', 'N')

    return out_name

def make_line(folder, sample, N14 = True):
    line_d = dict()
    #name
    line_d['name'] = make_name(sample, N14)
    #sample_type
    if 'Lys' in sample:
        line_d['sample_type'] = 'Lysate'
    else:
        line_d['sample_type'] = 'Supernatant'
    #Enriched & Type
    if "Uner" in sample:
        line_d['enriched'] = False
        line_d['probe'] = np.nan
    else:
        line_d['enriched'] = True
        if 'AOMK' in sample:
            line_d['probe'] = 'AOMK'
        if 'CMK' in sample:
            line_d['probe'] = 'CMK'
        if 'DMSO' in sample:
            line_d['probe'] = 'DMSO'
    #Technical Rep #
    line_d['technical'] = sample.split('_')[-1]
    #Collection Date
    tmp = folder.split('-')[-1]
    if tmp == 'Pool':
        line_d['col_date'] = '2015-11-18 2015-11-22'
    else:
        line_d['col_date'] = '2015-' + tmp[:2] + '-' + tmp[2:]
    #run date date
    date = folder.split('-')[0]
    line_d['run_date'] = date[:4] + '-' + date[4:6] + '-' + date[6:8]
    # N14 and path
    if N14:
        line_d['n15'] = False
        line_d['path'] = os.path.join(BASE, folder, "rawXtractor/dta-pfp-0.01/DTASelect-filter.txt")
    else:
        line_d['n15'] = True
        line_d['path'] = os.path.join(BASE, folder, "rawXtractor/n15_search/dta-pfp-0.01/DTASelect-filter.txt")
    # Determine if its an N14 N15 mix file for quant, and if so store paired DTA output
    if 'N14N15' in folder:
        line_d['quant'] = True
        line_d['l_dta'] = os.path.join(BASE, folder, "rawXtractor/dta-pfp-0.01/DTASelect-filter.txt")
        line_d['h_dta'] = os.path.join(BASE, folder, "rawXtractor/n15_search/dta-pfp-0.01/DTASelect-filter.txt")
        line_d['comb_dta'] = os.path.join(BASE, folder, "rawXtractor/combined_search/DTASelect-filter.txt")
        line_d['census'] = os.path.join(BASE, folder, "rawXtractor/combined_search/census-out.txt")
    else:
        line_d['quant'] = False
        line_d['l_dta'] = np.nan
        line_d['h_dta'] = np.nan
        line_d['comb_dta'] = np.nan
        line_d['census'] = np.nan
    return line_d

def get_category(row):
    """
    Categories:  Rag_Unenriched, RT_Unenriched, RAG_BioGlyCMK, RT_BioGlyCMK.
    Others can be added later, but for the use of this study, these are the important ones.

    Primary use of category will be for grouping in PCA plots.
    """
    if row['n15'] and not row['enriched']:
        return 'RAG Unenriched'
    elif not row['n15'] and not row['enriched']:
        return 'RT Unenriched'
    elif row['n15'] and row['probe'] == 'CMK':
        return 'RAG BioGlyCMK'
    elif not row['n15'] and row['probe'] == 'CMK':
        return 'RT BioGlyCMK'
    else:
        return np.nan

def main():
    folders = sorted([x for x in os.listdir() if '201' in x])
    samples = list(map(strip_date, folders))

    for sample in set(samples):
        i = 1
        while(sample in samples):
            idx = samples.index(sample)
            samples[idx] = samples[idx]+'_'+str(i)
            i += 1
    metadata = []

    for folder, sample in zip(folders, samples):
        if 'N14' in sample:
            metadata.append(make_line(folder,sample))
        if 'N15' in sample:
            ssplt = sample.split('_')
            sample = '_'.join(ssplt[:-1]) + '_n15_' + ssplt[-1]
            metadata.append(make_line(folder,sample, False))

    # put the metadata into a csv
    metadata = pd.DataFrame(metadata)
    metadata = metadata.set_index('name')
    metadata['category'] = metadata.apply(get_category, axis=1)
    metadata.to_csv('metadata.csv', index_label = '')


if __name__ == '__main__':
    main()
