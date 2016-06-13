import os



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

def make_line(line, sample, N14 = True):
    temp_line = ""
    #name
    temp_line += make_name(sample, N14) + ','
    #sample_type
    if 'Lys' in sample:
        temp_line += 'Lysate,'
    else:
        temp_line += 'Supernatant,'
    #Enriched & Type
    if "Uner" in sample:
        temp_line += 'FALSE,N/A,'
    else:
        temp_line += 'TRUE,'
        if 'AOMK' in sample:
            temp_line += 'AOMK,'
        if 'CMK' in sample:
            temp_line += 'CMK,'
        if 'DMSO' in sample:
            temp_line += 'DMSO,'
    #Technical Rep #
    temp_line += sample.split('_')[-1] + ','
    #Collection Date
    tmp = line.split('-')[-1]
    if tmp == 'Pool':
        temp_line += '2015-11-18 2015-11-22,'
    else:
        temp_line += '2015-' + tmp[:2] + '-' + tmp[2:] + ','
    #run date date
    date = line.split('-')[0]
    temp_line += date[:4] + '-' + date[4:6] + '-' + date[6:8] + ','
    # N14 and path
    if N14:
        temp_line += 'FALSE,'
        temp_line += (line + "/rawXtractor/dta-pfp-0.01/DTASelect-filter.txt")
    else:
        temp_line += 'TRUE,'
        temp_line += (line + "/rawXtractor/n15_search/dta-pfp-0.01/DTASelect-filter.txt")

    return temp_line


lines = sorted([x for x in os.listdir() if '201' in x])

samples = list(map(strip_date, lines))

for sample in set(samples):
    i = 1
    while(sample in samples):
        idx = samples.index(sample)
        samples[idx] = samples[idx]+'_'+str(i)
        i += 1

new_lines = []
for i, line in enumerate(lines):
    if 'N14' in line:
        new_lines.append(make_line(line,samples[i]))
    if 'N15' in line:
        ssplt = samples[i].split('_')
        samples[i] = '_'.join(ssplt[:-1]) + '_n15_' + ssplt[-1]
        new_lines.append(make_line(line,samples[i], False))


header = ',sample_type,enriched,probe,technical,col_date,run_date,n15,path'



with open("metadata.csv", 'w') as f:
    f.write(header + '\n')
    f.write("\n".join(new_lines))
