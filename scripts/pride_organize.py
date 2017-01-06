"""
Script to reorganize all data and rename.
Originally created for PRIDE submission.
"""
import os
import glob
import shutil


folder_names = ['20151207-UnerFecal-Lys-N14N15-1021',  '20151208-UnerFecal-Lys-N14N15-1019',
                '20151209-UnerFecal-Lys-N14N15-1016',  '20151213-UnerFecal-Lys-N14N15-1121',
                '20160122-UnerFecal-Lys-N14N15-1121',  '20160126-UnerFecal-Lys-N14N15-1121',
                '20160127-UnerFecal-Lys-N14N15-1121',  '20160216-UnerFecal-Lys-N14N15-1111',
                '20160219-UnerFecal-Lys-N14N15-1111',  '20160220-UnerFecal-Lys-N14N15-1111',
                '20160331-UnerFecal-Lys-N14N15-Pool',  '20160409-BioGlyCMK-Lys-N14N15-Pool',
                '20160410-DMSO-Lys-N14N15-Pool',  '20160506-BioGlyCMK-Lys-N14N15-Pool',
                '20160507-DMSO-Lys-N14N15-Pool',  '20160508-BioGlyCMK-Lys-N14N15-Pool',
                '20160509-DMSO-Lys-N14N15-Pool',  '20160528-BioGlyCMK-Lys-N14N15-Pool',
                '20160529-DMSO-Lys-N14N15-Pool',  '20160530-BioGlyCMK-Lys-N14N15-Pool']


def copy_files(glob_ex, out_fold):
    for file in glob.glob(glob_ex):
        shutil.copy2(file, out_fold)

def restructure_dirs():
    """
    Create new folders for samples used in study, and copy data to folders
    """

    BASE = 'PRIDE'
    os.mkdir(BASE)

    for folder in folder_names:
        out_fold = os.path.join(BASE, folder)

        # Make new Directories
        os.makedirs(os.path.join(out_fold,'14n-search'))
        os.makedirs(os.path.join(out_fold,'15n-search'))
        os.makedirs(os.path.join(out_fold,'census-result'))

        # Copy RAW, MS1, MS2 files
        copy_files(os.path.join(folder, '*.RAW'), out_fold)
        copy_files(os.path.join(folder, 'rawXtractor', '*.ms1'), out_fold)
        copy_files(os.path.join(folder, 'rawXtractor', '*.ms2'), out_fold)

        # Move 14N files
        copy_files(os.path.join(folder,'rawXtractor', '*.sqt'), os.path.join(out_fold,'14n-search'))
        copy_files(os.path.join(folder,'rawXtractor', 'dta-pfp-0.01', 'DTASelect-filter.txt'), os.path.join(out_fold,'14n-search'))

        # Move 15N Files
        copy_files(os.path.join(folder,'rawXtractor', 'n15_search', '*.sqt'), os.path.join(out_fold,'15n-search'))
        copy_files(os.path.join(folder,'rawXtractor', 'n15_search', 'dta-pfp-0.01', 'DTASelect-filter.txt'), os.path.join(out_fold,'15n-search'))

        # Move combined census Files
        copy_files(os.path.join(folder,'rawXtractor', 'combined_search', 'census-out.txt'), os.path.join(out_fold,'census-result'))
        copy_files(os.path.join(folder,'rawXtractor', 'combined_search', 'DTASelect-filter.txt'), os.path.join(out_fold,'census-result'))
        copy_files(os.path.join(folder,'rawXtractor', 'census-out_peptide.txt'), os.path.join(out_fold,'census-result'))

def rename_files():
    """
    Change DTASelect-filter.txt, census-out.txt, and *.sqt files so all have unique names.
    Append run-date (and search type) to the front of .txt files.
    Add '-n15-' to sqt files in 15n folder so names are unique.
    """

    # Change to the PRIDE folder
    os.chdir('PRIDE')

    dirs = ['14n-search', '15n-search', 'census-result']

    for fold in folder_names:

        print fold

        # Get date for appending to results files
        date = fold.split('-')[0]

        # Change DTASelect-filter.txt result files
        for d in dirs:
            # Get result type
            f_type = d.split('-')[0]
            shutil.move(os.path.join(fold, d, 'DTASelect-filter.txt'), os.path.join(fold, d, '-'.join([date, f_type, 'DTASelect-filter.txt'])))

            # Change n15 sqt files so they're uniquely named
            if d == dirs[1]: # '15n-search' folder

                f_type = 'n15'
                for file_loc in glob.glob(os.path.join(fold, d, '*.sqt')):

                    file = os.path.basename(file_loc)

                    # UnerFecal-Pool-1118-1122a-s4.sqt
                    # add '-n15-' to name right before s##.sqt
                    new_name = '-'.join(file.split('-')[:-1] + [f_type, file.split('-')[-1]])

                    shutil.move(file_loc, os.path.join(fold, d, new_name))
            # Change the census result filenames
            elif d == dirs[-1]: # 'census-results' folder
                for file_loc in glob.glob(os.path.join(fold, d, 'census*.txt')):
                    # Get filename
                    file = os.path.basename(file_loc)
                    shutil.move(os.path.join(fold, d, file), os.path.join(fold, d, '-'.join([date,file])))

if __name__ == '__main__':
    restructure_dirs()
    rename_files()
