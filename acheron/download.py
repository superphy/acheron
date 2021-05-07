import os
import time
import sys
import pandas as pd
import numpy as np
import re


from acheron.workflows.NCBI_antibiogram_downloader import query_to_df
from acheron.workflows.PATRIC_antibiogram_downloader import download_PATRIC

def are_equal_mic(mic1, mic2):
    """
    Determines if 2 MIC values are functionally equivalent,
    such as >32 and >= 32, or 0.12 and 0.125
    """
    # MIC values within 10% considered equal (likely rounding error)
    allowable_diff = 0.1

    mic1 = re.sub('[<>= ]', '', mic1)
    mic2 = re.sub('[<>= ]', '', mic2)

    if isinstance(mic1,str) and mic1[0]=='.':
        mic1 = '0'+mic1
    if isinstance(mic2,str) and mic2[0]=='.':
        mic2 = '0'+mic2

    try:
        mic_float_1 = float(mic1)
    except:
        print("\nValue: '{}' found in dataframe, unable to cast to MIC value\n".format(mic1))
        return ['uncastable', 'lhs']
    try:
        mic_float_2 = float(mic2)
    except:
        print("\nValue: '{}' found in dataframe, unable to cast to MIC value\n".format(mic2))
        return ['uncastable', 'rhs']

    ratio = mic_float_1/mic_float_2
    diff = abs(ratio-1)

    if diff > allowable_diff:
        # also forgive rounding errors:
        shortest_length = len(min([mic1,mic2], key=len))

        if mic1[:shortest_length] == mic2[:shortest_length]:
            return True
        else:
            return False
    else:
        return True

def is_empty(mic):
    """
    Returns true if MIC is nan or empty, not to be confused with invalid
    """
    if mic in ['','NaN']:
        return True
    try:
        if mic.isspace():
            return True
    except:
        # float nans fail this check
        pass
    try:
        if np.isnan(mic):
            return True
    except:
        # strings might fail this check
        pass

    return False

def merge_antibiogram(df1, df2):

    merged_df = df1.merge(df2, on=['BioSample'], how='outer')

    all_cols = merged_df.columns

    # cols that appear in both dataframes withh end in _X or _Y, we just need one
    dup_cols = [i for i in all_cols if i[-2:] in ['_x', '_x']]

    conflicts = {}
    fixed = []

    # for every duplicated columns, we need to merge them and set conflicts to NaN
    for dup in [i[:-2] for i in dup_cols if i[-2:] == '_x']:
        # add the new column
        #merged_df[dup] = pd.Series([i if i==j else np.nan for i,j in zip(merged_df[dup+'_x'],merged_df[dup+'_y'])])

        if 'MIC' in dup:
            is_mic = True
        else:
            is_mic = False

        new_col = []

        for itr, lhs in enumerate(merged_df[dup+'_x']):
            rhs = merged_df[dup+'_y'][itr]

            # if left and right hand sides are the same, we accept it as correct, or if only one side has a valid answer
            # if the values are not the same, set to 'invalid' NOT NaN!!, record error
            # Also keep in mind that NCBI data is formated with a space i.e. <= 32, whereas PATRIC is <=32

            if is_empty(lhs):
                new_col.append(rhs)
                continue
            elif is_empty(rhs):
                new_col.append(lhs)
                continue

            if is_mic:
                are_equal = are_equal_mic(lhs, rhs)
                if not isinstance(are_equal, bool):
                    if lhs == rhs:
                        raise Exception("Both dataframes contain uncastable value: {} under {}".format(lhs,dup))
                    elif are_equal[1] =='lhs':
                        new_col.append(rhs)
                    else:
                        new_col.append(lhs)
                elif are_equal:
                    new_col.append(lhs)
                else:
                    # in this case we have an unresolvedable conflict

                    # record conflict
                    if dup in conflicts.keys():
                        conflicts[dup] += 1
                    else:
                        conflicts[dup] = 1

                    # this HAS to be set to invalid, if its set to NaN,
                    # another dataframe having any value will override the NaN
                    new_col.append('invalid')
            else:
                if lhs == rhs:
                    new_col.append(lhs)
                else:
                    new_col.append('invalid')

        merged_df[dup] = pd.Series(new_col)
        merged_df = merged_df.drop(columns=[dup+'_x',dup+'_y'])
        fixed.append(dup)

    print("Automatically merging the following columns:", fixed)

    # put the important stuff first
    mic_cols = [i for i in merged_df.columns if 'MIC' in i]
    sir_cols = [i for i in merged_df.columns if 'SIR' in i]
    amr_cols = ['BioSample'] + mic_cols + sir_cols

    # other columns after
    mic_first_order = amr_cols + [i for i in merged_df.columns if i not in amr_cols]

    if len(conflicts.keys()) > 0:
        print('\nThe following columns had conflicting data between columns (and the number of removed values):')
        print(conflicts)

    return merged_df.reindex(columns = mic_first_order)

def download_antibiogram(database, pathogen, email, antimicrobial, path, use_local, check_date):
    antimicrobials = ['amoxicillin/clavulanic acid', 'ampicillin', 'azithromycin',
    'cefoxitin', 'ceftiofur', 'ceftriaxone', 'chloramphenicol', 'ciprofloxacin',
    'gentamicin', 'nalidixic acid', 'streptomycin', 'sulfisoxazole', 'tetracycline',
    'trimethoprim/sulfamethoxazole','kanamycin', 'clindamycin', 'erythromycin',
    'florfenicol', 'telithromycin']

    mics = ['AMC','AMP','AZM','FOX','TIO','CRO','CHL','CIP','GEN','NAL','STR','FIS',
    'TET','SXT','KAN', 'CLI','ERY', 'FLO', 'TEL']

    if(check_date):
        for db_check, db_path in [["PATRIC","data/PATRIC_genomes_AMR.txt"],["NCBI","data/NCBI_{}_antibiogram.csv".format(pathogen)]]:
            try:
                print("{} AMR data last pulled on".format(db_check),
                    time.ctime(os.path.getctime(db_path)))
            except:
                print("{} AMR data not found".format(db_check))
        return


    print("Looking in {} database(s):".format(len(database)),database)
    print("for {} antibiogram data using the email {}\n".format(pathogen, email))
    if(antimicrobial!= 'all'):
        raise Exception("not yet setup for individual abx")
    else:
        antimicrobial = antimicrobials

    if use_local is None:
        use_local = ['']

    mergeable_dfs = []

    # NCBI
    if 'NCBI' in database:
        if 'NCBI' in use_local:
            try:
                ncbi_df = pd.read_csv("data/NCBI_{}_antibiogram.csv".format(pathogen), dtype=str)
            except:
                print('NCBI AMR data not found, please download before passing `use_local`')
                raise
        else:
            query = "antibiogram[filter] AND {}[organism]".format(pathogen)
            from Bio import Entrez
            Entrez.email = email
            ncbi_df = query_to_df(query, mics, antimicrobials)
            ncbi_df.to_csv("data/NCBI_{}_antibiogram.csv".format(pathogen))
        mergeable_dfs.append(ncbi_df)

    # PATRIC
    if 'PATRIC' in database:
        if 'PATRIC' in use_local:
            try:
                patric_df = pd.read_csv("data/PATRIC_{}_antibiogram.csv".format(pathogen), dtype=str)
            except:
                print('PATRIC AMR data not found, please download before passing `use_local`')
                raise
        else:
            patric_df = download_PATRIC(pathogen, antimicrobial)
            patric_df.to_csv("data/PATRIC_{}_antibiogram.csv".format(pathogen))
        patric_df.rename(columns = {'biosample_accession':'BioSample'}, inplace = True)
        mergeable_dfs.append(patric_df)

    # Removed any unnamed columns from residual indeces
    mergeable_dfs = [j.drop(columns=[i for i in j.columns if 'Unnamed: ' in i]) for j in mergeable_dfs]

    # take the first df only
    df = mergeable_dfs[0]

    # if there is more than one df, merge them in 1 by 1
    if len(mergeable_dfs) != 1:
        for extra_abx_df in mergeable_dfs[1:]:
            df = merge_antibiogram(df, extra_abx_df)

    df.to_csv(path)


def download_genomes(databases, output, pathogen):

    # these datasheet might need to be loaded using something other
    # than biosample numbers, depending on how PATRICS ftp server
    # operates
    for database in databases:
        if database.upper() in ['NCBI','PATRIC']:
            os.system("snakemake -k -s acheron/workflows/{}_sequence_downloader.smk -j 1 --config databases={} out={} pathogen={}".format(
                database.upper(),'_'.join(databases),output,pathogen))
        else:
            raise Exception("Database {} not yet setup for sequence download".format(database))
