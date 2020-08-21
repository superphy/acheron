import os
import time
import sys
import pandas as pd


from acheron.workflows.NCBI_antibiogram_downloader import query_to_df
from acheron.workflows.PATRIC_antibiogram_downloader import download_PATRIC
#from acheron.workflows.antibiogram_tools import *

def merge_antibiogram(df1, df2):
    
    merged_df = df1.merge(df2, on=['BioSample'], how='outer')


    # TODO: filter matching rows i.e. MIC_X MIC_Y to check for conflicts before merging

    return merged_df


def download_antibiogram(database, pathogen, email, antimicrobial, path, use_local, check_date):
    antimicrobials = ['amoxicillin/clavulanic acid', 'ampicillin', 'azithromycin',
    'cefoxitin', 'ceftiofur', 'ceftriaxone', 'chloramphenicol', 'ciprofloxacin',
    'gentamicin', 'nalidixic acid', 'streptomycin', 'sulfisoxazole', 'tetracycline',
    'trimethoprim/sulfamethoxazole','kanamycin']

    mics = ['AMC','AMP','AZM','FOX','TIO','CRO','CHL','CIP','GEN','NAL','STR','FIS',
    'TET','SXT','KAN']

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
                ncbi_df = pd.read_csv("data/NCBI_{}_antibiogram.csv".format(pathogen))
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
                patric_df = pd.read_csv("data/PATRIC_{}_antibiogram.csv".format(pathogen))
            except:
                print('PATRIC AMR data not found, please download before passing `use_local`')
                raise
        else:
            patric_df = download_PATRIC(pathogen, antimicrobial)
            patric_df.to_csv("data/PATRIC_{}_antibiogram.csv".format(pathogen))
        patric_df.rename(columns = {'biosample_accession':'BioSample'}, inplace = True)
        mergeable_dfs.append(patric_df)

    # take the first df only
    df = mergeable_dfs[0]

    # if there is more than one df, merge them in 1 by 1
    if len(mergeable_dfs) != 1:
        for extra_abx_df in mergeable_dfs[1:]:
            df = merge_antibiogram(df, extra_abx_df)

    df.to_csv(path)


def download_genomes(input, output):
    print('genome download not yet setup')
    print("Downloading genomes missing from {} but found in {}".format(output, input))
