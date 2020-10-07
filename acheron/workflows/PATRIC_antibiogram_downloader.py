import wget
import os.path
import sys
import pandas as pd
from collections import Counter

def pull_PATRIC_ftp():
    """
    Gets new AMR and metadata information from the PATRIC ftp server
    """
    # do not write test for, file size is 50MB + 200MB
    PATRIC_url = "ftp://ftp.patricbrc.org/RELEASE_NOTES/PATRIC_genomes_AMR.txt"
    wget.download(PATRIC_url, out="data/PATRIC_genomes_AMR.txt")

    if not os.path.isfile("data/PATRIC_genome_metadata.tsv"):
        metadata_url = "ftp://ftp.patricbrc.org/RELEASE_NOTES/genome_metadata"
        wget.download(metadata_url, out="data/PATRIC_genome_metadata.tsv")

def filter_df(df, pathogen, antimicrobials):
    """
    Removes excess pathogens and antimicrobials
    """
    # use only MIC (mg/L), not disk diffusion (mm) values
    df = df[df['measurement_unit']=='mg/L']

    # filter for requested pathogen
    df = df[df['genome_name'].str.contains(pathogen, na=False)]

    for abx in antimicrobials:
        if abx not in df['antibiotic'].values:
            print("{} requested but not found in PATRIC datasheet".format(abx))

    abx_counts = Counter(df['antibiotic'].values)
    unused_abx = [i for i in abx_counts.keys() if i not in antimicrobials]
    print("\nSkipping the following antimicrobials because they were not requested (and how many times they are seen):")

    print(["{} ({})".format(i, abx_counts[i]) for i in unused_abx])

    # filter for requested antimicrobials
    df = df[df['antibiotic'].isin(antimicrobials)]

    return(df)

def columnize_tsv(df):
    """
    Currently in the following format:
    Row 1: pathogen, drug1
    Row 2: pathogen, drug2

    This function converts to:
    Row 1: pathogen, drug1, drug2
    """

    long_names = ['amoxicillin/clavulanic acid', 'ampicillin', 'azithromycin',
    'cefoxitin', 'ceftiofur', 'ceftriaxone', 'chloramphenicol', 'ciprofloxacin',
    'gentamicin', 'nalidixic acid', 'streptomycin', 'sulfisoxazole', 'tetracycline',
    'trimethoprim/sulfamethoxazole','kanamycin']

    three_letter_codes = ['AMC','AMP','AZM','FOX','TIO','CRO','CHL','CIP','GEN','NAL','STR','FIS',
    'TET','SXT','KAN']

    abx = set(df['antibiotic'])
    ids = set(df['genome_id'])

    cols = [i for i in df.columns if i!='antibiotic'] + ['MIC_'+three_letter_codes[long_names.index(i)] for i in abx]
    sorted_df = pd.DataFrame(columns=cols)

    duplicate_mics = 0
    conflicting_mics = 0
    conflicts = {}

    for i in ids:
        single_id_df = df[df["genome_id"]==i]

        to_add_df = pd.DataFrame(columns=cols)
        for seen_col in single_id_df.columns:
            if len(set(single_id_df[seen_col])) == 1:
                to_add_df[seen_col] = pd.Series([single_id_df[seen_col].iloc[0]])
            else:
                to_add_df[seen_col] = pd.Series(['mixed'])

        for available_abx in single_id_df['antibiotic'].values:
            short_code = 'MIC_'+three_letter_codes[long_names.index(available_abx)]

            # find rows matching available_abx
            matching_abx_row = single_id_df[single_id_df['antibiotic']==available_abx]

            # check if there are 2 values for this genome+abx combo, and are they just a rounding error
            if len(matching_abx_row['measurement']) > 1:
                in_conflict=False
                all_mics = matching_abx_row['measurement']
                for equality_smb in ['=','<','>']:
                    all_mics = [i.replace(equality_smb,'') for i in all_mics]

                fewest_digits = min(all_mics, key=len)

                # check all values to see if the first n digits match the shortest
                # so for example if 0.125 and 0.12 are both seen, they would be considered equal
                for extra_abx in all_mics:
                    if extra_abx[:len(fewest_digits)] != fewest_digits:
                        #print("For genome {}, there are multiple non equal MIC values for {}, this abx will be ignored".format(i,available_abx))
                        conflicting_mics +=1
                        in_conflict=True
                        if i in conflicts.keys():
                            conflicts[i] = conflicts[i] + [available_abx]
                        else:
                            conflicts[i] = [available_abx]
                    else:
                        duplicate_mics +=1
                if in_conflict:
                    continue
            # add the first (and most likely only) value to the main df
            to_add_df[short_code] = pd.Series([matching_abx_row['measurement'].iloc[0]])

        sorted_df = sorted_df.append(to_add_df)

    sorted_df = sorted_df.drop(columns=[
        "antibiotic","resistant_phenotype","measurement",
        "measurement_sign","measurement_value","measurement_unit"])

    print("\n{} duplicate but matching mic values were found".format(duplicate_mics))
    print("{} MIC values were ignored because there were multiple conflicting values, these were:".format(conflicting_mics))
    for k in conflicts.keys():
        print("{}: {}".format(k,conflicts[k]))

    return sorted_df



def add_metadata(amr_df, metadata_df):
    """
    Adds the extra information in the metadata tsv into the amr df
    """
    merged_df = amr_df.merge(metadata_df, on=['genome_id','genome_name','taxon_id'], how='left')
    return merged_df

def download_PATRIC(pathogen, antimicrobials):
    pull_PATRIC_ftp()

    # load the data, filter out unwanted pathogens and antimicrobials
    amr_df = pd.read_csv("data/PATRIC_genomes_AMR.txt", delimiter='\t',low_memory=False, dtype=str)
    df = filter_df(amr_df, pathogen, antimicrobials)

    # put all matching samples in a single row
    df = columnize_tsv(df)

    # adds metadata including biosample to be used as id
    metadata_df = pd.read_csv("data/PATRIC_genome_metadata.tsv",delimiter = '\t', dtype=str)
    df = add_metadata(df, metadata_df)

    return df
