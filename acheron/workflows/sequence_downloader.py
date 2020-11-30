import pandas as pd

def find_seq_sources(dfs):
    """
    Takes DataFrame list of DataFrames containing biosamples

    Returns dictionary {database1:[SAMN01, SAMN02],
                        database2:[SAMN02, SAMN03],
                        database3:[SAMN04]}
    """
    sources = {}

    for df_name in dfs.keys():
        if df_name.upper() == 'NCBI':
            sources[df_name] = dfs[df_name]['BioSample']
        elif df_name.upper() == 'PATRIC':
            #sources[df_name] = dfs[df_name]['biosample_accession']
            sources[df_name] = dfs[df_name]['BioSample']
        else:
            raise Exception("key not set for database {}".format(df_name))

    return sources

def download_sources(seq_sources):
    """
    Takes list of labels (including BioSamples) and list of allowed databases
    Returns a dictionary of which database to download each sample from
    Which looks like:
        {database1:[SAMN01, SAMN02],
         database2:[SAMN03],
         database3:[SAMN04]}
    """
    dl_lists = {}
    downloaded = []

    # As we assign samples to a database, save it to 'downloaded' so we
    # dont end up downloading the same sequence twice.
    for database in seq_sources.keys():
        for sample in seq_sources[database]:
            if sample in downloaded:
                # dont include, already attributed to a database
                pass
            else:
                # record we are applying it to this database
                downloaded.append(sample)

                # assign sample to database
                if database not in dl_lists.keys():
                    dl_lists[database] = []
                dl_lists[database].append(sample)

    return dl_lists
