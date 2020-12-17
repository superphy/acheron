import joblib
import sys
import os

from collections import Counter
from Bio import Seq, SeqIO

import pandas as pd
import numpy as np

def make_master_fasta(dataset):
    """
    Takes a directory of fastas and merges them into a single fasta
    which can then be turned into a blast database
    """
    genomes_path = "data/{}/wgs/clean".format(dataset)
    save_path = "data/{}/wgs/master.fasta".format(dataset)
    # remove master if it already exists, so we dont have duplicate entries
    if os.path.isfile(save_path):
        os.remove(save_path)

    fastas = [files for r,d,files in os.walk(path)][0]
    fastas = [genomes_path+'/'+i for i in fastas]

    print("merging {} fastas".format(len(fastas)))

    for fasta in fastas:
        for record in SeqIO.parse(fasta):
            # header is record.id, sequence is record.seq
            seq_name = fasta.split("/")[-1].split('.')[0]
            new_header = ">{}_{}".format(seq_name, record.id)

            with open(save_path, 'a') as fh:
                fh.write(new_header)
                fh.write("\n")
                fh.write(record.seq)
                fh.write("\n")


def find_top_feats(model, model_type, n):
    """
    Returns a list of tuples of the top n features from a model and their
    relative importance in the form [(feature, importance)]
    """

    # For XGBoost objects
    if model_type.upper() in ['XGB', 'XGBOOST']:
        feats = Counter(model.get_score(importance_type='gain')).most_common(n)
        feats.sort(key=lambda i: i[1], reverse=True)

        return feats

    else:
        raise Exception("Model type {} not supported for feature extraction".format(model_type))

# prep important k-mers for blast
def save_query(importances, save_path):
    seqs = [i[0] for i in importances]

    # remove query if it already exists
    if os.path.isfile(save_path):
        os.remove(save_path)

    # save seqs to fasta format
    with open(save_path, 'a') as fh:
        for seq in seqs:
            fh.write(">{}\n".format(seq))
            fh.write(seq+"\n")

def biosamples_to_genomeids(biosamples, pathogen):
    """
    Takes a list of biosamples, returns PATRIC genome id's
    """
    #patric_metadata = pd.read_csv("data/PATRIC_genome_metadata.tsv",delimiter='\t', dtype=str)
    patric_metadata = pd.read_csv("data/PATRIC_{}_antibiogram.csv".format(pathogen), dtype=str)
    patric_metadata = patric_metadata[[
        'genome_id','biosample_accession','assembly_method']]

    ids = []

    for sample in biosamples:
        row = patric_metadata[patric_metadata['biosample_accession']==sample]

        # if there is only one id for this biosample
        if len(row) == 1:
            id = row['genome_id'].values[0]

        # if there are ids for both the raw reads and the assembly
        # we want to return the id for the assembly
        else:
            row = row[row['assembly_method'] == 'CLC Genomics Workbench v. 6.0.2']
            try:
                id = row['genome_id'].values[0]
            except:
                print("No row for {} found in PATRIC AMR csv".format(sample))
                raise

        ids.append(str(id))

    return ids
