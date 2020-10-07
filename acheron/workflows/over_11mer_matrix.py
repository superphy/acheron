import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import os, sys
from itertools import repeat

def make_row(filename, dataset, kmer_length):
    from Bio import Seq, SeqIO
    """
    Given a genome file, create and return a row of kmer counts
    to be inserted into the kmer matrix.
    """
    relevant_feats = np.load("data/{}/wgs/master_{}mers/all_kmers.npy".format(
    dataset, kmer_length))

    cols_dict = { relevant_feats[i] : i for i in range(0, len(relevant_feats))}

    # Create a temp row to fill and return (later placed in the kmer_matrix)
    temp_row = [0]*len(relevant_feats)

    jf_path = "data/{}/wgs/{}mer_jellyfish_results/{}".format(
    dataset, kmer_length, filename)

    for record in SeqIO.parse(jf_path[:-3], "fasta"):
        # Retrieve the sequence as a string
        kmer_seq = record.seq
        kmer_seq = str(kmer_seq)
        #kmer_seq = kmer_seq._get_seq_str_and_check_alphabet(kmer_seq)

        kmer_count = int(record.id)
        if kmer_count>255:
            kmer_count = 255
        temp_row[cols_dict[kmer_seq]] = kmer_count

    return filename, temp_row


def make_matrix(dataset, kmer_length, matrix_dtype, num_threads, split_num):

    relevant_feats = np.load("data/{}/wgs/master_{}mers/all_kmers.npy".format(
    dataset, kmer_length))

    splits_array = np.load("data/{}/wgs/master_{}mers/matrix_splits.npy".format(
    dataset, kmer_length), allow_pickle=True)

    genomes = splits_array[int(split_num)-1]
    genomes = [i.split('/')[-1] for i in genomes]

    total = len(genomes)

    runs = [i.split('.')[0] for i in genomes]

    # declaring empty kmer matrix to fill
    kmer_matrix = np.zeros((len(genomes),len(relevant_feats)),dtype = matrix_dtype)

    # making dicts for faster indexing
    # note that rows dict is in filenames not genome/run names
    rows_dict = { genomes[i] : i for i in range(0, len(genomes))}
    cols_dict = { relevant_feats[i] : i for i in range(0, len(relevant_feats))}

    # Use concurrent futures to get multiple rows at the same time
    # Then place completed rows into the matrix and update the row dictionary
    with ProcessPoolExecutor(max_workers=min(16, num_threads)) as ppe:
        for genome_name,temp_row in ppe.map(make_row, genomes, repeat(dataset),
        repeat(kmer_length)):
            for i, val in enumerate(temp_row):
                kmer_matrix[rows_dict[genome_name]][i] = val

    df = pd.DataFrame(data = kmer_matrix, columns = relevant_feats, index = runs)
    return df
