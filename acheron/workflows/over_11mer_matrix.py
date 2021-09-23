import numpy as np
import pandas as pd
import gc
from concurrent.futures import ProcessPoolExecutor
import os, sys
from itertools import repeat
from sklearn.feature_selection import f_classif

def make_row(filename, dataset, kmer_length):
    from Bio import Seq, SeqIO
    """
    Given a genome file, create and return a row of kmer counts
    to be inserted into the kover_11mer_matrix.mer matrix.
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
    genomes = [i.split('/')[-1] for i in genomes if i != None]

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

def get_feat_splits(all_feats, RAM_denom):
    """
    Takes a list of features, splits it into a dictionary
    {0:[feat1,feat2,feat3], 1:[feat4,feat5......
    """
    assignments = []
    num_feats = len(all_feats)

    feats_per = num_feats//RAM_denom

    if num_feats % RAM_denom != 0:
        feats_per += 1

    for i in range(RAM_denom):
        assignments.append(all_feats[i*feats_per:(i+1)*feats_per])

    return assignments

def large_feat_select(file_paths, RAM_denom, num_feats, labels):
    """
    Takes a list of filepaths to slices of a feature matrix
    Does feature selection on a subset of data at a time as to not blow through the ram.

    If RAM_denom is set to 5, only 20% of each slice will be kept in ram.
    You are passing the denominator of a fraction, i.e. 1/2, 1/3, 1/4, 1/20
    i.e. if you pass 10, then it will reduce ram per slice to 1/10th

    Ram use will then be 0.2*(memory of 1 slice)*(number of slices loaded) + (memory of 1 slice)
    i.e., it will use more than 20% of the RAM as before, one part of the problem is reduced by 80%
    """

    f_values = []

    #  First load one slice, prime everything, then we will loop through the rest after
    slice1 = pd.read_pickle(file_paths[0])
    all_feats = slice1.columns
    del slice1

    feat_splits = np.array_split(all_feats, RAM_denom)

    for slice_num in range(RAM_denom):
        # as stated above, now we loop through the rest
        slices = []
        gc.collect()

        for i, file in enumerate(file_paths):
            if i == 0:
                slices = pd.read_pickle(file)[feat_splits[slice_num]]
            else:
                new_slice = pd.read_pickle(file)[feat_splits[slice_num]]
                slices = pd.concat([slices,new_slice])
                del new_slice
                gc.collect()

        f_vals, p_vals = f_classif(slices.values,labels)
        for i in range(len(f_vals)):
            f_values.append([f_vals[i],feat_splits[slice_num][i]])

    feats_to_keep = sorted(f_values, reverse=True)[:num_feats]
    feats_to_keep = [i[1] for i in feats_to_keep]

    slices = []
    for i, file in enumerate(file_paths):
        if i == 0:
            slices = pd.read_pickle(file)[feats_to_keep]
        else:
            new_slice = pd.read_pickle(file)[feats_to_keep]
            slices = pd.concat([slices,new_slice])
            del new_slice
            gc.collect()

    return slices


    # slices look like This
    """
                          AGCCTGGCTGGTCGGCTGTACAAAGACGAAT  ...  GTAGCATGAATGGGGGTAATCTGGAATGGAA
    SAMN02640777                                0  ...                                0
    SAMN02640827                                0  ...                                0
    SAMN02640878                                0  ...                                0
    SAMN02640928                                0  ...                                0
    SAMN02699192                                0  ...                                0
    ...                                       ...  ...                              ...
    SAMN05596337                                0  ...                                0
    SAMN05596788                                0  ...                                0
    SAMN05771752                                0  ...                                0
    SAMN05771803                                0  ...                                0
    SAMN05907766                                0  ...                                0

    [127 rows x 72504712 columns]
    """
