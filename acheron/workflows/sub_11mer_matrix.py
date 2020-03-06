#!/usr/bin/env python

from Bio import Seq, SeqIO
from pathlib import Path
import numpy as np
import pandas as pd
import os
import sys
import itertools
import re
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

from clean import get_files_to_analyze


def make_row(filename, num_cols, col_names):
    """
    Given a genome file, create and return a row of kmer counts
    to be inerted into the mer matrix.
    """
    print('calling make row')
    # Filepath
    thefile = str(filename[0])

    # Get the genome id from the filepath
    genomeid = filename[0].split('/')[-1]
    genomeid = genomeid.split('.')[-2]

    # Create a temp row to fill and return (later placed in the kmer_matrix)
    temp_row = [0]*num_cols

    # Walk through the file
    for record in SeqIO.parse(thefile, "fasta"):
        # Retrieve the sequence as a string
        kmerseq = record.seq
        kmerseq = kmerseq._get_seq_str_and_check_alphabet(kmerseq)

        # Retrieve the kmer count as an int
        kmercount = record.id
        kmercount = int(kmercount)
        if kmer_count>255:
            kmer_count = 255

        # Lookup the seq in the column list for the index
        col_index = col_names[kmerseq]

        # Put the kmercount in the right spot in the row
        temp_row[col_index] = kmercount

    return genomeid,temp_row


def make_matrix(kmer_length, matrix_dtype, num_threads, results_path, save_path):
    """
    Creates a matrix of kmer counts
    """
    col_names = {}
    row_names = {}

    # 11mers or shorter use precomputed lists of kmers
    chars = "ACGT"
    i = 0
    for item in itertools.product(chars, repeat=kmer_length):
        dna = "".join(item)
        revcomp = Seq.reverse_complement(dna)
        if revcomp < dna:
            dna = revcomp
        if not dna in col_names:
            col_names[dna] = i
            i += 1

    # Get a list of all files and reshape it to use with concurrent futures
    files = get_files_to_analyze(results_path)
    x = np.asarray(files)
    y = x.reshape((len(x),1))

    # Initialize the kmer matrix
    num_rows    = len(x)
    num_cols    = i
    kmer_matrix = np.zeros((num_rows,num_cols), dtype=matrix_dtype)

    # Use concurent futures to get multiple rows at the same time
    # Then place completed rows into the matrix and update the row dictionary
    row_index = 0
    with ProcessPoolExecutor(max_workers=num_threads) as ppe:
        for genomeid,temp_row in ppe.map(make_row, y, itertools.repeat(num_cols), itertools.repeat(col_names)):
            row_names[genomeid] = row_index
            kmer_matrix[row_index,:] = temp_row
            row_index += 1

    # Convert dict to array
    row_array = np.empty([num_rows], dtype='object')
    col_array = np.empty([num_cols], dtype='object')

    # Walk through row dictionary, place genome in correct index
    for key, index in row_names.items():
    	row_array[index] = key

    # Walk through col dictionary, place sequence in correct index
    for key, index in col_names.items():
    	col_array[index] = key

    df = pd.DataFrame(data = kmer_matrix, columns = col_array, index = row_array)

    return df
