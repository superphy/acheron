"""
Given a directory of jellyfish outputs, find the union of kmers
"""

import os, sys
import math
import glob
import time
import numpy as np
from Bio import Seq, SeqIO
from itertools import chain, repeat
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

def batch_genomes(genomes, num_batches):
    """
    populates 2D numpy array with len(rows)==num_batches in column major order
    """
    total_genomes = len(genomes)
    genomes_per_batch = math.ceil(total_genomes/num_batches)

    batches = np.empty([num_batches, genomes_per_batch], dtype=object)

    for i, genome in enumerate(genomes):
        batches[i%num_batches][i//num_batches] = genome

    return batches

def merge_masters(dataset, kmer_length, batches):

    feat_sets = {}
    for set_num in set_nums:
        feat_sets[set_num] = np.load(batches[set_num])

    all_feats = [feat_sets[i] for i in set_nums]

    master_mers = np.array(list(set(chain(*all_feats))))
    bool_mask = [len(master_mers[i]) for i in range(len(master_mers))]
    bool_mask = np.array([i == int(kmer_length) for i in bool_mask])
    master_mers = master_mers[bool_mask]

    return master_mers

def parse_fasta(genome, kmer_length, jf_path):
    jf_path = jf_path + genome
    current_multi_mers = []
    for record in SeqIO.parse(jf_path, "fasta"):
        kmer_seq = record.seq
        kmer_seq = kmer_seq._get_seq_str_and_check_alphabet(kmer_seq)
        if(len(kmer_seq) == int(kmer_length)):
            current_multi_mers.append(kmer_seq)
    return current_multi_mers

def find_union(kmer_length, dataset, batch_row, num_threads):
    # batch row is a single row of the batch_genomes array
    jf_path = "data/{}/wgs/{}mer_jellyfish_results/".format(dataset,kmer_length)
    genomes = ([files for r,d,files in os.walk(jf_path)][0])

    genomes = np.array(genomes)[np.array([i in batch_row for i in genomes])]

    multi_mers = []

    with ProcessPoolExecutor(max_workers = num_threads) as ppe:
        for current_multi_mer in ppe.map(parse_fasta, genomes,
        repeat(kmer_length), repeat(jf_path)):
            multi_mers.append(current_multi_mer)

    master_mers = np.array(list(set(chain(*multi_mers))))
    bool_mask = [len(master_mers[i]) for i in range(len(master_mers))]
    bool_mask = np.array([i == int(kmer_length) for i in bool_mask])
    master_mers = master_mers[bool_mask]

    return master_mers
