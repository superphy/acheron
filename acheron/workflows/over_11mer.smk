from clean import get_files_to_analyze
import numpy as np
import pandas as pd

import math
import os

# Location of the raw genomes
RAW_GENOMES_PATH = "data/{}/wgs/raw/".format(config['dataset'])

# Kmer length that you want to count
KMER_SIZE = config['kmer_length']

# Data type of the resulting kmer matrix. Use uint8 if counts are
# all under 256. Else use uint16 (kmer counts under 65536)
# note: counts should be under 256 anyways, so if thats not the case,
# check your sequences
MATRIX_DTYPE = 'uint8'

# it is assumed there is enough ram for the union tests to handle MAX_GENOMES
# This is used supercluster nodes with 1TB RAM, lower substantially
# if using a normal computer
MAX_GENOMES = 1000
GENOMES = get_files_to_analyze(RAW_GENOMES_PATH)
NUM_GENOMES = len(GENOMES)
UNION_SPLITS = [i+1 for i in range(math.ceil(NUM_GENOMES/MAX_GENOMES))]

MATRIX_SPLITS = [i+1 for i in range(math.ceil(NUM_GENOMES/128))]

ids, = glob_wildcards(RAW_GENOMES_PATH+"{id}.fasta")

# from sub_11mer.smk, we will use clean_fastas, count_kmers, dump_kmers
subworkflow kmer:
    workdir:
        os.getcwd()
    snakefile:
        "sub_11mer.smk"

rule all:
    input:
        "data/{}/features/{}mer_matrix.df".format(config['dataset'],
        config['kmer_length'])

# for subsets of genomes, adding this as a seperate step saves on ram
# as we cannot fit all genomes in memory at a single time, so we batch the kmers
# into groups of at most MAX_GENOMES, union those, and then union the unions

# builds 2D numpy array, each row is a batch with len(colums)<=MAX_GENOMES
rule batch_genomes_for_union:
    input:
        kmer(expand("data/"+config['dataset']+"/wgs/"+
        str(config['kmer_length'])+"mer_jellyfish_results/{id}.fa",id=ids))
    output:
        "data/{}/wgs/master_{}mers/batches.npy".format(config['dataset'],
        config['kmer_length'])
    run:
        from acheron.workflows import find_common_kmers
        batches = find_common_kmers.batch_genomes(GENOMES, len(UNION_SPLITS), 'col')
        np.save(output[0], batches)

# find union of kmers in each batch
rule batch_kmers:
    input:
        "data/{}/wgs/master_{}mers/batches.npy".format(config['dataset'],
        config['kmer_length'])
    output:
        "data/"+config['dataset']+"/wgs/master_"+str(config['kmer_length'])+
        "mers/kmer_subset_{set_num}_of_"+str(len(UNION_SPLITS))+'.npy'
    run:
        from acheron.workflows import find_common_kmers
        batches = np.load("data/{}/wgs/master_{}mers/batches.npy".format(
        config['dataset'], config['kmer_length']), allow_pickle=True)
        batch_row = batches[int(wildcards.set_num)-1]

        union_of_batch = find_common_kmers.find_union(config['kmer_length'],
        config['dataset'], batch_row, config['cores'])
        np.save(output[0], union_of_batch)

# merge all the unions from the above step
rule merge_kmer_batches:
    input:
        expand("data/"+config['dataset']+"/wgs/master_"+
        str(config['kmer_length'])+"mers/kmer_subset_{set_num}_of_"+
        str(len(UNION_SPLITS))+'.npy',set_num = UNION_SPLITS)
    output:
        "data/"+config['dataset']+"/wgs/master_"+str(config['kmer_length'])+
        "mers/all_kmers.npy"
    run:
        from acheron.workflows import find_common_kmers
        # im hoping here that input is a list of inputs that i can use,
        # gonna have to test this one
        master_mers = find_common_kmers.merge_masters(config['dataset'],
        config['kmer_length'], input, len(UNION_SPLITS))
        np.save(output[0], master_mers)


# now that we know what kmers we will eventually see, we can make parts
# of the dataframe

rule split_matrices:
    output:
        "data/{}/wgs/master_{}mers/matrix_splits.npy".format(config['dataset'],
        config['kmer_length'])
    run:
        from acheron.workflows import find_common_kmers
        batches = find_common_kmers.batch_genomes(GENOMES, len(MATRIX_SPLITS), 'col')
        np.save(output[0], batches)

# we start by building the matrix in sets, for easier distribution over
# compute cluster nodes
rule build_sub_matrices:
    input:
        "data/"+config['dataset']+"/wgs/master_"+str(config['kmer_length'])+
        "mers/all_kmers.npy",

        kmer(expand("data/"+config['dataset']+"/wgs/"+
        str(config['kmer_length'])+"mer_jellyfish_results/{id}.fa",id=ids)),

        "data/{}/wgs/master_{}mers/matrix_splits.npy".format(config['dataset'],
        config['kmer_length'])
    threads:
        16
    output:
        "data/"+config['dataset']+"/wgs/master_"+str(config['kmer_length'])+
        "mers/sub_df_{matrix_num}_of_"+str(len(MATRIX_SPLITS))+".df"
    run:
        #multi_mer_matrix
        from over_11mer_matrix import make_matrix
        df = make_matrix(config['dataset'], config['kmer_length'], MATRIX_DTYPE,
        config['cores'], wildcards.matrix_num)
        df.to_pickle(output[0])


rule merge_sub_matrices:
    input:
        expand("data/"+config['dataset']+"/wgs/master_"+
        str(config['kmer_length'])+"mers/sub_df_{matrix_num}_of_"
        +str(len(MATRIX_SPLITS))+".df", matrix_num = MATRIX_SPLITS)
    output:
        "data/{}/features/{}mer_matrix.df".format(config['dataset'],
        config['kmer_length'])
    run:
        # again, were hoping that the input is actually a list of inputs,(check)
        #final_matrix = pd.concat([pd.read_pickle(i) for i in input])

        pre_filter = False # toggle this to save RAM at the expensive of storage

        if not pre_filter:

            final_matrix = pd.read_pickle(input[0])
            for i, mat in enumerate(input):
                if i != 0:
                    new_slice = pd.read_pickle(mat)
                    final_matrix = pd.concat([final_matrix, new_slice])
                    del new_slice
                    gc.collect()
            final_matrix.to_pickle(output[0])

        else:
            # This is because we dont have labels during feature generation, but we
            # need it in this specific case to pre-filter poor quality k-mers
            print("This expansion is explicity written for MIC's, it will not work with other modules")
            samples = list(pd.read_pickle(input[0]).index)
            for drug in ['AMP','AMC','AZM','CHL','CIP','CRO','FIS','FOX','GEN','NAL','SXT','TET','TIO','STR','KAN']:
                label_df = pd.read_pickle("data/{}/labels/AMR_MIC.df".format(config['dataset']))[[drug]]
                label_df.reindex(samples)

                with open("data/{}/labels/AMR_MIC_encoder.pkl".format(config['dataset']), 'rb') as unpickler:
                    encoder = pickle.load(unpickler)
                labels = [encoder[i] for i in label_df[drug]]

                #             large_feat_select(file_paths, RAM_denom, num_feats, labels)
                drug_matrix = large_feat_select(input, 5, 10000000, labels)

                drug_matrix.to_pickle("data/{}/features/{}mer_matrix{}.df".format(config['dataset'],
                config['kmer_length'], drug))

                # NOTE: This is intented to fail the snakemake, you need to check manually for the output
