import pytest

from acheron import kmer
dataset = 'acheron_test_samples'
test_fasta = 'data/acheron_test_samples/wgs/raw/subSRR2407535.fasta'

def test_get_genomes():
    sequences = kmer.get_genomes(dataset)
    assert len(sequences) == 1
    assert sequences[0] == test_fasta

def test_build_kmer_matrix():
    import pandas as pd
    import numpy as np
    import os

    #11mer
    kmer.build_kmer_matrix(dataset, 11, 1, 'none', 'False',1,5,False,'AMR_MIC')
    df = pd.read_pickle("data/{}/features/11mer_matrix.df".format(dataset))
    assert len(df.values) == 1
    assert len(df.values[0]) == (4**11)/2
    assert np.sum(df.values) == 4977
    assert len(df.columns[0]) == 11
    assert df['AAAAAAAAAAA'].values[0] == 0
    assert df['AAAAAATACTT'].values[0] == 1
    os.system("rm data/acheron_test_samples/features/11mer_matrix.df")
    os.system("rm -r data/acheron_test_samples/wgs/11mer_jellyfish*")
    os.system("rm -r data/acheron_test_samples/wgs/clean/")

    #31mer
    kmer.build_kmer_matrix(dataset, 31, 1, 'none', 'False',1,5,False,'AMR_MIC')
    df = pd.read_pickle("data/{}/features/31mer_matrix.df".format(dataset))
    assert len(df.values) == 1
    assert len(df.values[0]) == 4294
    assert np.sum(df.values) == 4917
    assert len(df.columns[0]) == 31
    assert df['AGCTCCATGAGTAAATCACTTCACCTACATG'].values[0] == 1
    os.system("rm data/acheron_test_samples/features/31mer_matrix.df")
    os.system("rm -r data/acheron_test_samples/wgs/31mer_jellyfish*")
    os.system("rm -r data/acheron_test_samples/wgs/master_31mers/")
    os.system("rm -r data/acheron_test_samples/wgs/clean/")
