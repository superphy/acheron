import numpy as np
import pandas as pd

import pytest
import math

from acheron.workflows import supervised_model
from acheron import kmer, label

def setup_data():
    meta_path = "data/acheron_test_samples/labels/test_metadata.csv"

    label.build_module_label('acheron_test_samples','MIC', 'test_MIC',
    'AMP_MIC,CIP_MIC,SXT_MIC',meta_path, 'names')

    kmer.build_kmer_matrix('acheron_test_samples', 11, 1)

setup_data()

def test_make_mask():

    label_df = pd.read_pickle("data/acheron_test_samples/labels/test_MIC.df")
    label_df.at['extra_2','CIP'] = np.nan
    label_df.at['extra_2','CIP'] = math.nan
    label_df.at['extra_3','AMP'] = 'invalid'

    mask = supervised_model.make_mask(label_df, 2)

    for i, row in enumerate(mask.values):
        for j, col in enumerate(row):
            if j == 1:
                assert col == False
            elif j == 2 and i == 0:
                assert col == False
            elif j == 0 and i in [1,3]:
                assert col == False
            else:
                assert col == True

def test_make_split():
    valid_genomes = ['subSRR2407535','extra_1','extra_2','extra_3']
    data = np.array([[0,0],[0,1],[1,1],[1,0]])
    label_df = pd.DataFrame(data=data, index = valid_genomes, columns=['a','b'])
    mask = supervised_model.make_mask(label_df, 2)
    split_df = supervised_model.make_split(label_df, mask, 2, valid_genomes)

    for col in split_df.columns:
        assert(np.sum(split_df[col]) == 2)
