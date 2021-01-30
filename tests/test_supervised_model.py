import numpy as np
import pandas as pd

import pytest
import math

from acheron.workflows import supervised_model
from acheron import kmer, label

import joblib

def setup_data():
    meta_path = "data/acheron_test_samples/labels/test_metadata.csv"

    label.build_module_label('acheron_test_samples','MIC', 'test_MIC',
    'AMP_MIC,CIP_MIC,SXT_MIC',meta_path, 'names')

    kmer.build_kmer_matrix('acheron_test_samples', 11, 1, 'none')

setup_data()

def test_is_valid_type():
    assert supervised_model.is_valid_type(16) == True
    assert supervised_model.is_valid_type(16.0000) == True
    assert supervised_model.is_valid_type("16") == True
    assert supervised_model.is_valid_type("16.0000") == True
    assert supervised_model.is_valid_type(np.nan) == False
    assert supervised_model.is_valid_type('NaN') == False
    assert supervised_model.is_valid_type('invalid') == False
    assert supervised_model.is_valid_type([1,2,3]) == False

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

def test_apply_mask():
    features = pd.DataFrame(
        data = [[0,1,2],[3,4,5],[6,7,8]],
        index = ['SAMN00000001','SAMN00000002','SAMN00000003'],
        columns = ['first','second','third'])

    labels = pd.Series(
        data = ['<=1.0000',np.nan,'16.0000'],
        index = ['SAMN00000002','SAMN00000001','SAMN00000003'])

    mask = pd.Series(
        data = [True,False,True],
        index = ['SAMN00000002','SAMN00000001','SAMN00000003'])

    new_feats, new_labels = supervised_model.apply_mask(features, labels, mask)

    assert new_feats.shape == (2,3)
    assert new_feats['first']['SAMN00000002'] == 3
    assert new_feats['second']['SAMN00000003'] == 7

    assert len(new_labels) == 2
    assert new_labels['SAMN00000002'] == '<=1.0000'
    assert new_labels['SAMN00000003'] == '16.0000'

def test_select_features():
    # SelectKBest will have a temper tantrum deciding
    # between columns 2 and 3 here
    import warnings
    warnings.filterwarnings('ignore')

    features = pd.DataFrame(
        data = [
        [1,1,0],
        [0,1,0],
        [1,0,1]],
        index = ['SAMN00000001','SAMN00000002','SAMN00000003'],
        columns = ['first','second','third'])

    labels = [0,0,1]

    new_features = supervised_model.select_features(features, labels, 2, False)

    assert len(new_features.columns) == 2
    assert new_features.columns[0] == 'second'
    assert new_features.columns[1] == 'third'

    assert new_features.shape == (3,2)

    new_columns = new_features.columns

    newer_features = supervised_model.select_features(features,labels,2,new_columns)

    assert len(newer_features.columns) == 2
    assert newer_features.columns[0] == 'second'
    assert newer_features.columns[1] == 'third'

    assert newer_features.shape == (3,2)

def test_evaluate_model():
    # [1 ,1,2,2,8]
    prediction = [0,0,1,1,3]
    # [32,1,2,4,16]
    actual =     [5,0,1,2,4]

    test_encoder = joblib.load('data/acheron_test_samples/labels/test_MIC_encoder.pkl')
    results_df = supervised_model.evaluate_model(prediction, actual, 'XGB', [0,1], 'AMP', test_encoder['AMP'])

    cols = ['Accuracy','Within 1 Dilution','Very Major Error','Major Error','Non Major Error','Correct']
    accs = [0.4       ,0.8                ,0.2               ,0.0          ,0.4              ,0.4]

    for col, acc in zip(cols, accs):
        assert results_df[col][0] == acc
