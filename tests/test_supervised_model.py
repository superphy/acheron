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
    'AMP_MIC,CIP_MIC,SXT_MIC',meta_path, 'names', 'Test')

    kmer.build_kmer_matrix('acheron_test_samples', 11, 1, 'none', 'False',1,5,False,'AMR_MIC')

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


    dup_label = pd.DataFrame(
        data = [[0,0],
                [0,1],
                [1,0],
                [1,1],
                [0,0]],
        index = ['SAMN00000001','SAMN00000002','SAMN00000003','SAMN00000004','SAMN00000004'],
        columns = ['first','second'])

    mask = supervised_model.make_mask(dup_label, 2)
    split = supervised_model.make_split(dup_label, mask, 2, ['SAMN00000001','SAMN00000002','SAMN00000003','SAMN00000004'])

    assert len(split.index) == 4

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

    cols = ['Supports','Accuracy','Within 1 Dilution','Very Major Error','Major Error','Non Major Error','Correct']
    accs = [5         ,0.4       ,0.8                ,0.2               ,0.0          ,0.4              ,0.4]

    for col, acc in zip(cols, accs):
        assert results_df[col][0] == acc

def test_split_data():
    features = pd.DataFrame(
        data = [[0,1,2],[3,4,5],[6,7,8],[9,10,11],[12,13,14]],
        index = ['SAMN00000001','SAMN00000002','SAMN00000003','SAMN00000004','SAMN00000005'],
        columns = ['AAAAAAAAAAA','CCCCCCCCCCC','TTTTTTTTTTT'])

    labels = pd.DataFrame(
        data = [['1','2','4'],['8','8','2'],['8','4','1'],['32','32','16'],['4','2','1']],
        index = ['SAMN00000001','SAMN00000002','SAMN00000003','SAMN00000004','SAMN00000005'],
        columns = ['AMP','AMC','TET'])

    split = pd.DataFrame(
        data = [[0,0,0],[1,1,1],[2,2,2],[3,3,3],[4,4,4]],
        index = ['SAMN00000001','SAMN00000002','SAMN00000003','SAMN00000004','SAMN00000005'],
        columns = ['AMP','AMC','TET'])

    attribute = 'AMP'
    fold = 1
    cv_folds = 5

    # cv split
    x_train, y_train, x_test, y_test = supervised_model.split_data(features, labels[attribute], split, attribute, False, fold, cv_folds)
    assert x_train.shape == (4, 3)
    assert y_train.shape == (4,)
    assert x_test.shape == (1, 3)
    assert y_test.shape == (1,)

    assert 'SAMN00000001' in x_test.index

    # nested cv split
    x_train, y_train, x_val, y_val, x_test, y_test = supervised_model.split_data(features, labels[attribute], split, attribute, True, fold, cv_folds)
    assert x_train.shape == (3, 3)
    assert y_train.shape == (3,)
    assert x_test.shape == (1, 3)
    assert y_test.shape == (1,)
    assert x_val.shape == (1, 3)
    assert y_val.shape == (1,)

    assert 'SAMN00000001' in x_test.index
    assert 'SAMN00000002' in x_val.index
    assert 'SAMN00000003' in x_train.index

def test_mean_summaries():
    cols = ["Supports", "Accuracy", "Error Rate"]

    t1 = pd.DataFrame(
    data=[[10, 0.8, 0.5]],
    columns= cols)

    t2 = pd.DataFrame(
    data=[[20, 0.5, 0.2]],
    columns= cols)

    t3 = supervised_model.mean_summaries([t1,t2])

    assert len(t3.columns) == 3
    for i in cols:
        assert i in t3.columns

    for c,v in zip(cols, [30,0.6,0.3]):
        assert t3[c][0] == v

    t4 = supervised_model.mean_summaries([t1])
    assert len(t4.columns) == 3

    for c,v in zip(cols, [10,0.8,0.5]):
        assert t4[c][0] == v

def test_mean_prec_recall():
    cols = ["Precision", "Recall", "F-Score", "Supports"]
    mics = ['<=1.0000','2.0000']

    df1 = pd.DataFrame(
    data = [[0.8, 0.9, 0.6, 200],
            [0.4, 0.5, 0.6, 20]],
    columns = cols,
    index = mics)

    df2 = pd.DataFrame(
    data = [[0.5, 0.6, 0.3, 400],
            [0.2, 0.25, 0.3, 20]],
    columns = cols,
    index = mics)

    df3 = supervised_model.mean_prec_recall([df1,df2])

    assert df3.shape == (2,4)
    assert len(df3.index) == 2

    for mic in mics:
        assert mic in df3.index

    for col, res in zip(cols,[0.6,0.7,0.4,600]):
        assert df3.loc['<=1.0000', col] == res

    for col, res in zip(cols,[0.3,0.375,0.45,40]):
        assert df3.loc['2.0000', col] == res
