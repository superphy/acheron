import pytest
from acheron.workflows import sequence_tools
import joblib
from collections import Counter

import xgboost as xgb
import numpy as np

def test_find_top_feats():
    """
    importances of test booster look like:
    {
    'AACCGACAGTG': 5.705407776666667,
    'AAAATGCACGA': 0.854700804,
    'ACGTAATGACG': 3.0228514700000004,
    'CGTATGCCGTC': 2.11447525,
    'ATAACAAAAGG': 2.3529501,
    'CGGAACTCTTC': 1.58123231
    }
    """

    bst = xgb.Booster({'nthread': 1})
    bst.load_model("data/acheron_test_samples/test.bst")
    feats = np.load("data/acheron_test_samples/feats.npy")

    importances = sequence_tools.find_top_feats(bst,feats,'xgb',3)

    assert len(importances)==3
    assert abs((1 - (importances[0][1] / 5.70540)))<0.0001
    assert abs((1 - (importances[1][1] / 3.02285)))<0.0001
    assert abs((1 - (importances[2][1] / 2.35295)))<0.0001

def test_save_query():
    bst = xgb.Booster({'nthread': 1})
    bst.load_model("data/acheron_test_samples/test.bst")
    feats = np.load("data/acheron_test_samples/feats.npy")
    importances = sequence_tools.find_top_feats(bst,feats,'xgb',3)

    path = "data/acheron_test_samples/test.query"

    sequence_tools.save_query(importances, path)

    query = open(path,'r')
    lines = query.readlines()

    assert lines[0].strip() == ">AACCGACAGTG"
    assert lines[1].strip() == "AACCGACAGTG"
    assert lines[2].strip() == ">ACGTAATGACG"
    assert lines[3].strip() == "ACGTAATGACG"

# dont add this to circleci
def test_biosamples_to_genomeids():
    biosamples = ['SAMN02640950','SAMN02699189','SAMN03894094','SAMN02911931']

    ids = sequence_tools.biosamples_to_genomeids(biosamples, 'Salmonella')

    assert ids == ['590.15667','82689.4','1620419.4','28901.767']
