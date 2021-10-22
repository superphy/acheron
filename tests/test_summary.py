import pytest
from acheron import summary
import pandas as pd

def test_load_steinkey2021():
    res = summary.load_steinkey2021()

    assert res.shape[0] > 4700
    assert len(res.columns) == 10

def test_add_results():
    try:
        pd.read_pickle("results/model=ANN_train=grdi_test=none_validate=none_feats=1000_type=11mer_hyp=False_cvfolds=5_attribute=AMC_trial=0/summary.df")
    except:
        pytest.skip("Skipping test that requires results data")

    full_res = summary.add_results(summary.load_steinkey2021())

    print(tests)
    assert False
