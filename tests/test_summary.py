import pytest
from acheron import summary
import pandas as pd

def test_load_steinkey2021():
    try:
        pd.read_pickle("results/model=ANN_train=grdi_test=none_validate=none_feats=1000_type=11mer_hyp=False_cvfolds=5_attribute=AMC_trial=0/summary.df")
    except:
        pytest.skip("Skipping test that requires results data")
    res = summary.load_steinkey2021()

    assert res.shape[0] > 4700
    assert len(res.columns) == 10

def test_add_results():
    try:
        pd.read_pickle("results/model=ANN_train=grdi_test=none_validate=none_feats=1000_type=11mer_hyp=False_cvfolds=5_attribute=AMC_trial=0/summary.df")
    except:
        pytest.skip("Skipping test that requires results data")

    full_res = summary.add_results(summary.load_steinkey2021())

    assert full_res.shape[0] > 4700

def test_weigh_jobs_per_slurm():
    df = pd.DataFrame(
        columns=['fake1','SLURM_JOBID','Time (m)'],
        data=[
            [0,1,64],
            [0,2,64],
            [0,2,64],
            [0,3,64],
            [0,3,64],
            [0,3,64],
            [0,3,64]])

    df = summary.weigh_jobs_per_slurm(df)

    for i,j in zip (df['Time (m)'],[64,32,32,16,16,16,16]):
        assert abs(i/j-1)<0.01

def test_summarize():
    try:
        pd.read_pickle("results/model=ANN_train=grdi_test=none_validate=none_feats=1000_type=11mer_hyp=False_cvfolds=5_attribute=AMC_trial=0/summary.df")
    except:
        pytest.skip("Skipping test that requires results data")
    results = summary.load_steinkey2021()
    results = summary.add_results(results)
    results = summary.weigh_jobs_per_slurm(results)

    summ = summary.summarize(results,'steinkey2021', 'table')
