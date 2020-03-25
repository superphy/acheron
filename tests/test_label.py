import pytest
import pandas as pd
from acheron import label
import os

meta_path = "data/acheron_test_samples/labels/test_metadata.csv"

def test_build_custom_label():
    label.build_custom_label('acheron_test_samples','test_SIR', 'AMP_SIR,CIP_SIR',
    meta_path, 'names')

    df = pd.read_pickle("data/acheron_test_samples/labels/test_SIR.df")

    assert len(df.index) == 4
    assert len(df.columns) == 2
    assert df.index[0] == 'subSRR2407535'
    assert df['AMP_SIR']['subSRR2407535'] == 'S'

    os.system("rm data/acheron_test_samples/labels/test_SIR.df")

def test_build_module_label():
    label.build_module_label('acheron_test_samples','MIC', 'test_MIC',
    'AMP_MIC,CIP_MIC,SXT_MIC',meta_path, 'names')

    df = pd.read_pickle("data/acheron_test_samples/labels/test_MIC.df")

    assert len(df.index) == 4
    assert len(df.columns) == 3
    assert df['AMP']['subSRR2407535'] == '>=32.0000'
    assert df['AMP']['extra_1'] == 'invalid'
    assert df['AMP']['extra_2'] == '>=32.0000'
    assert df['AMP']['extra_3'] == 'invalid'
    assert df['SXT']['subSRR2407535'] == '<=0.1250'
    assert df['SXT']['extra_1'] == '>=64.0000'

    os.system("rm data/acheron_test_samples/labels/test_MIC.df")
