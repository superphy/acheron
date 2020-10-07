import pytest

from acheron import download
import numpy as np
import pandas as pd

def test_are_equal_mic():
    mics1 = ['32' ,'0.125','4','.01','>64','>= 64','chicken','32']
    mics2 = ['>32','0.12' ,'8','0.015','> 64','64','32','chicken']
    results = [True,True,False,True,True,True,['uncastable', 'lhs'],['uncastable', 'rhs']]

    for mics in zip(mics1,mics2,results):
        assert mics[2] == download.are_equal_mic(mics[0], mics[1])

def test_is_empty():
    mics = ['64','>32','',' ','\t','  ',np.nan,'NaN']
    results = [False,False,True,True,True,True,True,True]

    for mic in zip(mics,results):
        assert mic[1] == download.is_empty(mic[0])

def test_merge_antibiogram():
    df1 = pd.DataFrame(data=[['k1','a','1'],['k2','b','2']],
        columns=['BioSample','letter','number'])

    df2 = pd.DataFrame(data=[['k2','c','2'],['k3','c','3']],
        columns=['BioSample','letter','number'])

    merged_df = download.merge_antibiogram(df1,df2)

    expected_df = pd.DataFrame(data=[['k1','a','1'],['k2','invalid','2'],['k3','c','3']],
        columns=['BioSample','letter','number'])

    assert(merged_df.equals(expected_df))
