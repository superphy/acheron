import pytest
import math
import numpy as np
import pandas as pd
import os

from acheron.workflows import labeler
meta_path = "data/acheron_test_samples/labels/test_metadata.csv"

def test_load_metadata():
    test_df = labeler.load_metadata(meta_path)

    no_ext = meta_path.split('.')[0]
    test_df.to_csv(no_ext+'.tsv', sep='\t')
    test_df.to_excel(no_ext+'.xlsx')
    test_df.to_pickle(no_ext+'.df')

    csv_df = labeler.load_metadata(no_ext+'.csv')
    tsv_df = labeler.load_metadata(no_ext+'.tsv')
    excel_df = labeler.load_metadata(no_ext+'.xlsx')
    pandas_df = labeler.load_metadata(no_ext+'.df')

    for df_type in [csv_df, tsv_df, excel_df, pandas_df]:
        for header in ['names','AMP_MIC','CIP_MIC','AMP_SIR','CIP_SIR']:
            assert len(df_type[header]) == 4
            assert header in df_type.columns

        assert np.sum([math.isnan(i) for i in df_type['CIP_MIC']]) == 2

    for ext in ['.df','.tsv','.xlsx']:
        os.system("rm "+no_ext+ext)

def test_get_valid_columns():
    columns = ['col_1','col_2','col_3']

    npy_path = "data/acheron_test_samples/labels/columns.npy"
    np.save(npy_path, columns)

    as_npy = labeler.get_valid_columns(npy_path)
    as_list = labeler.get_valid_columns(','.join(columns))

    for col_type in [as_npy, as_list]:
        assert len(columns) == len(col_type)
        for i in range(len(columns)):
            assert columns[i] == col_type[i]

    os.system('rm '+npy_path)

def test_build_encoder():
    dataset = 'acheron_test_samples'
    name = 'mics'
    mic_encoder = labeler.build_encoder(dataset,name,'MIC', 'Salmonella')
    assert len(mic_encoder.keys()) == 3
    assert mic_encoder['AMP']['<=1.0000'] == 0
    assert mic_encoder['CIP']['0.0600'] == 2

    name = 'test_SIR'
    general_encoder = labeler.build_encoder(dataset,name,'none', 'Salmonella')
    assert len(general_encoder.keys()) == 2
    assert general_encoder['AMP_SIR']['S'] == 0
    assert general_encoder['CIP_SIR']['R'] == 1
