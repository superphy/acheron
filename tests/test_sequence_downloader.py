import pandas as pd
import numpy as np

from acheron.workflows import sequence_downloader

def test_find_seq_source():
    # NCBI
    df1 = pd.DataFrame(
        data=[['SAMN02640809','4','2'],
              ['SAMN02367860','8','16']],
        columns=['BioSample','MICA','MICB'])
    # PATRIC
    df2 = pd.DataFrame(
        data=[['SAMN02640806','4','2'],
              ['SAMN02640792','8','16']],
        columns=['biosample_accession','MICA','MICB'])

    df_dict = {'NCBI':df1,'PATRIC':df2}

    seq_sources = sequence_downloader.find_seq_sources(df_dict)

    assert list(seq_sources.keys()) == ['NCBI','PATRIC']

    assert list(seq_sources['NCBI']) == ['SAMN02640809','SAMN02367860']
    assert list(seq_sources['PATRIC']) == ['SAMN02640806','SAMN02640792']

def test_download_sources():
    # NCBI
    df1 = pd.DataFrame(
        data=[['SAMN02640809','4','2'],
              ['SAMN02367860','8','16']],
        columns=['BioSample','MICA','MICB'])
    # PATRIC
    df2 = pd.DataFrame(
        data=[['SAMN02640806','4','2'],
              ['SAMN02367860','8','16']],
        columns=['biosample_accession','MICA','MICB'])

    df_dict = {'NCBI':df1,'PATRIC':df2}

    seq_sources = sequence_downloader.find_seq_sources(df_dict)
    dl_lists = sequence_downloader.download_sources(seq_sources)

    assert list(seq_sources.keys()) == ['NCBI','PATRIC']

    assert list(dl_lists['NCBI']) == ['SAMN02640809','SAMN02367860']
    assert list(dl_lists['PATRIC']) == ['SAMN02640806']
