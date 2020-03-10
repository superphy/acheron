import pytest
import numpy as np
import os

from acheron.workflows import find_common_kmers

test_list = ['a','b','c','d','e','f','g']
test_set = 'acheron_test_samples'

def test_batch_genomes():
    col_batches = find_common_kmers.batch_genomes(test_list, 3, 'col')
    assert np.array_equal(col_batches, np.array(
    [['a','d','g'],
    ['b','e',None],
    ['c','f',None]]))

    row_batches = find_common_kmers.batch_genomes(test_list, 3, 'row')
    assert np.array_equal(row_batches, np.array(
    [['a','b','c'],
    ['d','e','f'],
    ['g',None,None]]))

    long_batches = find_common_kmers.batch_genomes(test_list, 2, 'row')
    assert np.array_equal(long_batches, np.array(
    [['a','b','c','d'],
    ['e','f','g',None]]))

def test_merge_masters():
    paths = ["data/{}/seq_1.npy".format(test_set),
    "data/{}/seq_2.npy".format(test_set)]

    np.save("data/{}/seq_1.npy".format(test_set),
    np.array(['AAAAAAAAAAA','CCCCCCCCCCC']))
    np.save("data/{}/seq_2.npy".format(test_set),
    np.array(['AAAAAAAAAAA','trash']))

    master_mers = find_common_kmers.merge_masters(test_set, 11, paths,
    len(paths))

    assert len(master_mers) == 2
    assert 'AAAAAAAAAAA' in master_mers
    assert 'CCCCCCCCCCC' in master_mers
    assert isinstance(master_mers[0], str)

    for path in paths:
        os.system("rm "+path)


def test_parse_fasta():
    multi_mers = find_common_kmers.parse_fasta('test',4,
    "data/{}/wgs/tests/".format(test_set))
    assert len(multi_mers) == 3
    for seq in ['seq1', 'seq2','seq3']:
        assert seq in multi_mers
    assert isinstance(multi_mers[0], str)

"""
def test_find_union():
"""
