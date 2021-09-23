import pytest
from acheron.workflows import over_11mer_matrix

def test_get_feat_splits():
    feats = ['a','b','c','d','e','f']
    splits = over_11mer_matrix.get_feat_splits(feats, 3)
    assert splits[0] == ['a','b']
    assert splits[1] == ['c','d']
    assert splits[2] == ['e','f']

    feats = ['a','b','c','d','e','f','g']
    splits = over_11mer_matrix.get_feat_splits(feats, 3)
    assert splits[0] == ['a','b','c']
    assert splits[1] == ['d','e','f']
    assert splits[2] == ['g']


def test_large_feat_select():
    paths = []

    for i in range(3):
        paths.append("data/acheron_test_samples/features/sub{}_of_3.df".format(i+1))

    final_df = over_11mer_matrix.large_feat_select(paths,3,2,[0,0,2,2,1,1])

    assert list(final_df.columns) == ['c','d']

    assert list(final_df.index) == ['one','two','three','four','five','six']
