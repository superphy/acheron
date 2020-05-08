#!/usr/bin/env python

import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from collections import Counter
import math

def is_valid_type(val):
    # covers core floats, numpy floating and complex floating
    if isinstance(val, (float, np.inexact)):
        if math.isnan(val) or np.isnan(val):
            return False
        else:
            return True
    # covers core ints, numpy ints (signed unsigned), longs, shorts, bytes
    # TODO check for bytes later, as they do NOT like math eqns when exponential
    elif isinstance(val, (int, np.integer)):
        return True
    # covers core and np strings, excludes unicode
    elif isinstance(val, (str, np.str_)):
        if val == 'invalid':
            return False
        else:
            return True
    else:
        return False

def make_mask(label_matrix, k):
    """
    Takes in a label matrix and saves a mask of which labels are valid
    for each index, this allows the model to keep a sample when only a single
    columns is invalid, instead of discarding the sample.
    If doing kfold cross validation (k != 1), there must be at least k samples
    in each class for that class to be considered valid
    """
    mask = np.zeros((label_matrix.values.shape), dtype='bool')

    invalid_classes = {}

    # if any class has fewer than k samples, it is marked as invalid
    for attribute in label_matrix.columns:
        class_counter = Counter(label_matrix[attribute])
        invalid_classes[attribute] = []
        for key in class_counter:
            if class_counter[key] < k:
                invalid_classes[attribute].append(key)

    for i, row in enumerate(label_matrix.values):
        for j, col in enumerate(row):
            if not is_valid_type(col):
                continue
            if col in invalid_classes[label_matrix.columns[j]]:
                continue
            mask[i,j] = True

    return pd.DataFrame(data=mask, columns = label_matrix.columns,
        index = label_matrix.index)

def make_split(label_matrix, mask, k, samples):
    """
    Takes a label matrix, splits it according to k fold cross validation
    but only for valid samples. Produces a random split each time (stratified)
    """
    assert(2<=k<=255) # if exceeding 255, change dtype below to uint16
    split_df =  pd.DataFrame(data=np.zeros(label_matrix.shape),
        columns=label_matrix.columns, index = label_matrix.index, dtype='uint8')

    for col in label_matrix.columns:
        # which labels are valid in this specific column
        valid_labels = label_matrix[col].values[mask[col].values]
        # matching sample name for each i in valid_labels
        valid_samples = label_matrix.index.values[mask[col].values]

        if len(valid_samples) == 0:
            print("All samples in column "+col+" are invalid, skipping split")
            continue

        # we also need to factor in that we only have the samples in /samples,
        # where a datasheet might have thousands of valid, but extra datapoints
        seen_bool_mask = np.array([i in samples for i in valid_samples])
        final_labels = valid_labels[seen_bool_mask]
        final_samples = valid_samples[seen_bool_mask]

        # at this point, we only have labels and samples that are eligible
        # for machine learning
        skf = StratifiedKFold(n_splits=k, shuffle=True)

        num_samples = len(final_samples)
        splits = enumerate(skf.split(np.zeros((num_samples,k)),final_labels))
        for i, split in splits:
            # pull the array of values assigned to the testing set,
            # label these genomes as per the fold they belong to
            for sample in final_samples[split[1]]:
                split_df.at[sample,col] = i

    return split_df


def load_data():
    """
    load requested dataset, call to apply masks
    """

def make_model():
    """
    """

if __name__ == "__main__":
    """
    If we know that the inputs are satisfied (can be determined using
    the dry_run flag of acheron build model), we can manually call this script
    """
    pass
