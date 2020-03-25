#!/usr/bin/env python

import pandas as pd
import numpy as np

def make_mask(label_matrix, num_splits):
    """
    Takes in a label matrix and saves a mask of which labels are valid
    for each index, this allows the model to keep a sample when only a single
    columns is invalid, instead of discarding the sample
    """

def make_split(dataset, mask):
    """
    Goes through samples in wgs/raw, removes any considered invalid for the
    attribute, then splits and saves which genomes belong to which split
    """

def get_columns(label_df):
    """
    Takes column df and returns which attributes will be seen
    """

def apply_mask(df, mask):
    """
    Takes a dataset and mask, returns a dataset with only valid samples
    """

def load_data():
    """
    load requested dataset, call to apply masks
    """

def make_model():
    """

if __name__ == "__main__":
    """
    If we know that the inputs are satisfied (can be determined using
    the dry_run flag of acheron build model), we can manually call this script
    """
    pass
