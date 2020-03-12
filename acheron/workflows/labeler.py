import numpy as np
import pandas as pd

def load_metadata(path):
    if path[-4:] == '.csv':
        df = pd.read_csv(path)
    elif path[-4:] == '.tsv':
        df = pd.read_csv(path, sep='\t')
    elif path[-4:] == 'xlsx':
        df = pd.read_excel(path)
    elif path[-4:] == '.pkl' or path[-3:] == '.df':
        df = pd.read_pickle(path)
    else:
        raise Exception("Label metadata file must end in .csv, .tsv, xlsx, or \
        .pkl, .df if pickled, {} was given".format(path.split('/')[-1]))
    return df


def get_valid_columns(columns):
    # check if what the use passed for the columns argument is valid

    is_short = False
    is_numpy = False
    try:
        potential_extension = columns[-4:]
        if potential_extension == '.npy':
            is_numpy = True
    except:
        is_short = True

    if is_numpy:
        headers = np.load(columns)
        if len(headers.shape)!=1:
            raise Exception("numpy passed as columns needs to be 1 dimensional")
    else:
        try:
            headers = np.array(columns.split(','))
        except:
            raise Exception("columns need to be a in comma seperated format, \
            e.g. 'col_1,col_2,col_3'")

    return headers
