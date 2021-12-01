import pandas as pd
import numpy as np
from collections import Counter
import skbio
import sys
import pickle

drugs = ['AMC','AMP','AZM','FOX','TIO','CRO','CHL','CIP','GEN','NAL','STR','FIS',
    'TET','SXT','KAN']

def make_simpsons(datasets):

    res_df = pd.DataFrame(data = np.zeros((15,2), dtype = 'float'), columns = ['salm_amr','grdi'],index = drugs)
    for dataset in datasets:
        df = pd.read_pickle("data/{}/labels/AMR_MIC.df".format(dataset))

        with open("data/{}/labels/AMR_MIC_encoder.pkl".format(dataset),'rb') as fh:
            order_dict = pickle.load(fh)

        for drug in drugs:
            if dataset == "grdi" and drug == 'FIS':
                continue
            counts = Counter(df[drug])
            ordered_counts = {}
            for mic in order_dict[drug]:
                ordered_counts[str(mic)] = counts[mic]

            res_df[dataset][drug] = skbio.diversity.alpha.simpson(list(ordered_counts.values()))
    return res_df

if __name__ == "__main__":
    """
    Saves a simpsons diversity df's for 2 datasets
    for use with the figure generation in `acheron summary --media figure`
    """
    df = make_simpsons(['salm_amr','grdi'])
    df.to_pickle("data/simpsons_diversity.df")
