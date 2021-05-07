import pandas as pd
import numpy as np
import logging
import pickle
import yaml

import os,sys

from mic import MICPanel
from acheron.workflows import labeler

def transform(input, log, output,
    slice_cols=["run", "MIC_AMP", "MIC_AMC", "MIC_FOX", "MIC_CRO",
    "MIC_TIO", "MIC_GEN", "MIC_FIS", "MIC_SXT", "MIC_AZM",
    "MIC_CHL", "MIC_CIP", "MIC_NAL", "MIC_TET"]):

    logging.basicConfig(filename=log[0],level=logging.DEBUG)
    logging.info('MIC binning')

    micsdf = labeler.load_metadata(input[0])

    if slice_cols[0] == 'all':
        key = slice_cols[1]
        all_seen_cols = micsdf.columns
        all_seen_mic = [i for i in all_seen_cols if 'MIC' in i]

        slice_cols = [key]+all_seen_mic

    micsdf = micsdf[slice_cols]
    micsdf = micsdf.set_index(config['key'])

    logging.debug('Using classes defined in {}'.format(input[1]))
    with open(input[1], 'r') as infh:
        mic_class_labels = yaml.load(infh)

        classes = {}
        class_orders = {}
        for col in micsdf:
            logging.debug('Creating MIC panel for {}'.format(col))
            if col[:4] == 'MIC_':
                drug = col.replace('MIC_', '')
            elif col[-4:] == '_MIC':
                drug = col.replace('_MIC', '')
            else:
                drug = col

            mic_series = micsdf[col]
            class_list = []
            panel = MICPanel()

            # Initialize MIC bins for this drug
            panel.build(mic_class_labels[drug])

            # Iterate through MIC values and assign class labels
            logging.debug('MIC values will be mapped to: {}'.format(panel.class_labels))

            i = 0
            for m in mic_series:
                i += 1
                mlabel, isna = panel.transform(m)

                if isna:
                    class_list.append(np.nan)
                else:
                    class_list.append(mlabel)

            logging.info("Invalid MICs found in dataset: {}".format(', '.join(panel.invalids)))

            classes[drug] = pd.Series(class_list, index=micsdf.index)
            class_orders[drug] = panel.class_labels

            logging.debug("Final MIC distribution:\n{}".format(classes[drug].value_counts()))

        c = pd.DataFrame(classes)

        cfile = output[0]
        cofile = output[1]

        c.to_pickle(cfile)
        pickle.dump(class_orders, open(cofile, "wb"))


rule all:
  input:
    "data/{}/labels/{}.df".format(config['dataset'],config['name']),
    "data/label_modules/mic/{}_mic_class_order_dict.pkl".format(config['pathogen'])

rule make_mic_df:
    input:
        config['path'],
        "data/label_modules/mic/{}_class_ranges.yaml".format(config['pathogen'])
    log:
        "logs/{}/{}_mic.log".format(config['dataset'], config['name'])
    output:
        "data/{}/labels/{}.df".format(config['dataset'],config['name']),
        "data/label_modules/mic/{}_mic_class_order_dict.pkl".format(config['pathogen'])
    run:
        if config['columns'].lower() == 'all':
            headers =['all',config['key']]
        else:
            headers = labeler.get_valid_columns(config['columns'])
            headers = np.concatenate([headers, [config['key']]])
        transform(input, log, output, slice_cols=headers)
