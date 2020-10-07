import numpy as np
import pandas as pd
import os
import pickle

from acheron.workflows import labeler

def build_module_label(dataset, module, name, columns, path, key):
    # call to snakemake subworkflow in workflows
    print("Building labels named {} for dataset {} for columns {} in {} on \
    key={} based on module {}".format(
        name, dataset, columns, path, key, module))

    if module =='MIC':
        workflow_smk = "acheron/workflows/mic_labels.smk"
    else:
        raise Exception("Only supported modules are [MIC], or switch to custom")


    os.system("snakemake -s {} -j {} \
    --config path={} dataset={} columns={} key={} name={}".format(
    workflow_smk, 1, path, dataset, columns, key, name))

    # TODO skip on continuous labels
    encoder = labeler.build_encoder(dataset,name,module)
    with open("data/{}/labels/{}_encoder.pkl".format(dataset, name),'wb') as pickler:
        pickle.dump(encoder, pickler, protocol=pickle.HIGHEST_PROTOCOL)

def build_custom_label(dataset, name, columns, path, key):
    print("Building custom labels named {} for dataset {} for columns {} in {} \
    on key={}".format(
        name, dataset, columns, path, key))

    data = labeler.load_metadata(path)
    headers = labeler.get_valid_columns(columns)
    headers = np.concatenate([headers, [key]])

    data = data[headers]
    data = data.set_index(key)

    data.to_pickle("data/{}/labels/{}.df".format(dataset, name))

    # TODO skip on continuous labels
    encoder = labeler.build_encoder(dataset,name,'none')
    with open("data/{}/labels/{}_encoder.pkl".format(dataset, name),'wb') as pickler:
        pickle.dump(encoder, pickler, protocol=pickle.HIGHEST_PROTOCOL)
