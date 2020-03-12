import numpy as np
import pandas as pd

from acheron.workflows import labeler

def build_module_label(dataset, module, name, columns, path, key):
    # call to snakemake subworkflow in workflows
    print("Building labels named {} for dataset {} for columns {} in {} on \
    key={} based on module {}".format(
        name, dataset, columns, path, key, module))

    """
        os.system("snakemake -s {0} -j {1} \
        --config kmer_length={2} dataset={3} cores={1}".format(
        workflow_smk, cores, kmer_length, dataset))
    """

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
