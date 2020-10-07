import pandas as pd
import numpy as np
import joblib
import glob

from acheron.workflows import supervised_model

# set to 5 for 5 fold cross-validation
k_folds = config['folds']

# if we are cross validating we only need one feature dataframe, otherwise
# we need to check if the other datasets are also ready
features = [] # paths to feature matrices
labels = [] # paths to label matrices
masks = [] # paths to validity masks

# if a dataset has been declared, add it to above paths
for dataset in ['train','test','validation']:
    if config[dataset] not in ['None','none']:
        features += ["data/{}/features/{}_matrix.df".format(
            config[dataset],config['type'])]
        labels += ["data/{}/labels/{}.df".format(
            config[dataset],config['label'])]
        masks += ["data/{}/features/masks/{}_{}.df".format(
            config[dataset],config['type'],config['label'])]

splits = [] # path to splits df, for when only train is passed

if config['test'] == 'none' and config['validation'] == 'none':
    for trial in range(config['trial']):
        splits += ["data/{}/masks/split{}_{}_{}_{}xCV.df".format(
            config['train'],trial,config['type'],config['label'],config['cv'])]

print("this needs to be set to build an array of attributes, if multiple or 'all' is passed")
columns = [config['attribute']]
trials = [i for i in range(config['trial'])]

rule all:
    input:
        expand("results/model={}_train={}_test={}_validate={}_feats={}_hyp={}_cvfolds={}".format(
        config['model'], config['train'], config['test'], config['validation'],
        config['num_features'], config['hyperparam'], config['cv'])+'_attribute={atb}_trial={trl}.df',atb=columns, trl=trials)

rule mask:
    input:
        labels
    output:
        masks
    run:
        for dataset_num in range(len(labels)):
            label = pd.read_pickle(labels[dataset_num])

            mask = supervised_model.make_mask(label, len(splits))
            mask.to_pickle(output[dataset_num])

# splits will only run when cross validating (i.e. test == None)
rule split_samples:
    input:
        masks
    output:
        splits
    run:
        samples = glob.glob("data/{}/wgs/raw/*.fasta".format(config['train']))

        # keep in mind, we only have splits when we are cross validating
        mask = pd.read_pickle(masks[0])
        for trial in trials:
            label_matrix = pd.read_pickle("data/{}/labels/{}.df".format(config['train'],config['label']))
            split = supervised_model.make_split(label_matrix, mask, config['cv'], samples)
            split.to_pickle("data/{}/masks/split{}_{}_{}_{}xCV.df".format(
                config['train'],trial,config['type'],config['label'],config['cv']))

rule build_model:
    input:
        features+splits+masks
    output:
        "results/model={}_train={}_test={}_validate={}_feats={}_hyp={}_cvfolds={}".format(
        config['model'], config['train'], config['test'], config['validation'],
        config['num_features'], config['hyperparam'], config['cv'])+'_attribute={atb}_trial={trl}.df'
    threads:
        8
    run:
        model, results = supervised_model.make_model(config['model'],config['train'],
            config['test'],config['validation'],config['label'],config['type'],
            config['attribute'],config['num_features'],config['hyperparam'],config['cv'],wildcards.trl)
        joblib.dump(model, output[0].split['.'][0]+'.bst')
        results.to_pickle(output[0])
