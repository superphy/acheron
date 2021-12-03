import pandas as pd
import numpy as np
import joblib
import os
import glob

from itertools import chain
from acheron.predict import find_genomes

path = config['path']
module = config['module']
out = config['out']
cores = config['cores']
pathogen = config['pathogen']

if os.path.isfile(path):
    if path[-6:] != '.fasta':
        raise Exception("{} is not a path to a fasta file".format(path))
    genome_paths = [path]
elif os.path.isdir(path):
    genome_paths = find_genomes(path)
else:
    raise Exception("Only files or directories allowed, {} was unexpected".format(path))

if len(genome_paths) <= 0:
    raise Exception("No valid fastas found at {}".format(path))

fasta_ids = [i.split('/')[-1].split('.')[0] for i in genome_paths]
genome_dir = '/'.join(genome_paths[0].split('/')[:-1])

def get_top_feats(module,pathogen):
    if module.upper() in ['MIC','ABX']:
        master_feats = []
        if pathogen.upper() == 'SALMONELLA':
            drugs = ['AMP','AMC','AZM','CHL','CIP','CRO','FOX','GEN','NAL','SXT','TET','TIO','STR','KAN','FIS']
        elif pathogen.upper() == 'CAMPYLOBACTER':
            drugs = ['AMP','AZM','CIP','GEN','NAL','TET','CLI','ERY','FLO','TEL']
        else:
            raise Exception("pathogen {} not defined in predictor.smk".format(pathogen))

        for drug in drugs:
            bst = joblib.load("data/predict/models/{}/{}/{}.bst".format(module.upper(),pathogen.lower(),drug))
            master_feats.append(bst.feature_names)

        top = list(chain(*master_feats))
        return top
    else:
        raise Exception("Module {} not supported, use `acheron help` for more information".format(module))

rule all:
    input:
        expand("data/predict/jellyfish_results/"+pathogen+"/{id}.fa", id=fasta_ids)

rule get_top:
    output:
        "data/predict/models/MIC/{}/top.fasta".format(pathogen)
    run:
        top_feats = get_top_feats(module,pathogen)
        np.save(output[0],top_feats)
        with open("data/predict/models/MIC/{}/top.fasta".format(pathogen),'a') as fh:
            for feat in top_feats:
                fh.write(">{}\n".format(feat))
                fh.write(feat+"\n")

rule clean:
    input:
        genome_dir+"/{id}.fasta"
    output:
        "data/predict/cleaned_genomes/"+pathogen+"/{id}.fasta"
    shell:
        "python acheron/workflows/clean.py {input} data/predict/cleaned_genomes/"+pathogen

rule count:
    input:
        "data/predict/cleaned_genomes/"+pathogen+"/{id}.fasta"
    output:
        temp("data/predict/jellyfish_results/"+pathogen+"/{id}.jf")
    threads:
        2
    shell:
        "jellyfish count -C -m 11 -s 100M -t {threads} {input} -o {output}"

rule dump:
    input:
        "data/predict/jellyfish_results/"+pathogen+"/{id}.jf",
        "data/predict/models/MIC/{}/top.fasta".format(pathogen)
    output:
        "data/predict/jellyfish_results/"+pathogen+"/{id}.fa"
    shell:
        "jellyfish query {input[0]} -s data/predict/models/MIC/"+pathogen+"/top.fasta > {output}"
