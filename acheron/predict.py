import os, sys
import glob
import joblib
import pickle

import numpy as np
import pandas as pd

from Bio import Seq, SeqIO
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat, chain

from acheron.workflows.supervised_model import predict, select_features

def find_genomes(path):
    """
    Takes a path to a directory of genomes,
    Returns a list of genome paths
    """
    paths = []
    for fasta in glob.glob("{}/*.fasta".format(path)):
        paths.append(fasta)
    return paths

def make_row_from_query(jf_path, relevant_feats):
    """
    Given a genome file, create and return a row of kmer counts
    to be inserted into the kmer matrix.
    """

    cols_dict = { relevant_feats[i] : i for i in range(0, len(relevant_feats))}

    # Create a temp row to fill and return (later placed in the kmer_matrix)
    temp_row = [0]*len(relevant_feats)

    # place kmer counts in correct order
    with open(jf_path) as file:
        for line in file:
            line = line.rstrip()
            kmer, count = line.split()
            if(int(count)>255):
                count = 255
            temp_row[cols_dict[kmer]] = int(count)

    id = jf_path.split('/')[-1].split('.')[0]

    return id, temp_row

def build_matrix(relevant_feats, jf_list, cores):
    """
    Builds the feature matrix, e.g. k-mer counts
    """
    genomes = [i.split('/')[-1] for i in jf_list]
    runs = [i.split('.')[0] for i in genomes]
    total = len(runs)

    # declaring empty kmer matrix to fill
    kmer_matrix = np.zeros((total,len(relevant_feats)),dtype = 'uint8')

    rows_dict = { runs[i] : i for i in range(0, len(runs))}
    cols_dict = { relevant_feats[i] : i for i in range(0, len(relevant_feats))}

    with ProcessPoolExecutor(max_workers=cores) as ppe:
        for genome_name,temp_row in ppe.map(make_row_from_query, jf_list, repeat(relevant_feats)):
            kmer_matrix[rows_dict[genome_name],:] = temp_row

    df = pd.DataFrame(data = kmer_matrix, columns = relevant_feats, index = runs)
    return df

def load_module(module):
    """
    Loads everything about module from storage into memory
    Everything past this point is volatile
    """
    if module.upper() in ["MIC", "ABX"]:
        boosters = []
        mics = []
        master_feats = []
        for drug in ['AMP','AMC','AZM','CHL','CIP','CRO','FOX','GEN','NAL','SXT','TET','TIO','STR','KAN','FIS']:
            bst = joblib.load("data/predict/models/{}/{}.bst".format(module.upper(),drug))
            boosters.append(bst)
            mics.append(drug)
            master_feats.append(bst.feature_names)
        top = list(set(chain(*master_feats)))

        return top, boosters, mics

    else:
        raise Exception("Module {} not defined".format(module))

def decoder(preds, attr, module):
    """
    Takes a list of encoded predictions, converts them back into
    human readable labels
    """
    with open("data/predict/models/{}/encoder.pkl".format(module.upper()),'rb') as fh:
        encoder = pickle.load(fh)

    decoder = {v:k for k,v in encoder[attr].items()}
    return [decoder[i] for i in preds]

def make_predictions(path,module,out,cores,cluster):
    """
    Path can be either a directory, or a single file
    """
    if module.upper() not in ["MIC", "ABX"]:
        raise Exception("Currently only supporting prediction of abx MIC values, not {}".format(module))

    if os.path.isfile(path):
        if path[-6:] != '.fasta':
            raise Exception("{} is not a path to a fasta file".format(path))
        genome_paths = [path]
    elif os.path.isdir(path):
        genome_paths = find_genomes(path)
    else:
        raise Exception("Only files or directories allowed, {} was unexpected".format(path))

    RAM = (cores*2)+2

    """
    The purpose of the snakemake is to count the kmers, the rest happens below
    """
    if cluster.upper() != "NONE":
        raise Exception("Cluster support not yet ready, please wrap acheron yourself")

    #if cluster.upper()=='NONE':
    os.system("snakemake -s acheron/workflows/predictor.smk -j {3} \
    --config path={0} module={1} out={2} cores={3}".format(path,module,out,cores))
    #elif cluster.upper()== "SLURM":
    #    os.system("sbatch -c {3} acheron/workflows/predictor.smk --mem {4}G snakemake -s acheron/workflows/predictor.smk -j {3} \
    #    --config path={0} module={1} out={2} cores={3}".format(path,module,out,cores,RAM))

    # at this point, we can assume that all the jellyfish files have been created
    # the matrix has not been created as we want it in memory and not saved

    top_feats, models, attributes = load_module(module)

    ids = [i.split('/')[-1].split('.')[0] for i in genome_paths]
    jf_list = ["data/predict/jellyfish_results/{}.fa".format(i) for i in ids]
    matrix = build_matrix(top_feats, jf_list, cores)

    rows = []
    for model, attr in zip(models,attributes):
        # only thing relevant to pass to select_features is the feature matrix and the feature_names
        features = select_features(matrix,'none',1000,model.feature_names)

        # select_features puts the features in the order required by the model
        # Error will be thrown anyways if features arent in correct order
        res = predict(model, features, 'XGB')
        labels = decoder(res, attr, module)
        rows.append(labels)

    prediction_df = pd.DataFrame(data=np.transpose(rows), columns=attributes, index=ids)

    if out == 'stdout':
        print(prediction_df)
    elif out[-5:] == '.xlsx':
        prediction_df.to_excel(out)
    elif out[-4:] == '.csv':
        prediction_df.to_csv(out)
    elif out[-4:] == '.tsv':
        prediction_df.to_csv(out,sep='\t')
    elif out[-4:] == '.pkl':
        prediction_df.to_pickle(out)
    else:
        raise exception("File extension {} not supported".format(out[-4:]))
