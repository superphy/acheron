import joblib
import os,sys
import pandas as pd
import numpy as np
import skbio.io

from acheron.workflows import sequence_tools

if config['model'].upper() in ['XGB','XGBOOST']:
    model_ext = '.bst' #TODO add in model name
else:
    raise Exception("Model type {} not supported for feature extraction".format(config['model']))

model_path = "results/model={}_train={}_test={}_validate={}_\
feats={}_type={}_hyp={}_cvfolds={}_attribute={}_trial={}/".format(
    config['model'],config['train'],config['test'],config['validation'],
    config['feats'],config['type'],config['hyperparam'],config['cv'],
    config['attribute'],config['trial'])

rule all:
    input:
        # TODO


rule make_blast_db:
    output:
        "data/{}/wgs/master.db".format(config['train'])
    run:
        sequence_tools.make_master_fasta(config['train'])
        master_path = "data/{}/wgs/master.fasta".format(config['train'])

        shell("makeblastdb -in {} -parse_seqids -dbtype nucl -out {}".format(
        master_path, output[0]))

# Find top n features based on their importance to the model
rule blast_top_feats:
    input:
        model_path+model_ext,
        "data/{}/wgs/master.db".format(config['train'])
    output:
        model_path+'top_feats.blast'
    run:
        # load model
        model = joblib.load(input[0])

        # pull top n feats and their importances
        importances = sequence_tools.find_top_feats(
            model, config['model'],config['num_top'])

        # search for locations of kmers

        query_path = output[0][:-5]+"query"
        sequence_tools.save_query(importances, query_path)

        blast_command = "blastn -task blastn-short -db data/{}/wgs/blast_database.db\
        -query {} -ungapped -perc identity 100 -dust no -word_size {}\
        -max_target_seqs 50000 -evalue 100000 -outfmt 6 -out {}".format(
        config[train], query_path, kmer_size, output[0])

        shell(blast_command)

# Takes the top features as determined above and finds their locations
# Also returns nearby genes and their distances to the k-mer
rule find_hits:
    input:
        model_path+'top_feats.blast'
    output:
        "annotations/"+config['train']+"_{drug}_"+config['type']+"/hits_df.pkl" # TODO Fix path
    run:
        from acheron import gene_finder

        if config['train'] = 'grdi' and wildcards.drug == 'FIS':
            # This drug doesnt exist in the GRDI dataset, so we skip
            shell("touch {output}")
            sys.exit(0)

        #top_feats =  TODO read in top feats as npy 1d arr

        hits = gene_finder.find_hits(top_feats,input)

        hits.to_pickle(output)


rule hit_summary:
    input:
        expand("annotations/"+config['train']+"_{drug}_"+config['type']+"/hits_df.pkl",drug=drugs) #TODO fix path
    output:
        "annotation/{}_hit_summary.csv".format(config['type'])
    run:
        from acheron import gene_finder
        gene_finder.hit_summary(config['train'],output)

rule score_summary:
    input:
        "annotation/{}_hit_summary.csv".format(config['type'])
    output:
        "annotation/{}_score_summary.csv".format(config['type'])

rule select_best_hits:

rule add_card_hits:
