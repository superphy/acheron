import joblib

from acheron.workflows import sequence_tools

if config['model'].upper() in ['XGB','XGBOOST']:
    model_ext = '.bst'
else:
    raise Exception("Model type {} not supported for feature extraction".format(config['model']))

model_path = "model={}_train={}_test={}_validate={}_\
feats={}_hyp={}_cvfolds={}_attribute={}_trial={}".format(
    config['model'],config['train'],config['test'],config['validation'],
    config['feats'],config['hyperparam'],config['cv'],config['attribute'],
    config['trial'])

rule all:
    input:
        model_path+'.blast'

        
rule make_blast_db:
    output:
        "data/{}/wgs/master.db".format(config['train'])
    run:
        sequence_tools.make_master_fasta(config['train'])
        master_path = "data/{}/wgs/master.fasta".format(config['train'])

        shell("makeblastdb -in {} -parse_seqids -dbtype nucl -out {}".format(
        master_path, output[0]))

# Find top n features based on their importance to the model
# and find their location in the genomes
rule blast_top_feats:
    input:
        model_path+model_ext,
        "data/{}/wgs/master.db".format(config['train'])
    output:
        model_path+'.blast'
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

"""
rule find_hits:
    input:
        "data/.....mer_blast_hits/words.blast",
        "gffpandas/{id}.pkl"
    output:
        "hits.pkl"
    run:
        from sequence_tools import find_hits

rule hit_summary:
    input:
        "hits.pkl"
    output:
        "mer_summary.csv"
    run:
        from sequence_tools import hit_summary

rule score_summary:
    input:
        "mer_summary.csv"
    output:
        "score_summary.csv"

rule select_best_hits:

rule add_card_hits:
"""
