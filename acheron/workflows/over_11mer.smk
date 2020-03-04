from clean import get_files_to_analyze
import math

# Location of the raw genomes
RAW_GENOMES_PATH = "data/{}/wgs/raw/".format(config['dataset'])

# Kmer length that you want to count
KMER_SIZE = config['kmer_length']

# Data type of the resulting kmer matrix. Use uint8 if counts are
# all under 256. Else use uint16 (kmer counts under 65536)
# note: counts should be under 256 anyways, so if thats not the case,
# check your sequences
MATRIX_DTYPE = 'uint8'

# it is assumed there is enough ram for the union tests to handle MAX_GENOMES
# This is used supercluster nodes with 1TB RAM, lower substantially
# if using a normal computer
MAX_GENOMES = 1000
GENOMES = get_files_to_analyze(RAW_GENOMES_PATH)
NUM_GENOMES = len(GENOMES)
UNION_SPLITS = [i+1 for i in range(math.ceil(NUM_GENOMES/MAX_GENOMES))]

ids, = glob_wildcards(RAW_GENOMES_PATH+"{id}.fasta")

# from sub_11mer.smk, we will use clean_fastas, count_kmers, dump_kmers
subworkflow kmer:
    snakefile:
        "sub_11mer.smk"

rule all:
    input:
        "data/{}/features/{}mer_matrix.df".format(config['dataset'],
        config['kmer_length'])"

# for subsets of genomes, adding this as a seperate step saves on ram
# as we cannot fit all genomes in memory at a single time, so we batch the kmers
# into groups of at most MAX_GENOMES, union those, and then union the unions

# builds 2D numpy array, each row is a batch with len(colums)<=MAX_GENOMES
rule batch_genomes_for_union:
    input:
        kmer(expand("data/"+config['dataset']+"/wgs/"+
        str(config['kmer_length'])+"mer_jellyfish_results/{id}.fa",id=ids))
    output:
        "data/{}/wgs/master_{}mers/batches.npy".format(config['dataset'],
        config['kmer_length'])
    run:
        from acheron.workflows import find_common_kmers
        batches = find_common_kmers.batch_kmers(genomes, len(UNION_SPLITS))
        np.save(output[0], batches)

rule batch_kmers:
    input:
        "data/{}/wgs/master_{}mers/batches.npy".format(config['dataset'],
        config['kmer_length'])
    output:
        "data/"+config['dataset']+"/wgs/master_"+config['kmer_length']+
        "mers/kmer_subset_{set_num}_of_"+str(len(UNION_SPLITS))+'.npy'
    run:
        from acheron.workflows import find_common_kmers
        union_of_batch = find_common_kmers.find_union(config['kmer_length'],
        config['dataset', batch_row, config['cores']])
        np.save(output[0], union_of_batch)

rule merge_kmer_batches:
    input:
        expand("data/"+config['dataset']+"/wgs/master_"+config['kmer_length']+
        "mers/kmer_subset_{set_num}_of_"+str(len(UNION_SPLITS))+'.npy',
        set_nums = UNION_SPLITS)
    output:
        "data/"+config['dataset']+"/wgs/master_"+config['kmer_length']+
        "mers/all_kmers.npy'
    run:
        from acheron.workflows import find_common_kmers
        # im hoping here that input is a list of inputs that i can use,
        # gonna have to test this one
        master_mers = find_common_kmers.merge_masters(config['dataset'],
        config['kmer_length'], input)
        np.save(output[0], master_mers)


# now that we know what kmers we will eventually see, we can make parts
# of the dataframe
rule build_sub_matrices:
    input:
        "data/"+config['dataset']+"/wgs/master_"+config['kmer_length']+
        "mers/all_kmers.npy'
    output:
        #each split
    run:
        #multi_mer_matrix
rule merge_sub_matrices:
    input:
        #each split
    output:
        #master_df
    run:
        #multi_mer_merge
