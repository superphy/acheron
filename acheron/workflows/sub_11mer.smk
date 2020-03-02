# Location of the raw genomes
RAW_GENOMES_PATH = "data/{}/wgs/raw/".format(config['dataset'])

# Kmer length that you want to count
KMER_SIZE = config['kmer_length']

# Data type of the resulting kmer matrix. Use uint8 if counts are
# all under 256. Else use uint16 (kmer counts under 65536)
# note: counts should be under 256 anyways, so if thats not the case,
# check your sequences
MATRIX_DTYPE = 'uint8'

ids, = glob_wildcards(RAW_GENOMES_PATH+"{id}.fasta")

rule all:
    input:
        "data/{}/features/{}mer_matrix.df".format(config['dataset'], config['kmer_length'])

rule clean:
    input:
        RAW_GENOMES_PATH+"{id}.fasta"
    output:
        "data/"+config['dataset']+"/wgs/clean/{id}.fasta"
    run:
        from clean import get_files_to_analyze, format_files
        all_files = get_files_to_analyze(input[0])
        format_files(all_files, "data/"+config['dataset']+"/wgs/clean/")

rule count:
    input:
        "data/"+config['dataset']+"/wgs/clean/{id}.fasta"
    output:
        temp("data/"+config['dataset']+"/wgs/jellyfish_results/{id}.jf")
    threads:
        2
    shell:
        "jellyfish count -C -m {KMER_SIZE} -s 100M -t {threads} {input} -o {output}"

rule dump:
    input:
        "data/"+config['dataset']+"/wgs/jellyfish_results/{id}.jf"
    output:
        "data/"+config['dataset']+"/wgs/jellyfish_results/{id}.fa"
    shell:
        "jellyfish dump {input} > {output}"

rule sub11_matrix:
    input:
        expand("data/"+config['dataset']+"/wgs/jellyfish_results/{id}.fa", id=ids)
    output:
        "data/{}/features/{}mer_matrix.df".format(config['dataset'], config['kmer_length'])
    run:
        from sub_11mer_matrix import make_matrix
        df = make_matrix(config['kmer_length'], MATRIX_DTYPE, config['cores'], "data/"+config['dataset']+"/wgs/jellyfish_results/", "data/"+config['dataset']+"/features/" )
        df.to_pickle(output[0])
