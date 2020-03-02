import os
import glob

def get_genomes(dataset):
    print("looking in data/{}/wgs/raw/*.fasta".format(dataset))
    genomes = []
    for fasta in glob.glob("data/{}/wgs/raw/*.fasta".format(dataset)):
        genomes.append(fasta)
    return genomes

def print_genome_names(genomes):
    for genome in genomes:
        print(genome)

def build_kmer_matrix(dataset, kmer_length, cores):
    print("kmer.py:")
    print("Building {}-mer matrix of {} sequences".format(kmer_length, dataset))
    os.system("snakemake -s acheron/workflows/sub_11mer.smk -j {0} --config kmer_length={1} dataset={2} cores={0}".format(
    cores, kmer_length, dataset))
    #os.system("snakemake -s acheron/workflows/sub_11mer.smk --config kmer_length={}".format(kmer_length))
