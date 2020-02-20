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

def build_kmer_matrix(dataset, kmer_length):
    print("Building {}-mer matrix of {} sequences".format(kmer_length, dataset))
