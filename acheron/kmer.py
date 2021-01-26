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

def build_kmer_matrix(dataset, kmer_length, cores, cluster):
    """
    SLURM NODE INFORMATION
    If running out of memory, lower cores (stay in multiples of 16), keep RAM high

    For example,
    k = 11 : -c 144 --mem 1000G
    k = 31 : -c 16  --mem 1000G
    """
    SLURM_RAM = 1000

    print("Building {}-mer matrix of {} sequences".format(kmer_length, dataset))
    if kmer_length <= 11:
        workflow_smk = "acheron/workflows/sub_11mer.smk"
    else:
        workflow_smk = "acheron/workflows/over_11mer.smk"

    if cluster.upper()=='NONE':
        os.system("snakemake -s {0} -j {1} \
        --config kmer_length={2} dataset={3} cores={1}".format(
        workflow_smk, cores, kmer_length, dataset))
    elif cluster.upper()== "SLURM":
        os.system("sbatch -c {1} --mem {4}G snakemake -s {0} -j {1} \
        --config kmer_length={2} dataset={3} cores={1}".format(
        workflow_smk, cores, kmer_length, dataset, SLURM_RAM))
    else:
        raise Exception("cluster option {} not recognized, only slurm or none allowed".format(cluster))
