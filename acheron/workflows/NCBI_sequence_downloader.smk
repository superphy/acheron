import os, sys

import pandas as pd
import numpy as np

import click

from acheron.workflows import sequence_downloader

databases = config['databases']
databases = databases.split('_')

output_path = config['out']
pathogen = config['pathogen']

data = {}
for database in databases:
        data[database] = pd.read_csv("data/{}_{}_antibiogram.csv".format(database,pathogen))

seq_sources = sequence_downloader.find_seq_sources(data)
dl_lists = sequence_downloader.download_sources(seq_sources)

#List of ids to be downloaded
try:
    ids = list(dl_lists['NCBI'])
except:
    print("No samples to be downloaded from NCBI")
    sys.exit()
# THIS SECTION DEFAULTS YES AND SKIPS FOR SLURM COMPATABILITY
try:
    cont = click.confirm("You are about to download {} sequences from the NCBI, would you like to continue?".format(len(ids)),
        default = True)
except:
    cont = True
if not cont:
    print("Exiting download at user request")
    sys.exit()
"""
user_input = ''
while user_input.upper() not in ['Y','N','YES','NO','YE']:
    print("You are about to download {} sequences from the NCBI, would you like to continue? (y/n)".format(len(ids)))
    user_input = input()

if user_input.upper() in ['NO','N']:
    sys.exit()
"""

# Path where you want the assemblies to be saved
ASMBL = output_path+'/'

rule all:
     input:
        expand(ASMBL+"{id}.fasta", id=ids)

# Get input file list
rule ncbi_fna_list:
    output:
        temp(ASMBL+"{id}_esearch.txt")
    shell:
        """
        esearch -db biosample -query {wildcards.id} \
          | elink -target assembly \
          | esummary \
          | grep "FtpPath_GenBank" \
          | sed -r "s|.+>(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCA_.+)<.+|\\1\\2/\\2_genomic.fna.gz|" \
          > {output}
        """

# Download the files from the list
# unzip into .fasta
rule ncbi_dl:
    input:
        ASMBL+"{id}_esearch.txt"
    output:
        comp = temp(ASMBL+"{id}/{id}.fna.gz")
    shell:
        "wget -O {output} -i {input}"


rule unzip:
    input:
        ASMBL+"{id}/{id}.fna.gz"
    output:
        ASMBL+"{id}/{id}.fna"
    shell:
        "gunzip {input}"

rule move:
    input:
        ASMBL+"{id}/{id}.fna"
        #dir = ASMBL+"{id}"
    output:
        ASMBL+"{id}.fasta"
    shell:
        """
        mv {input} {output}
        """
# rm -r {input.dir}
