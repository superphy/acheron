import os, sys

import pandas as pd
import numpy as np

import click

from acheron.workflows import sequence_downloader
from acheron.workflows import sequence_tools

databases = config['databases']
databases = databases.split('_')

output_path = config['out']
pathogen = config['pathogen']

data = {}
for database in databases:
        data[database] = pd.read_csv("data/{}_{}_antibiogram.csv".format(database,pathogen))

seq_sources = sequence_downloader.find_seq_sources(data)
dl_lists = sequence_downloader.download_sources(seq_sources)

# List of biosamples to be downloaded
biosamples = list(dl_lists['PATRIC'])
biosamples = [i for i in biosamples if len(str(i)) > 1]

# THIS SECTION DEFAULTS YES AND SKIPS FOR SLURM COMPATABILITY
try:
    cont = click.confirm("You are about to download {} sequences from the PATRIC, would you like to continue?".format(len(biosamples)),
        default = True)
except:
    cont = True
if not cont:
    print("Exiting download at user request")
    sys.exit()

# PATRIC wants genomd id's instead of biosample numbers. We need to convert
ids = sequence_tools.biosamples_to_genomeids(biosamples, pathogen)
bio2id = dict(zip(biosamples,ids))

rule all:
    input:
        expand(output_path+'/'+"{bio}.fasta", bio=biosamples)

rule download_and_rename:
    output:
        output_path+'/'+"{bio}.fasta"
    run:
        id = bio2id[wildcards.bio]
        dl_path = "ftp://ftp.patricbrc.org/genomes/{0}/{0}.fna".format(id)
        shell("wget {dl_path} -O {output}")
"""
rule download:
    output:
        output_path+'/'+"{id}.fna"
    shell:
        "wget -P {output_path} ftp://ftp.patricbrc.org/genomes/{wildcard.id}/{wildcard.id}.fna"

rule rename:
    input:
        output_path+'/'+bio2id[wildcards.bio]+"{}.fna"
    output:
        output_path+'/'+"{bio}.fasta"
    shell:
        "mv {input} {output}"
"""
