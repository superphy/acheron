# Acheron was built with Salmonella in mind, as such, there are slight modifications to be done for other bacteria

# The following is an example of using Acheron for analyses on a bacteria not already described, in this case, campylobacter

## Overview

Step 0: Setup Acheron (see other tutorial for setup guide)
Step 1: Download AMR data
Step 2: Download sequence data
Step 2.5: Add your own data, should you have any
Step 3: Declare MIC ranges
Step 4: Convert AMR data in ML ready labels
Step 5: Convert sequence data into ML ready features
Step 6: Train models

## Step 1: Download AMR data
sbatch -c 1 --mem 4G --partition NMLResearch --wrap='acheron download antibiogram -db NCBI -db PATRIC --pathogen Campylobacter --email rylanj.steinkey@gmail.com --path data/public_campy_amr.csv'

## Step 2: Download sequence data
sbatch -c 1 --mem 4G --partition NMLResearch --wrap='acheron download genomes --pathogen Campylobacter -db NCBI --output data/public_campy/wgs/raw'

## Step 2.5: Add your own data, should you have any
Throw your sequences in data/pick_name/wgs/raw
You will have to do steps 4 and 5 twice, once for each dataset

## Step 3: Declare MIC ranges
Class ranges were defined in data/mic/label_modules/mic/Campylobacter_class_ranges.yaml
Change these if you want

## Step 4: Convert AMR data in ML ready labels
sbatch -c 1 --mem 4G --partition NMLResearch --wrap='acheron build label -m MIC --dataset public_campy --name AMR_MIC --key BioSample -p data/public_campy_amr.csv --columns MIC_AMP,MIC_AZM,MIC_CIP,MIC_GEN,MIC_NAL,MIC_TET,MIC_CLI,MIC_ERY,MIC_FLO,MIC_TEL --pathogen=Campylobacter'

## Step 5: Convert sequence data into ML ready features
acheron build feature -c 144 -t kmer -k 11 -d public_campy --cluster slurm


## Step 6: Train models
acheron build model -x public_campy -l AMR_MIC -f 1000 -m XGB -c 8 -a AMP -t 11mer --cluster slurm

see tutorials/steinkey_et_al_2021.md for examples on building models and training on multliple datasets
