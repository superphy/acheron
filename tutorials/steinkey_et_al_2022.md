## The following are all of the commands used to generated the data used in steinkey et al, 2021

## Setup Acheron
The only requirement prior to these commands is installing conda (either anaconda or miniconda). The instructions to do this can be found at: https://docs.conda.io/projects/conda/en/latest/user-guide/install/

### `git clone https://github.com/superphy/acheron.git`
This downloads all code used in this project.

### `cd acheron`
You will need to be at the top of the acheron directory for this to work, for example, your terminal can look like: rylan@rylans_computer:~/Desktop/acheron$
(this doesnt need to be on the desktop, but acheron does need to be at the end)

If this step fails, try running `conda update conda` first.
If it still doesnt work, use the following docker image: continuumio/anaconda3:5.0.1
The entire project is tested to run on this image with every commit.

### `conda env create -f data/envi.yaml`
This will install all dependancies and packages this projected depends on into a conda environment.

### `conda activate ach.0.8`
This activates the conda environment we just created, the name of the environment changes and 'ach.0.8' may need to be replaced with the name at the top of data/envi.yaml

### `pip install -e .`
Dont forget the period at the end of this command, it is required to install acheron into the conda environment.

### `pytest tests/` *OPTIONAL*
This is optional, but if you want to make sure everything is working you can run the tests, takes 0.5-2 minutes.

### `mkdir -p data/salm_amr/wgs/raw`
Create a folder for your project, including a folder for genomes sequences to go into. In this case, I have named the dataset salm_amr and will contain publically available data from both the National Center for Biotechnology Information (NCBI) and the Pathosystems Resource Integration Center (PATRIC).

## Download Data

### `acheron download antibiogram -db NCBI -db PATRIC --pathogen Salmonella --email rylanj.steinkey@gmail.com --path data/PATRIC_and_NCBI_merged_abx.csv`
First we need the antimicrobial resistance data. Once we have a list of genomes that we have AMR data for, we can then download those genomes after. This downloads all AMR data from the NCBI and the PATRIC repositories and formats them into a single datasheet.

### `acheron download genomes --pathogen Salmonella -db NCBI -db PATRIC --output data/salm_amr/wgs/raw`
Download sequence data from the NCBI and PATRIC repositories. It is important the genomes go into a data/{name}/wgs/raw directory. You will be prompted to approve the download, in this case 5488 sequences from the NCBI and 907 sequences from the PATRIC. Sequences available through either database will be downloaded from whichever is passed first, in this case the NCBI. This command actually stops with an error, because some sequencespython acheron/helpers/add_resources.py dont have assembly data available (60 will fail to download). This leaves us with 6335 downloaded sequences.

## Setup Labels

Supported extensions for label files are .csv, .tsv, .xlsx, and .pkl/.df (pandas pickle).

### `acheron build label -m MIC --dataset salm_amr --name AMR_MIC --key BioSample -p data/PATRIC_and_NCBI_merged_abx.csv --columns MIC_AMP,MIC_AMC,MIC_AZM,MIC_CHL,MIC_CIP,MIC_CRO,MIC_FIS,MIC_FOX,MIC_GEN,MIC_NAL,MIC_SXT,MIC_TET,MIC_TIO,MIC_STR,MIC_KAN`

This converts the antimicrobial resistance data into machine learning ready labels. The minimum inhibitory concentration module (-m MIC) bins equivalent MIC values into the same label. For example, == 32 mg/L, > 32 mg/L, >= 32 mg/L are all encoded as the same label. Without this module, these mic's would all count as seperate labels. If the MIC panel range is from 1 mg/L up to 32 mg/L in doubling dilutions (1,2,4,8,16,32) invalid MIC's such as > 8 will be filtered out.

### `acheron build label -m MIC --dataset grdi --name AMR_MIC --key SANumber -p data/GRDI_antibiogram.xlsx --columns MIC_AMC,MIC_AMP,MIC_AZM,MIC_CHL,MIC_CIP,MIC_CRO,MIC_FOX,MIC_GEN,MIC_KAN,MIC_NAL,MIC_STR,MIC_SXT,MIC_TET,MIC_TIO`

Same as above, but converting MIC's collected as part of the Canadian Genomics Research and Development Initiative (GRDI) into machine learning ready format. Note there are difference in which antimicrobials are available.

## Setup Features

Features are being setup on supercomputer nodes with 144 core CPU's and 1TB of RAM each. For analyses when k-mers are longer than 11 nucleotides long, the work load is distributed over multiple nodes to comply with RAM limits. Lowering the MAX_GENOMES global variable in acheron/workflows/over_11mer.smk will cause greater parallelization (will use more nodes or take longer if more nodes arent available, but will use less RAM per node).

Normal computers can be used, but these tests will run in sequence instead of in parallel, taking much longer.

### If using a Slurm Cluster (supercomputer)
### `acheron build feature -c 144 -t kmer -k 11 -d salm_amr --cluster slurm`
### `acheron build feature -c 16 -t kmer -k 31 -d salm_amr --cluster slurm`

### `acheron build feature -c 144 -t kmer -k 11 -d grdi --cluster slurm`

### If using a Workstation/Desktop Computer
### `acheron build feature -c 16 -t kmer -k 11 -d salm_amr --cluster slurm`
### `acheron build feature -c 16 -t kmer -k 31 -d salm_amr --cluster slurm`

The above commands will turn our whole-genome sequence data into machine learning ready 11-mer and 31-mer count matrices. If using a slurm controller, please use the first 2 commands. If you are using a normal computer, use the second 2 commands, but know that running a 31-mer matrix count can exceed your RAM capacity very fast (about 50MB per genome).

If getting RLIMIT_NPROC or pthread_create errors with slurm, run `export OPENBLAS_NUM_THREADS=1` and `export OMP_NUM_THREADS=16`

## Train Models

### Slurm
### `acheron build model -x salm_amr -y grdi -l AMR_MIC -f 1000 -m XGB -c 8 -a AMP -t 11mer --cluster slurm`

### Desktop Computer
### `acheron build model -x salm_amr -y grdi -l AMR_MIC -f 1000 -m XGB -c 8 -a AMP -t 11mer`

These are examples of building a model on publically available data and then using it to predict on the canadian GRDI dataset. This is done to test that both datasets are valid.


### `for drug in AMP AMC AZM CHL CIP CRO FIS FOX GEN NAL SXT TET TIO STR KAN; do acheron build model -x salm_amr -l AMR_MIC -f 1000 -m XGB -c 8 -a $drug -t 11mer --cluster slurm --trial 10; done`

This command will build 5-fold cross-validated xgboost models on the NCBI & PATRIC dataset, 10 times. If you had already ran 1 iteration, this command will only run the next 9. If you already have 10, and want another ten, then you need to pass `--trial 20`


### `for drug in AMP AMC AZM CHL CIP CRO FOX GEN NAL SXT TET TIO STR KAN; do acheron build model -x grdi -l AMR_MIC -f 1000 -m XGB -c 8 -a $drug -t 11mer --cluster slurm --trial 10; done`

Same as above, but for the grdi dataset, not the grdi dataset does not have sulfisoxazole data


### `for drug in AMP AMC AZM CHL CIP CRO FOX GEN NAL SXT TET TIO STR KAN; do acheron build model -x salm_amr -y grdi -l AMR_MIC -f 1000 -m XGB -c 8 -a $drug -t 11mer --cluster slurm --trial 10; done`

Now we build models trained on public NCBI & PATRIC data, then use the trained models to predict how well they preform predicting on canadian grdi data.


### `for drug in AMP AMC AZM CHL CIP CRO FOX GEN NAL SXT TET TIO STR KAN; do acheron build model -x grdi -y salm_amr -l AMR_MIC -f 1000 -m XGB -c 8 -a $drug -t 11mer --cluster slurm --trial 10; done`

We can then do the opposite, training on grdi and using those models to predict resistance to NCBI and PATRIC data.

### python acheron/helpers/rerun_missed.py

This will run all tests as seen in the steinkey_et_al_2021 paper.

### python acheron/helpers/add_resources.py

If running using slurm, you can pull the wall time and the max RAM from sacct. This command will add Elapsed time and max ram to the summary file in the results folder, alongside accuracy, various error rates. Its important to note that ram is sampled in intervals, not continously (which would have preformance implications). So if there was a temporary spike in RAM that did not exceed the allocated memory of that specific job, MaxRSS might not catch that spike. Also, if a job failed and then was finished in a 2nd job, that test will have 2 sacct files. i.e. trials 1-4 might have a different jobid as trials 5-10. Therefore, if generating average RAM use, trials 1-4 will have to be divided by 4, and trials 5-10 will have to be divided by 6.

###

Now that we have finished models, we can see how well they did.

There are predefined summaries that can be generated, see the wiki.

###

Now that you have finished models, as well as identified top k-mers, we can annotate the genomes and get an idea of where these k-mers are and what they might do.
