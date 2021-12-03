# This document details the steps taken to generate models for use in `acheron predict --pathogen campylobacter`

Note that none of this is wrapped for slurm, you need to add the appropriate flags or wrap it yourself.

## `git clone https://github.com/superphy/acheron.git`

## `cd acheron`

## `conda env create -f data/envi.yaml`
This command downloads all the packages that acheron uses.

## `conda activate ach.X.X`
This will put you into an environment that uses the above packages. For example, if you use pandas in this environement, it will run the version acheron specifies. Note that the version changes and the X's need to be replaced with the name at the top of the file named data/envi.yaml.

## `pip install -e .`
This installs acheron into the environement we just made. Don't forget the peroid at the end of the command.

## `pytest tests/` OPTIONAL
This is optional, but if you want to make sure everything is working you can run the tests, takes 0.5-2 minutes.

## `acheron download antibiogram -db NCBI -db PATRIC --pathogen Campylobacter --email rylanj.steinkey@gmail.com --path data/campy_MICs.csv`

This command searches for campylobacter antimicrobial resistance data in the NCBI and PATRIC databases. This command found 945 samples an saved the MIC, SIR, and metadata (source, assembly methods, etc) to data/data/campy_MICs.csv

## `acheron download genomes --pathogen Campylobacter -db NCBI -db PATRIC --output data/public_campy/wgs/raw`

This command looks at which samples were found using the last command, and then downloades the accompying whole genome sequence. Don't expect all of them to be found, in this case, 938 (99.3%) were found. This step might fail at the end, thats fine, it just couldn't resolved the other 7 sequences.

## `acheron build label -m MIC --dataset public_campy --name AMR_MIC --key BioSample -p data/campy_MICs.csv --columns MIC_AMP,MIC_AZM,MIC_CIP,MIC_GEN,MIC_NAL,MIC_TET,MIC_CLI,MIC_ERY,MIC_FLO,MIC_TEL --pathogen=Campylobacter`

This takes the antimicrobial resistance information we downloaded and converts it to machine learning-ready formats.
Labels can be found at data/public_campy/labels/AMR_MIC.df

## `acheron build feature -c 144 -t kmer -k 11 -d public_campy`

This turns the folder of whole-genome sequences into an 11-mer matrix. This can found at data/public_campy/features/11mer_matrix.df

## `cp -r data/public_campy/ data/campy_copy`

If we only pass one dataset to acheron, it will automatically cross-validate. When we train these prediction models, we want the entirety of the data to see seen by the model. This is not allowed as most times it means the user is making a mistake and it needs to be corrected. By copying our dataset, we can tell acheron it is different, bypassing this restriction.

## open acheron/acheron/helpers/model_evaluators.py

Add in any missing abx, the values dont matter, but it will allow acheron to run.


## `for drug in AMP AZM CIP GEN NAL TET CLI ERY FLO TEL; do acheron build model -x public_campy -y campy_copy -l AMR_MIC -f 1000 -m XGB -c 8 -a $drug -t 11mer --cluster slurm; done`

This command instructs acheron to train (without-cross validation) on the public_campy set, and then test it using the campy_copy. The actual results will be useless (hence why acheron doesnt allow this), but the model's outputted will be fully trained on all available campylobacter data.

## cp results/{...}/model0.joblib data/predict/models/MIC/campylobacter/DRUG.bst

Copy the trained model into the predictions folder. Repeat for each drug.

## cp data/public_campy/labels/AMR_MIC_encoder.pkl data/predict/models/MIC/campylobacter/encoder.pkl

Move the encoder that was made earlier into the predictions folder, this converts predictions (0,1,2,3,4) into their actually class (e.g. <=0.15 mg/L)

# All done!
