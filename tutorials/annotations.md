## This document details the steps of identifying import features

We will be starting off in the 'setup features' section seen in steinkey_et_al_2022.md

## acheron build model -x salm_amr -l AMR_MIC -t 31mer -f 1000 -a AMP --cluster slurm --trial 10

This will fail (it is out of order) because we havent made the feature matrices yet, but it will go through
and assign splits to each sample before failing, which we can then use to make the features.

This can be cancelled when `ls -1 data/salm_amr/splits/split*_31* | wc -l` returns 10

## acheron build feature -c 16 -t kmer -k 31 -d salm_amr --cluster slurm --prefiltering --trials 10

This is similar to whats seen before, except we will build the matrix on 31mers, instead of 11mers.

To do this, we use prefiltering, where features are scored for variance when they are in memory.
This means that when it comes time to load the data for training, we can load selectively based on
the recorded variance, instead of having to do it at runtime. This wouldn't be possible on ~500 GB matrices
without prefiltering.

## acheron annotate -d ncbi_salm -c 8

todo
