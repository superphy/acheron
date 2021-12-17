# Acheron

[![CircleCI](https://circleci.com/gh/superphy/acheron.svg?style=svg)](https://circleci.com/gh/superphy/acheron)

Acheron is a command line based workflow manager that allows for large scale, fast, and memory efficient machine-learning analyses of whole-genome sequence data.

This repository is under active development, see wiki for instructions.

## Acheron lets you:

### Download

- Download antibiogram data from different databases, check for inconsistencies, and save all the data in a common format.
- Download whole-genome sequence data from different databases, downloading only missing genomes you might want to include in your analyses.

### Convert to machine learning ready formats
- Turn AMR data into ready to use, ML friendly labels, including binning of MIC values into predefined ranges.
- Turn whole-genome sequence data into k-mer count matrices ready for ML. Analyses using k-mer's longer than 11-mer's are automatically batched into subgroups to comply with user-defined memory restrictions.
- Frequency counts of 11-mer's up to 256 for 6,000 genomes can be stored for 12 GB with no loss in information.
- Frequency counts of 31-mer's up to 256 for 6,000 genomes can be stored for 320 GB with no loss in information and computed using high performance parallel computing clusters with less than 1TB of RAM.

### Train machine learning models
- Choose between sci-kit learn, keras, or XGBoost models.
- Enable nested cross validation for hyperparameter optimizations.
- Train 11-mer models in 10-15 minutes for 6,000 genomes on consumer grade computers.
- Train 31-mer models in 1-2 hours for 6,000 genomes on high performance parallel computing clusters (w/o GPU's, w/ 1TB RAM).

### Annotate and Identify
- Annotate your genomes to determine which genes are where.
- Extract the most important k-mers from your models and identify regions of the genome important in your predictions.
- Identify if the genes with the most importance contain known antimicrobial properties.


## Future Plans:
- Train models based on virulence data (labels).
- Train models based on area under the curve (AUC) omnilog 96-well plates (features).
- Train models based on genes instead of k-mers, for those with extremely limited compute ability.
