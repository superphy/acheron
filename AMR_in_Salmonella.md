The following commands make up the entire of the AMR prediction for Salmonella pipeline

# first download metadata
[IN]$ acheron download antibiogram -db NCBI -db PATRIC --pathogen Salmonella --email "rylanj.steinkey@gmail.com" --path "data/PATRIC_and_NCBI_merged_abx.csv"

[OUT]
Looking in 2 database(s): ['NCBI', 'PATRIC']
for Salmonella antibiogram data using the email rylanj.steinkey@gmail.com

Number of samples found in the NCBI database:  5490
100% [........................................................................] 50965564 / 50965564Skipping the following antimicrobials because they were not requested (and how many times they are seen):
['levofloxacin (15)', 'ceftazidime (6)', 'tigecycline (13)', 'amikacin (2)', 'ampicillin/sulbactam (3)', 'aztreonam (2)', 'cefazolin (2)', 'cefepime (3)', 'cefotaxime (2)', 'colistin (2)', 'doripenem (2)', 'ertapenem (2)', 'imipenem (3)', 'meropenem (3)', 'piperacillin/tazobactam (2)', 'polymyxin B (2)', 'tobramycin (2)', 'nitrofurantoin (1)']

89077 duplicate but matching mic values were found
84 MIC values were ignored because there were multiple conflicting values, these were:
590.16607: ['amoxicillin/clavulanic acid', 'amoxicillin/clavulanic acid', 'ampicillin', 'ampicillin', 'azithromycin', 'azithromycin', 'cefoxitin', 'cefoxitin', 'ciprofloxacin', 'ciprofloxacin']
590.16602: ['ampicillin', 'ampicillin', 'azithromycin', 'azithromycin', 'ceftiofur', 'ceftiofur', 'nalidixic acid', 'nalidixic acid', 'sulfisoxazole', 'sulfisoxazole']
590.16826: ['gentamicin', 'gentamicin', 'nalidixic acid', 'nalidixic acid']
590.17199: ['azithromycin', 'azithromycin', 'chloramphenicol', 'chloramphenicol', 'ciprofloxacin', 'ciprofloxacin']
590.15977: ['chloramphenicol', 'chloramphenicol', 'gentamicin', 'gentamicin', 'nalidixic acid', 'nalidixic acid']
590.17817: ['cefoxitin', 'cefoxitin', 'ceftriaxone', 'ceftriaxone', 'chloramphenicol', 'chloramphenicol', 'gentamicin', 'gentamicin']
590.16622: ['amoxicillin/clavulanic acid', 'amoxicillin/clavulanic acid', 'cefoxitin', 'cefoxitin', 'ceftiofur', 'ceftiofur', 'gentamicin', 'gentamicin', 'nalidixic acid', 'nalidixic acid', 'sulfisoxazole', 'sulfisoxazole']
590.15701: ['cefoxitin', 'cefoxitin']
590.16559: ['azithromycin', 'azithromycin', 'gentamicin', 'gentamicin', 'nalidixic acid', 'nalidixic acid']
590.16986: ['cefoxitin', 'cefoxitin', 'ceftriaxone', 'ceftriaxone', 'gentamicin', 'gentamicin', 'nalidixic acid', 'nalidixic acid']
590.16092: ['azithromycin', 'azithromycin', 'chloramphenicol', 'chloramphenicol', 'gentamicin', 'gentamicin', 'nalidixic acid', 'nalidixic acid']
590.15668: ['cefoxitin', 'cefoxitin', 'chloramphenicol', 'chloramphenicol']
Automatically merging the following columns: ['MIC_AMC', 'MIC_AMP', 'MIC_AZM', 'MIC_FOX', 'MIC_TIO', 'MIC_CRO', 'MIC_CHL', 'MIC_CIP', 'MIC_GEN', 'MIC_NAL', 'MIC_STR', 'MIC_FIS', 'MIC_TET', 'MIC_SXT', 'MIC_KAN', 'isolation_source', 'serovar', 'collection_date', 'strain']

The following columns had conflicting data between columns (and the number of removed values):
{'MIC_AMP': 44, 'MIC_AZM': 59, 'MIC_FOX': 84, 'MIC_TIO': 41, 'MIC_CRO': 28, 'MIC_CHL': 64, 'MIC_CIP': 36, 'MIC_GEN': 84, 'MIC_NAL': 51, 'MIC_STR': 72, 'MIC_FIS': 81, 'MIC_TET': 18}

# next move in / download genomes
[IN]$ ls data/ncbi_salm/wgs/raw/

[OUT]
SAMN02367860.fasta  SAMN02640777.fasta  SAMN02640780.fasta  SAMN02640784.fasta  SAMN02640792.fasta  SAMN02640799.fasta  SAMN02640803.fasta  SAMN02640806.fasta
SAMN02367879.fasta  SAMN02640778.fasta  SAMN02640781.fasta  SAMN02640788.fasta  SAMN02640797.fasta  SAMN02640800.fasta  SAMN02640804.fasta  SAMN02640807.fasta
SAMN02367906.fasta  SAMN02640779.fasta  SAMN02640783.fasta  SAMN02640789.fasta  SAMN02640798.fasta  SAMN02640801.fasta  SAMN02640805.fasta  SAMN02640809.fasta

# make 11-mer and/or 31-mer matrices

[IN]$ acheron build feature -t kmer -k 11 -d ncbi_salm

[OUT]
Building 11-mer matrix of ncbi_salm sequences
Provided cores: 7
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	all
	1	build_sub11mer_matrix
	24	clean_fastas
	24	count_kmers
	24	dump_kmers
	74

.... follow by tons of snakemake progress output ....

Finished job 21.
72 of 74 steps (97%) done

rule build_sub11mer_matrix:
    input: data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640806.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640798.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640789.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640799.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640805.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640781.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640779.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640792.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640809.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640807.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640783.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640803.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02367879.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640797.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640780.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640778.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640784.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02367860.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640800.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640777.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640804.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640801.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02367906.fa, data/ncbi_salm/wgs/11mer_jellyfish_results/SAMN02640788.fa
    output: data/ncbi_salm/features/11mer_matrix.df
    jobid: 1

    Finished job 1.
    73 of 74 steps (99%) done

    localrule all:
        input: data/ncbi_salm/features/11mer_matrix.df
        jobid: 0

    Finished job 0.
    74 of 74 steps (100%) done

# now build the label matrices

[IN]$ acheron build label --module MIC --name PATRIC_and_NCBI --columns all --key BioSample -p data/PATRIC_and_NCBI_merged_abx.csv -d ncbi_salm

[OUT]
Building labels named PATRIC_and_NCBI for dataset ncbi_salm for columns all in data/PATRIC_and_NCBI_merged_abx.csv on     key=BioSample based on module MIC
Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	all
	1	make_mic_df
	2

rule make_mic_df:
    input: data/PATRIC_and_NCBI_merged_abx.csv, data/label_modules/mic/class_ranges.yaml
    output: data/ncbi_salm/labels/PATRIC_and_NCBI.df, data/label_modules/mic/mic_class_order_dict.pkl
    log: logs/ncbi_salm/PATRIC_and_NCBI_mic.log
    jobid: 1

... MIC panel info printouts ....

Finished job 1.
1 of 2 steps (50%) done

localrule all:
    input: data/ncbi_salm/labels/PATRIC_and_NCBI.df
    jobid: 0

Finished job 0.
2 of 2 steps (100%) done

# build model
