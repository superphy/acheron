
dataset = config['dataset']
seq_path = "data/{}/wgs/raw/".format(dataset)

ids, = glob_wildcards(seq_path+"{id}.fasta")

print('testing only 4, delete this')
ids = ids[:4]

rule all:
    input:
        expand("data/"+dataset+"/wgs/annotations/{id}/{id}.pkl", id = ids)

# Clean samples if you havent already run other acheron modules
rule clean:
    input:
        "data/"+dataset+"/wgs/raw/{id}.fasta"
    output:
        "data/"+dataset+"/wgs/clean/{id}.fasta"
    run:
        from clean import get_files_to_analyze, format_files
        all_files = get_files_to_analyze(input[0])
        format_files(all_files, "data/"+dataset+"/wgs/clean/")

# Identify which genes are in the genome, and where they are located
rule prokka:
    input:
        "data/"+dataset+"/wgs/clean/{id}.fasta"
    output:
        "data/"+dataset+"/wgs/annotations/{id}/{id}.ffn"
    threads:
        1
    shell:
        "prokka {input} --outdir data/"+dataset+"/wgs/annotations/{wildcards.id} --prefix {wildcards.id} --cpus {threads} --force --compliant"

# stored results in a pandas dataframe for easier access
rule prokka_to_dataframe:
    input:
        "data/"+dataset+"/wgs/annotations/{id}/{id}.ffn"
    output:
        "data/"+dataset+"/wgs/annotations/{id}/{id}.pkl"
    run:
        shell("export OPENBLAS_NUM_THREADS=1")
        import pandas as pd
        import gffpandas.gffpandas as gffpd
        anno = gffpd.read_gff3("data/"+dataset+"/wgs/annotations/{0}/{0}.gff".format(wildcards.id))
        anno.df.to_pickle(str(output))
