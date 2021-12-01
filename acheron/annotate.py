import os

def annotate_genomes(dataset, cores):
    print("Annotating genomes in data/{}/wgs/raw".format(dataset))

    smk = "acheron/workflows/annotate_genomes.smk"

    os.system("snakemake -s {0} -j {1} --config dataset={2} cores={1}".format(
    smk, cores, dataset))

def identify_important_regions(arguments):
    # First call the annotate snakemake, to make sure the genomes are annotated
    dataset = getattr(arguments,'train')
    cores = getattr(arguments,'cores')

    annotate_genomes(dataset, cores)

    # Prepare arguments for indentifier snakemake
    model_args = ''
    for arg in vars(arguments):
        model_args += " {}={}".format(arg, getattr(arguments,arg))

    # Call to region identify
    os.system("snakemake -s acheron/workflows/identify_regions.smk -j {} \
    --config{}".format(cores, model_args))
