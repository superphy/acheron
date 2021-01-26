import os
MAX_RAM = '1000G'
def build_model(arguments):
    model_args = ''
    for arg in vars(arguments):
        attr = getattr(arguments,arg)
        model_args += " {}={}".format(arg, attr)

        if arg == 'hyperparam':
            # hyperparameter optimizations require additional RAM
            # subject to change
            if attr == True:
                RAM = '250G'
            else:
                RAM = '125G'

    # if using >11mers, max out the ram
    if getattr(arguments,'type')[-3:] == 'mer':
        if int(getattr(arguments,'type')[:2]) > 11:
            RAM = MAX_RAM


    cluster = getattr(arguments,'cluster').upper()

    # when not using cluster
    if cluster == 'NONE':
        os.system("snakemake -s acheron/workflows/modeler.smk -j {} \
        --config{}".format(getattr(arguments,'cores'), model_args))

    # when using slurm cluster management
    elif cluster == 'SLURM':
        os.system("sbatch -c {0} --mem {1} snakemake -s acheron/workflows/modeler.smk -j {0} \
        --config{2}".format(getattr(arguments,'cores'), RAM, model_args))

    else:
        raise Exception("cluster config {} not supported, use slurm or none".format(cluster))
