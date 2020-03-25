import os

def build_model(arguments):

    model_args = ''
    for arg in vars(arguments):
        model_args += " {}={}".format(arg, getattr(arguments,arg))

    os.system("snakemake -s acheron/workflows/modeler.smk -j {} \
    --config{}".format(getattr(arguments,'cores'), model_args))
