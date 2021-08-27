import os
import glob

def find_genomes(path):
    """
    Takes a path to a directory of genomes,
    Returns a list of genome paths
    """
    paths = []
    for fasta in glob.glob("{}/*.fasta".format(path)):
        paths.append(fasta)
    return paths

def make_predictions(path,module,out,cores,cluster):
    """
    Path can be either a directory, or a single file
    """
    if module.upper() not in ["MIC", "ABX"]:
        raise Exception("Currently only supporting prediction of abx MIC values, not {}".format(module))

    if os.path.isfile(path):
        if path[-6:] != '.fasta':
            raise Exception("{} is not a path to a fasta file".format(path))
    elif os.path.isdir(path):
        pass
    else:
        raise Exception("Only files or directories allowed, {} was unexpected".format(path))

    RAM = (cores*2)+2

    """
    The purpose of the snakemake is to count the kmers, the rest happens below
    """
    if cluster.upper() != "NONE":
        raise Exception("Cluster support not yet ready, please wrap acheron yourself")

    #if cluster.upper()=='NONE':
    os.system("snakemake -s acheron/workflows/predictor.smk -j {3} \
    --config path={0} module={1} out={2} cores={3}".format(path,module,out,cores))
    #elif cluster.upper()== "SLURM":
    #    os.system("sbatch -c {3} acheron/workflows/predictor.smk --mem {4}G snakemake -s acheron/workflows/predictor.smk -j {3} \
    #    --config path={0} module={1} out={2} cores={3}".format(path,module,out,cores,RAM))

    # at this point, we can assume that all the jellyfish files have been created
    # the matrix has not been created as we want it in memory and not saved

    print("ending prediction")

    # remaining steps
    # generate k-mer matrix
    # predict on each row of k-mer matrix
    # offer accuracy testing if fed a truth table?
