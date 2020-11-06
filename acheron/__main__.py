import sys
import argparse

from multiprocessing import cpu_count

# Build Functions
from .kmer import build_kmer_matrix
from .genes import build_genes_matrix
from .abricate import build_abricate_matrix
from .omnilog import build_omnilog_matrix

from .label import build_module_label
from .label import build_custom_label

from .download import download_antibiogram
from .download import download_genomes

from .model import build_model

from .annotate import annotate_genomes
from .annotate import identify_important_regions

# Result Functions
# TODO

# Summary Caller
# TODO

def main():
    arguments = parse_arguments()

    # Build Caller
    if arguments.action_command == 'build':
        if arguments.build_command == 'feature':
            if arguments.type == 'kmer':
                build_kmer_matrix(arguments.dataset, arguments.kmer_length, arguments.cores)
            elif arguments.type == 'genes':
                build_genes_matrix(arguments.dataset)
            elif arguments.type == 'abricate':
                build_abricate_matrix(arguments.dataset, arguments.database)
            elif arguments.type == 'omnilog':
                build_omnilog_matrix(arguments.dataset)
            else:
                raise argparse.ArgumentError(arguments.type, "unexpected build --type: {}".format(argument.type))

        elif arguments.build_command == 'label':
            if arguments.module != 'custom':
                build_module_label(
                    arguments.dataset, arguments.module, arguments.name,
                    arguments.columns, arguments.path, arguments.key)
            else:
                build_custom_label(
                    arguments.dataset, arguments.name, arguments.columns,
                    arguments.path, arguments.key)

        elif arguments.build_command == 'model':
            build_model(arguments)

        else:
            raise argparse.ArgumentError(arguments.build_command, "acheron build requires another positional argument from the list above")

    # Annotation Caller
    elif arguments.action_command == 'annotate':
        annotate_genomes(arguments.dataset, arguments.cores)

    elif arguments.action_command == 'identify':
        identify_important_regions(arguments)

    # Result Caller
    elif arguments.action_command == 'result':
        print('finish me')

    # Summary Caller
    elif arguments.action_command == 'summary':
        print('finish me ')

    # Download Caller
    elif arguments.action_command == 'download':
        if arguments.download_command == 'antibiogram':
            download_antibiogram(arguments.database, arguments.pathogen, arguments.email, arguments.antimicrobial, arguments.path, arguments.use_local, arguments.check_date)
        elif arguments.download_command == 'genomes':
            download_genomes(arguments.input, arguments.output)
        else:
            raise argparse.ArgumentError(arguments.download_command, "Only downloading genomes and antibiogram data is supported right now")

    else:
        raise argparse.ArgumentError(arguments.action_command, "acheron requires one of the positional arguments listed above")

def parse_arguments():
    # parser to be inherited
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--dry_run', default=False, action='store_true',
                    help="Checks if all requirements of the given command are satisfied")
    parent_parser.add_argument('-c', '--cores', default = max(1,cpu_count()-1),
                    help="Number of cores to use, defaults to number of system cores minus 1")
    parent_parser.add_argument('-o', '--out', default = 'stdout',
                    help="Output path/file to save results")

    test_params = argparse.ArgumentParser(add_help=False)
    test_params.add_argument('-x', '--train', required=True,
                    help="Name of dataset to tests models on")
    test_params.add_argument('-y', '--test', default='none',
                    help="Name of dataset to tests models on, not passing this results sets cv=True")
    test_params.add_argument('-v', '--validation', default='none',
                    help="name of dataset to validate hyperparameters with")
    test_params.add_argument('-f', '--num_features',
                    help="Number of features to keep past feature selection, not passing will skip feature selection")
    # labels required for supervised only
    test_params.add_argument('-l', '--label', required=True,
                    help="what labels to base the models on, created by `acheron build label ...` ")
    test_params.add_argument('--columns', default='none',help="subset of columns in label")
    test_params.add_argument('-m', '--model', default='XGB', choices = ['XGB','SVM','ANN','kSNP3'],
                    help="The model you would like to build")
    test_params.add_argument('-p', '--hyperparam', default=False,
                    help="Enable hyperparameter optimizations and nest the cross validations, will use training set to validate hyperparams")
    test_params.add_argument('-a', '--attribute',required=True,
                    help="Which attribute to train the model on (column of label)")
    test_params.add_argument('-t', '--type', required=True,
                    help="which features the model is based on (i.e. 11mer or AMR)")
    test_params.add_argument('--trial', default=1,
                    help="to run the same test multiple times, change trial number")
    test_params.add_argument('--cv', default=5, help="number of folds in cross validation")

    download_params = argparse.ArgumentParser(add_help=False)
    download_params.add_argument('-db', '--database', required=True, action='append',
                    help="Which database or databases in [NCBI, PATRIC] to download from. Add multiple with `acheron download -db NCBI -db PATRIC`")
    download_params.add_argument('--pathogen', required=True, help="Which pathogen for which to download the data")

    # main parser
    root_parser = argparse.ArgumentParser()

    action_subparsers = root_parser.add_subparsers(title='action', dest='action_command',
                    help="What action you would like to take in ['build','result','annotate','identify','summary','download']")

    # Download subparser
    download_parser = action_subparsers.add_parser('download',
                    help="For downloading metadata and loading it into the correct format")
    download_subparsers = download_parser.add_subparsers(title='download', dest="download_command")

    antibiogram_parser = download_subparsers.add_parser('antibiogram',
                    help="For downloading/updating antibiogram data", parents=[download_params,parent_parser])
    antibiogram_parser.add_argument('--email', help="NCBI requires an email address", required=True)
    antibiogram_parser.add_argument('-abx', '--antimicrobial', default='all')
    antibiogram_parser.add_argument('-p','--path', help="Save path and file name for resulting dataframe", required=True)
    antibiogram_parser.add_argument('--use_local', help="Same usage as -db, will use local copy for merging.\
                    For example, if you wanted to use local PATRIC but pull a new NCBI, you would pass `--use_local PATRIC`",
                    action='append')
    antibiogram_parser.add_argument('--check_date', default=False, action='store_true',
                    help="checks the date the antibiograms were pulled")
    # add a use_local path

    genome_parser = download_subparsers.add_parser('genomes', help="For downloading missing genomes" )
    genome_parser.add_argument('--input', help="Path to antibiogram sheet containing biosamples", required=True)
    genome_parser.add_argument('--output', help="Location to save to and check for pre-existing genomes", required=True)

    # Build Subparser
    build_parser = action_subparsers.add_parser('build',
                    help="For building feature matrices, data labels, or models")
    build_subparsers = build_parser.add_subparsers(title='build', dest='build_command')

    feature_parser = build_subparsers.add_parser('feature',
                    help="For building feature matrices", parents=[parent_parser])
    feature_parser.add_argument('-t', '--type', required = True, choices = ['kmer','genes','abricate','omnilog'])
    feature_parser.add_argument('-k', '--kmer_length', type = int, default=11,
                    help="Length of kmer to use, note k > 11 takes substantial resources, see docs")
    feature_parser.add_argument('-db', '--database', choices = ['AMR','VF'], default='VF',
                    help="Choose between building AMR or VF with abricate")
    feature_parser.add_argument('-d', '--dataset', required = True,
                    help="Name of dataset, what the name of the folder containing sequences is named")


    label_parser = build_subparsers.add_parser('label',
                    help="For building labels for the data",
                    parents = [parent_parser])
    label_parser.add_argument('-m', '--module', default = 'custom', choices = ['MIC'],
                    help="Specify the pre-built module you would like to build from")
    label_parser.add_argument('-n', '--name', required=True,
                    help="name of the labels you are creating, will be used later")
    label_parser.add_argument('--columns', required=True,
                    help="comma seperated listed of columns headers, or path to numpy list")
    label_parser.add_argument('--key', required=True,
                    help="columns header containing the id/name/sequence filename")
    label_parser.add_argument('-p', '--path', required=True,
                    help="path to xlsx, csv, or tsv containing label data")
    label_parser.add_argument('-d', '--dataset', required = True,
                    help="Name of dataset, what the name of the folder containing sequences is named")


    model_parser = build_subparsers.add_parser('model', parents=[parent_parser,test_params],
                    help="For building machine learning models")
    model_parser.add_argument('-k', '--folds', default = 5,
                    help="How many folds for k fold cross validation, default=5")

    # Annotate subparser
    # takes same arguments as model builder, analyses results instead of builds model
    annotate_parser = action_subparsers.add_parser('annotate', parents=[parent_parser],
                    help="annotates whole-genome sequences")
    annotate_parser.add_argument('-d', '--dataset', required = True,
                    help="Name of dataset, what the name of the folder containing sequences is named")

    # Identify subparser
    # takes same arguments as model builder, analyses results instead of builds model
    identify_parser = action_subparsers.add_parser('identify', parents=[parent_parser,test_params],
                    help="Identifies important regions as determined by a trained model")
    identify_parser.add_argument("--num_top", default=5,
                    help="How many of the top features to search for")

    # Result Subparser
    # need to specify each parameter, for a single model
    result_parser = action_subparsers.add_parser('result', parents=[parent_parser,test_params],
                    help="Prints the results of a single test, i.e. xgboost model at 1000 features")

    # Summary Subparser
    # give range of results to summarize, and figure
    summary_parser = action_subparsers.add_parser('summary', parents=[parent_parser,test_params],
                    help="For general results and figures, i.e. compare 3 model types at 50 different feature sizes")
    # TODO: sub parser for acheron summary table, acheron summary figure

    args = root_parser.parse_args()

    if len(sys.argv) < 2:
        root_parser.print_help(sys.stderr)
    elif len(sys.argv) < 3:
        if args.action_command == 'build':
            build_parser.print_help(sys.stderr)
        elif args.action_command == 'result':
            result_parser.print_help(sys.stderr)
        elif args.action_command == 'annotate':
            annotate_parser.print_help(sys.stderr)
        elif args.action_command == 'identify':
            identify_parser.print_help(sys.stderr)
        elif args.action_command == 'summary':
            summary_parser.print_help(sys.stderr)
        elif args.action_command == 'download':
            download_parser.print_help(sys.stderr)
        else:
            raise Exception('uncaught erroneous action_command')
    return args

if __name__ == "__main__":
    main()
