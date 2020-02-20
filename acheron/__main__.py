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

from .model import build_model

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
                build_kmer_matrix(arguments.dataset, arguments.kmer_length)
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

    # Result Caller
    elif arguments.action_command == 'result':
        print('finish me')

    # Summary Caller
    elif arguments.action_command == 'summary':
        print('finish me ')

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
    parent_parser.add_argument('-d', '--dataset', required = True,
                    help="Name of dataset, what the name of the folder containing sequences is named")

    test_params = argparse.ArgumentParser(add_help=False)
    test_params.add_argument('-x', '--train', required=True,
                    help="Name of dataset to tests models on")
    test_params.add_argument('-y', '--test', required=False,
                    help="Name of dataset to tests models on, not passing this results sets cv=True")
    test_params.add_argument('-f', '--num_features',
                    help="Number of features to keep past feature selection, not passing will skip feature selection")
    # labels required for supervised only
    test_params.add_argument('-l', '--label',
                    help="what labels to base the models on, created by `acheron build label ...` ")
    test_params.add_argument('-m', '--model', default='XGB', choices = ['XGB','SVM','ANN','kSNP3'],
                    help="The model you would like to build")
    test_params.add_argument('-p', '--hyperparam',
                    help="Enable hyperparameter optimizations and nest the cross validations, will use training set to validate hyperparams")

    # main parser
    root_parser = argparse.ArgumentParser()

    action_subparsers = root_parser.add_subparsers(title='action', dest='action_command',
                    help="What action you would like to take in ['build','result','summary']")

    # Build Subparser
    build_parser = action_subparsers.add_parser('build',
                    help="For building feature matrices, data labels, or models")
    build_subparsers = build_parser.add_subparsers(title='build', dest='build_command')
    #build_subparsers.required=True

    feature_parser = build_subparsers.add_parser('feature',
                    help="For building feature matrices", parents=[parent_parser])
    feature_parser.add_argument('-t', '--type', required = True, choices = ['kmer','genes','abricate','omnilog'])
    feature_parser.add_argument('-k', '--kmer_length', type = int, default=11,
                    help="Length of kmer to use, note k > 11 takes substantial resources, see docs")
    feature_parser.add_argument('-db', '--database', choices = ['AMR','VF'],
                    help="Choose between building AMR or VF with abricate")

    label_parser = build_subparsers.add_parser('label',
                    help="For building labels for the data")
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

    model_parser = build_subparsers.add_parser('model', parents=[parent_parser,test_params],
                    help="For building machine learning models")

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
        elif args.action_command == 'summary':
            summary_parser.print_help(sys.stderr)
        else:
            raise Exception('uncaught erroneous action_command')
    return args

if __name__ == "__main__":
    main()
