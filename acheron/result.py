import pandas as pd
pd.set_option("display.max_columns",500)
pd.set_option("display.width", 1000)

def get_result(arguments):
    """
    Reads the corresponding results file,
    prints out the summary dataframe
    """
    first = "results/model={}_train={}_test={}_validate={}_feats={}_type={}_hyp={}".format(
    getattr(arguments,'model'),
    getattr(arguments,'train'),
    getattr(arguments,'test'),
    getattr(arguments,'validation'),
    getattr(arguments,'num_features'),
    getattr(arguments,'type'),
    getattr(arguments,'hyperparam'))


    second = "_cvfolds={}_attribute={}_trial={}/summary.df".format(
    getattr(arguments,'cv'),
    getattr(arguments,'attribute'),
    getattr(arguments,'trial'))

    df = pd.read_pickle(first+second)

    return df


def print_results(arguments):
    print(get_result(arguments))
