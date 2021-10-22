import pandas as pd
import sys

def load_steinkey2021():
    from acheron.helpers.rerun_missed import steinkey2021
    tests = steinkey2021()

    print("Looking for {} results".format(len(tests)))

    """ tests in form
    model,train,test,validate,feats,hyp,cvfolds,attribute,trial
    ['XGB', 'salm_amr', 'none', 'none', '100', 'False', 5, 'AMP', '0']
    """
    cols = ['model','train','test','validate','feats','hyp','cvfolds','attribute','trial']

    res = []
    for test in tests:
        for i in range(10):
            res.append(test[:-1]+[str(i)])

    results = pd.DataFrame(data = res, columns = cols)
    results['type'] = pd.Series(['11mer' for i in res])

    return results

def add_results(tests):
    """
    Takes a df of tests, such as from load_steinkey2021,
    then adds to the dataframe the results of those tests
    """
    # res order for MIC is
    cols = ['Supports', 'Accuracy', 'Within 1 Dilution', 'Very Major Error',
    'Major Error', 'Non Major Error', 'Correct', 'SLURM_JOBID', 'Time (m)',
    'Max Ram (GiB)']

    stats = []
    for test in tests.values:
        res_path = "results/model={}_train={}_test={}_validate={}_feats={}_type=11mer_hyp={}_cvfolds={}_attribute={}_trial={}/summary.df".format(*test)
        res = pd.read_pickle(res_path)

        try:
            assert list(res.columns) == cols
        except:
            print("The following file is missing columns")
            print(res_path)
            print(res.columns)
            raise
        stats.append(res.values)

    for i, col in enumerate(cols):
        tests[col] = pd.Series(stats[:,i])

    return tests

def make_summary(subset, out, media):
    if subset == 'steinkey2021':
        summary = load_steinkey2021()
        summary = add_results(summary)
    else:
        raise Exception("Summaries need to be defined, either make your own or call one in summary.py")

    if media == 'table':
            if out == 'stdout':
                print(summary)
            elif out[-5:] == '.xlsx':
                summary.to_excel(out)
            elif out[-4:] == '.csv':
                summary.to_csv(out)
            elif out[-4:] == '.tsv':
                summary.to_csv(out,sep='\t')
            elif out[-4:] == '.pkl':
                summary.to_pickle(out)
            else:
                raise exception("File extension {} not supported".format(out[-4:]))

    elif media = 'figure':
        pass

    else:
        raise Exception("Summary media not in ['table','figure']")
