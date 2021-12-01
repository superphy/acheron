import pandas as pd
import numpy as np
import sys
from collections import Counter

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
        """
        for i in range(10):
            res.append(test[:-1]+[str(i)])
        """
        res.append(test)

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
    failed = []
    for test in tests.values:
        res_path = "results/model={}_train={}_test={}_validate={}_feats={}_type=11mer_hyp={}_cvfolds={}_attribute={}_trial={}/summary.df".format(*test)
        res = pd.read_pickle(res_path)


        if 'Time (m)' not in res.columns:
            failed.append(res_path)
            continue

        for i, ele in enumerate(res.columns):
            assert ele == cols[i]

        stats.append(res.values[0])

    if len(failed) > 0:
        import warnings
        warnings.warn("{} files did not have valid Time and Ram information, continuing anyways".format(len(failed)))

    stats = np.array(stats)
    for i, col in enumerate(cols):
        tests[col] = pd.Series(stats[:,i])

    return tests

def weigh_jobs_per_slurm(summary):
    """
    Scales the ram and time requirements according to how many jobs ran in that
    slurm instance. For example. If a slurm instance runs 6 jobs, the time for that
    instance is divided by 6.
    """
    num_jobs = dict(Counter(summary['SLURM_JOBID']))
    job_times = zip(summary['SLURM_JOBID'],summary['Time (m)'])
    new_times = [i[1]/num_jobs[i[0]] for i in job_times]
    summary['Time (m)'] = pd.Series(new_times)

    return summary

def convert(summary):
    """
    Takes a df and converts each column to the relevant datatype.
    e.g. feats is cast from a string to a float
    """

    summary[['feats','cvfolds']] = summary[['feats','cvfolds']].apply(pd.to_numeric)
    return summary

def get_stats(model_df):
    n = np.sum(model_df['Supports'])

    acc = np.mean(model_df['Accuracy'])
    oned_acc = np.mean(model_df['Within 1 Dilution'])
    vme = np.mean(model_df['Very Major Error'])
    me = np.mean(model_df['Major Error'])
    nme = np.mean(model_df['Non Major Error'])
    time = np.mean(model_df['Time (m)'])
    ram = np.mean(model_df['Max Ram (GiB)'])

    acc_stdev = np.std(model_df['Accuracy'])
    oned_acc_stdev = np.std(model_df['Within 1 Dilution'])
    vme_stdev = np.std(model_df['Very Major Error'])
    me_stdev = np.std(model_df['Major Error'])
    nme_stdev = np.std(model_df['Non Major Error'])
    time_stdev = np.std(model_df['Time (m)'])
    ram_stdev = np.std(model_df['Max Ram (GiB)'])

    return n,acc,oned_acc,vme,me,nme,time,ram,acc_stdev,oned_acc_stdev,vme_stdev,me_stdev,nme_stdev,time_stdev,ram_stdev

def summarize(results, subset, media):
    """
    Turns the table of thousands of results into simple, reportable stats
    i.e. 95.2% (+/- 2.45%)
    cols = ['model','train','test','validate','feats','hyp','cvfolds','attribute','trial']
    """

    if subset =='steinkey2021':
        datasets = [['salm_amr','none','none'],['grdi','none','none'],['salm_amr','grdi','none'],['grdi','salm_amr','none']]
        models = ['XGB','SVM','ANN']

        # If doing a summary of all drugs in a model, then set attribute to 'all_diverse'
        abx = ['AMP','AMC','AZM','CHL','CIP','CRO','FIS','FOX','GEN','NAL','SXT','TET','TIO','STR','KAN']
        types = ['11mer','31mer'] # then manually set feats based on type

        sum_cols=['model','train','test','validate','feats','type','hyp','cvfolds','attribute']
        res_cols=['n', 'Accuracy', 'Accuracy stdev', 'Within 1 Dilution', 'Within 1 Dilution stdev',
        'Very Major Error','Very Major Error  stdev','Major Error','Major Error stdev',
        'Non Major Error', 'Non Major Error stdev','Time (m)', 'Time (m) stdev',
        'Max Ram (GiB)','Max Ram (GiB) stdev']
        """
        res_cols=['Supports', 'Accuracy', 'Within 1 Dilution', 'Very Major Error',
        'Major Error', 'Non Major Error', 'Correct', 'SLURM_JOBID', 'Time (m)',
        'Max Ram (GiB)']
        """
        summs = []

        if media == 'table':
            for ds in datasets:
                # TODO hyperparam'd models
                dataset_df = results[(results['train']==ds[0]) & (results['test']==ds[1]) & (results['validate']==ds[2]) & (results['hyp']=='False')]
                for type in types:
                    if type=='31mer':
                        feats = '1000000'
                    else:
                        feats = '1000'
                    type_df = dataset_df[(dataset_df['type']==type) & (dataset_df['feats']==feats)]
                    for model in models:
                        model_df = type_df[type_df['model']==model]
                        # as of this point, there are 150 samples in the df, n=10 for each of 15 abx

                        stats = get_stats(model_df)
                        summs.append([model,ds[0],ds[1],ds[2],feats,type,False,5,'all',*stats])
                        print("{} predicted with 1D accuracy of {}% (+/- {}) using {} {}s xyv:[{},{},{}]".format(model,stats[1],stats[8],feats,type,ds[0],ds[1],ds[2]))


                        for drug in abx:
                            # These are for accuracies specific to each abx, rather than overall
                            drug_df = model_df[model_df['attribute']==drug]
                            drug_stats = get_stats(drug_df)
                            summs.append([model,ds[0],ds[1],ds[2],feats,type,False,5,drug,*drug_stats])
        elif media == 'figure':
            for ds in datasets:
                # TODO hyperparam'd models
                dataset_df = results[(results['train']==ds[0]) & (results['test']==ds[1]) & (results['validate']==ds[2]) & (results['hyp']=='False')]
                for type in types:
                    type_df = dataset_df[(dataset_df['type']==type)]
                    for feat in set(type_df['feats']):
                        feat_df = type_df[type_df['feats']==feat]
                        for model in models:
                            model_df = feat_df[feat_df['model']==model]

                            stats = get_stats(model_df)
                            summs.append([model,ds[0],ds[1],ds[2],feats,type,False,5,'all',*stats])

                            for drug in abx:
                                # These are for accuracies specific to each abx, rather than overall
                                drug_df = model_df[model_df['attribute']==drug]
                                drug_stats = get_stats(drug_df)
                                summs.append([model,ds[0],ds[1],ds[2],feats,type,False,5,drug,*drug_stats])
        else:
            raise Exception("media {} not defined in summarize".format(media))
        summs_df = pd.DataFrame(data=summs, columns=sum_cols+res_cols)
        return summs_df
    else:
        raise Exception("Summaries need to be defined, either make your own or call one in summary.py")

def make_summary(subset, out, media):
    if subset == 'steinkey2021':
        results = load_steinkey2021()
        results = add_results(results)
        results = weigh_jobs_per_slurm(results)
        results = convert(results)
    else:
        raise Exception("Summaries need to be defined, either make your own or call one in summary.py")

    if media == 'table':
        summary = summarize(results, subset, media)
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

    elif media == 'figure':
        from acheron import figures
        for i in range(5):
            if i == 3:
                figures.group_figures(subset, results, out, i)

    else:
        raise Exception("Summary media not in ['table','figure']")
