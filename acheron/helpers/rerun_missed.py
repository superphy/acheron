"""
Declare which tests should be run,
tells you which are missing

This is so you dont have to rerun all the snakemakes,
as well as informs you which has failed.

"""
import pandas as pd
import sys, os
import itertools


def steinkey2022():
    """
    Returns all the tests that were run with the steinkey et al., 2022 paper
    return list is in the order [model,train,test,validate,feats,type,hyp,cvfolds,attribute,trial]
    """
    attrs = ['AMP','AMC','AZM','CHL','CIP','CRO','FIS','FOX','GEN','NAL','SXT','TET','TIO','STR','KAN']
    mdls = ['SVM'] # only put non XGB-ANN models here
    grdi_models = ['SVM','XGB','ANN']
    #currently leaving out the 1mil feat set
    full_feats = ['100','1000','10000','100000']+[str(i) for i in range(10,100,10)]
    part_feats = ['100','1000','10000']
    trials = [str(i) for i in range(10)]

    tests = []
    for i in attrs:
        for l in trials:
            for k in full_feats:
                tests.append(['XGB','salm_amr','none','none',k,'11mer','False',5,i,l])
            for k in part_feats:
                tests.append(['XGB','salm_amr','none','none',k,'11mer','True',5,i,l])
                tests.append(['ANN','salm_amr','none','none',k,'11mer','False',5,i,l])
                tests.append(['ANN','salm_amr','none','none',k,'11mer','True',5,i,l])
            for j in mdls:
                tests.append([j,'salm_amr','none','none','1000','11mer','False',5,i,l])
            # grdi
            if i != 'FIS':
                for j in grdi_models:
                    tests.append([j,'salm_amr','grdi','none','1000','11mer','False',5,i,l])
                    tests.append([j,'grdi','none','none','1000','11mer','False',5,i,l])
                    tests.append([j,'grdi','salm_amr','none','1000','11mer','False',5,i,l])
            # for figure, run all in 100-3k feats
            """
            for mdl in ['SVM','XGB','ANN']:
                ft = [str(i) for i in range (200,3000,100)]
                ft.remove('1000')
                for f in ft:
                    tests.append([mdl,'salm_amr','none','none',f,'11mer','False',5,i,l])

            #for drug in FIS AMP AMC AZM CHL CIP CRO FOX GEN NAL SXT TET TIO STR KAN; do for model in XGB SVM ANN; do for feat in {200..2900..100}; do acheron build model -x salm_amr -y none -l AMR_MIC -f $feat  -m $model -c 8 -a $drug -t 11mer --cluster slurm --trial 1; done: done; done
            """

    return tests

def find_missing_results(tests):
    """
    Checks which files are missing (and valid) and which need to be rerun.
    """
    #print("results/model={}_train={}_test={}_validate={}_feats={}_hyp={}_cvfolds={}_attribute={}_trial={}/summary.df".format(*tests[0]))
    missing = []
    missing_dirs = []
    for test in tests:
        # using try to make sure the file can be read and isnt corrupt
        try:
            file_path ="results/model={}_train={}_test={}_validate={}_feats={}_type={}_hyp={}_cvfolds={}_attribute={}_trial={}/summary.df".format(*test)
            results = pd.read_pickle(file_path)
        except:
            missing.append(test)
            missing_dirs.append(file_path)

    return missing, missing_dirs

def rerun_missing(missing, missing_dirs):
    """
    Executes acheron snakemake cluster jobs to build the missing files
    """
    # remove tests that are the same, just different trials
    # as acheron will run all missing trials from that specific test
    missing = [i[:-1] for i in missing]
    missing.sort()
    missing = list(i for i,_ in itertools.groupby(missing))

    missing = [i+['10'] for i in missing]

    print("running {} jobs".format(len(missing)))

    for miss, dir in zip(missing, missing_dirs):
        model,train,test,validate,feats,type,hyp,cvfolds,attribute,trial = miss
        if hyp == "True":
            command = "acheron build model -x {} -y {} -l AMR_MIC -f {} -m {} -c 8 -a {} -p -t 11mer --cluster slurm --trial 10".format(train,test,feats,model,attribute)
        else:
            command = "acheron build model -x {} -y {} -l AMR_MIC -f {} -m {} -c 8 -a {} -t 11mer --cluster slurm --trial 10".format(train,test,feats,model,attribute)
        print(command)
        # remove which jobs have already been completed, as to not corrupt the time and mem tests
        # actually, on 2nd thought, we can just check how many results came from each job and math accordingly
        # os.system("rm -rf {}".format(dir[:-11]))
        os.system(command)

if __name__ == "__main__":
    tests = steinkey2022()
    print("{} tests found".format(len(tests)))
    missing, missing_dirs = find_missing_results(tests)
    print("{} tests results were missing".format(len(missing)))
    rerun_missing(missing, missing_dirs)
