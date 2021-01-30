#!/usr/bin/env python

import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import precision_recall_fscore_support
from collections import Counter
import math
import xgboost as xgb
import pickle
import sys, os

from acheron.helpers import model_evaluators

def is_valid_type(val):
    # covers core floats, numpy floating and complex floating
    if isinstance(val, (float, np.inexact)):
        if math.isnan(val) or np.isnan(val):
            return False
        else:
            return True
    # covers core ints, numpy ints (signed unsigned), longs, shorts, bytes
    # TODO check for bytes later, as they do NOT like math eqns when exponential
    elif isinstance(val, (int, np.integer)):
        return True
    # covers core and np strings, excludes unicode
    elif isinstance(val, (str, np.str_)):
        if val.upper() in ['INVALID', 'NAN']:
            return False
        else:
            return True
    else:
        return False

def make_mask(label_matrix, k):
    """
    Takes in a label matrix and saves a mask of which labels are valid
    for each index, this allows the model to keep a sample when only a single
    columns is invalid, instead of discarding the sample.
    If doing kfold cross validation (k != 1), there must be at least k samples
    in each class for that class to be considered valid
    """
    mask = np.zeros((label_matrix.values.shape), dtype='bool')

    invalid_classes = {}

    # if any class has fewer than k samples, it is marked as invalid
    for attribute in label_matrix.columns:
        class_counter = Counter(label_matrix[attribute])
        invalid_classes[attribute] = []
        for key in class_counter:
            if class_counter[key] < k:
                invalid_classes[attribute].append(key)

    for i, row in enumerate(label_matrix.values):
        for j, col in enumerate(row):
            if not is_valid_type(col):
                continue
            if col in invalid_classes[label_matrix.columns[j]]:
                continue
            mask[i,j] = True

    return pd.DataFrame(data=mask, columns = label_matrix.columns,
        index = label_matrix.index)

def make_split(label_matrix, mask, k, samples):
    """
    Takes a label matrix, splits it according to k fold cross validation
    but only for valid samples. Produces a random split each time (stratified)
    """
    assert(2<=k<=255) # if exceeding 255, change dtype below to uint16
    split_df =  pd.DataFrame(data=np.zeros(label_matrix.shape),
        columns=label_matrix.columns, index = label_matrix.index, dtype='uint8')

    for col in label_matrix.columns:
        # which labels are valid in this specific column
        valid_labels = label_matrix[col].values[mask[col].values]
        # matching sample name for each i in valid_labels
        valid_samples = label_matrix.index.values[mask[col].values]

        if len(valid_samples) == 0:
            print("All samples in column "+col+" are invalid, skipping split")
            continue

        # we also need to factor in that we only have the samples in /samples,
        # where a datasheet might have thousands of valid, but extra datapoints
        seen_bool_mask = np.array([i in samples for i in valid_samples])
        final_labels = valid_labels[seen_bool_mask]
        final_samples = valid_samples[seen_bool_mask]

        # at this point, we only have labels and samples that are eligible
        # for machine learning
        skf = StratifiedKFold(n_splits=k, shuffle=True)

        num_samples = len(final_samples)
        splits = enumerate(skf.split(np.zeros((num_samples,k)),final_labels))
        for i, split in splits:
            # pull the array of values assigned to the testing set,
            # label these genomes as per the fold they belong to
            for sample in final_samples[split[1]]:
                split_df.at[sample,col] = i

    return split_df


def load_data(dataset_name, label_name, trial, type, num_feats, k, attribute):
    """
    load requested dataset, mask, and split
    """
    features = pd.read_pickle("data/{}/features/{}_matrix.df".format(dataset_name, type))
    labels = pd.read_pickle("data/{}/labels/{}.df".format(dataset_name,label_name))[attribute]
    mask = pd.read_pickle("data/{}/features/masks/{}_{}.df".format(dataset_name,type,label_name))[attribute]

    if k!=1:
        split = pd.read_pickle("data/{}/masks/split{}_{}_{}_{}xCV.df".format(dataset_name,trial,type,label_name,k))
    else:
        split = []

    features, labels = apply_mask(features, labels, mask)

    return features, labels, mask, split

def train_model(features, label, model_type):
    """
    Converts feature and label matrices in a trained model
    """
    num_classes = len(Counter(label).keys())

    if model_type.upper() in ['XGB', 'XGBOOST']:
        if num_classes == 2:
            objective = 'binary:logistic'
        else:
            objective = 'multi:softmax'

        # this is is probably going to suck for memory, so lets revist XGBClassifier
        # if we explode past our ram usage on this step

        xgb_matrix = xgb.DMatrix(features.values, label, feature_names=features.columns)
        params = {'objective':objective, 'num_class': num_classes}
        #params = {'objective':objective}
        booster = xgb.train(params, xgb_matrix)
        return booster

    else:
        raise Exception("model type {} not defined".format(model_type))

def predict(model, features, model_type):
    """
    Takes a model and a feature set, returns an label like array of predictions
    """
    if model_type.upper() in ['XGB', 'XGBOOST']:
        xgb_matrix = xgb.DMatrix(features.values, feature_names = features.columns)
        return [round(i) for i in model.predict(xgb_matrix, validate_features=True)]
    else:
        raise Exception("model type {} not defined".format(model_type))

def evaluate_model(predicted, actual, model_type, dilutions, attribute, encoder):
    """
    Evaluates how well a model did (accuracy)
    For mic modules, also reports off-by-one accuracy and error rates
    Takes encoded class labels (0,1,2) not decoded values (2,4,8,16)
    """
    # this df will eventually get all info about the test
    direct_accuracy  = np.sum([predicted[i]==actual[i] for i in range(len(predicted))])/len(predicted)

    dilutional_accuracies = {}
    find_errors = False

    if len(dilutions) > 0:
        find_errors = True
        for dilution in dilutions:
            total = 0
            correct = 0
            for i in range(len(predicted)):
                total +=1
                if abs(predicted[i]-actual[i]) <= dilution:
                    correct +=1
            dilutional_accuracies[dilution] = correct/total

    data = [direct_accuracy]
    columns = ['Accuracy']

    for dilution in dilutions:
        if str(dilution) == '0':
            continue
        else:
            data.append(dilutional_accuracies[dilution])
            columns.append("Within {} Dilution".format(dilution))

    if find_errors:
        decoder = {v:k for k,v in encoder.items()}
        pred_decoded = [decoder[i] for i in predicted]
        act_decoded = [decoder[i] for i in actual]
        errors = [model_evaluators.find_error_type(i[0],i[1], attribute) for i in zip(pred_decoded, act_decoded)]

        error_counts = Counter(errors)

        error_types = ["Very Major Error", "Major Error", "Non Major Error", "Correct"]
        total_errors = 0

        for error_type in error_types:
            total_errors += error_counts[error_type]

            percent = error_counts[error_type]/len(predicted)

            data.append(percent)
            columns.append(error_type)

        try:
            assert len(predicted) == total_errors
        except:
            print('Number of Errors+Correct does not equal number of predictions')
            raise

    results_df = pd.DataFrame(data=[data], columns=columns)

    return results_df

def apply_mask(features, labels, mask):
    """
    Takes in a pandas dataframe or series with a mask
    and returns only valid samples
    Mask series looks like:
                    AMC
    BioSample
    SAMN00000001    False
    SAMN00000002    True
    """
    # its important to note that the mask is for everything in the label df, but
    # when we mask the features, that df might not have it
    # therefore we reduce the mask to samples that are seen
    seen = list(features.index)

    mask = mask[[i in seen for i in mask.index]]
    labels = labels[[i in seen for i in labels.index]]

    # prior to reindexing, we need to make sure there are no duplicates
    mask = mask[~mask.index.duplicated(keep='first')]
    labels = labels[~labels.index.duplicated(keep='first')]

    # reorder dataframes to make sure they are in the same order as the features
    mask = mask.reindex(seen)
    labels = labels.reindex(seen)

    # remove samples that we dont have data for
    # (these appear as nans in the mask)
    mask = pd.Series(index = mask.index,
        data = [False if np.isnan(i) else i for i in mask])

    # double check that the mask biosample order matches the feature biosample order
    for i, biosample in enumerate(mask.index):
        assert biosample == seen[i]

    if isinstance(features, pd.Series) or isinstance(features, pd.DataFrame):
        labels = labels[list(mask)]
        features = features[list(mask)]

    else:
        raise Exception("Masking of type {} is not defined".format(type(df)))

    return features, labels

def select_features(features, labels, num_feats, pre_selection):
    """
    Reduces feature matrix down to the top num_feats features,
    returns filtered feature matrix

    If an array of features if provided in pre_selection,
    returns matrix with only columns that match pre_selection, and in that order
    """
    # we are choosing the top features
    if isinstance(pre_selection, bool) and pre_selection == False:
        # save the non values, which are lost in SelectKBest
        indeces = features.index
        cols = features.columns

        sk_obj = SelectKBest(f_classif, k=num_feats)

        # reduce features to top x features
        features = sk_obj.fit_transform(features,labels)

        # remove the columns that were removed from features,
        # this needs to be saved for when we filter the testing set
        cols = sk_obj.transform([cols])[0]

        features = pd.DataFrame(data=features,columns=cols,index=indeces)

        return features

    else:
        # we are using pre-determined features
        # needs to be list or pandas.columns
        assert isinstance(pre_selection, list) or isinstance(pre_selection, pd.core.indexes.base.Index)

        # Only keep important data (this also reorders)
        features = features[pre_selection]

        return features



def make_model(model_type,train,test,validation,label_name,type,attribute,num_feats,do_hyp,cv_folds,trial):
    """
    Loads data, removes invalid data, splits data, feeds into ML model.

    NOTE: for hyperparameter optimization, the model is trained using the TRAINING set,
    the hyperparameters are optimized and validated using the VALIDATION set, and the final
    trained model is tested using the TESTING set. You CANNOT use the testing set anywhere else.
    You may split the training set for validation by not passing validation.

    Possible configs for NON hyperparameter optimization:
    train is passed: will cross-validate by splitting into 80% train, 20% test
    train/test is passed: trains using 100% of train, tests using 100% of test.
    validation passed: Will throw error, as validation set cant be used.


    Possible configs for hyperparameter optimization: (nested cross-validation)
    train is passed: Splits in 60%-20%-20% for train/validation/testing respectively
    train/test passed: splits train into 80% train, 20% validation, uses testing set for testing
    train/validation: splits off training data for testtraining
    train/test/validation passed: trains on train, validates hyperparams using validation set,
    testing of the final model is done using the testing set.
    error if test set is the same as train or validation
    """
    with open("data/{}/labels/{}_encoder.pkl".format(train, label_name), 'rb') as unpickler:
        train_encoder = pickle.load(unpickler)

    if 'MIC' in label_name:
        dilutions = [1]
    else:
        dilutions = []

    if not do_hyp:
        # for non-nested cross validation[attribute]
        if validation not in ['none','None']:
            raise Exception("Validation set can only be used when doing hyperparameter optimizations")
        if test == '':
            # k-fold cv split training set into train and test
            # x is data, y is labels, splits are which samples belong in each split for each cv
            # dataset_name, label_name, trial, type, num_feats, k
            features, labels, mask, split = load_data(train,label_name,trial,type,num_feats,k)
            for fold in cv_folds:

                pass
                # TODO do split stuff here, apply mask, etc


        else:
            # call with supplied train and test
            # splits in this case will be a matrix of 1's, as there is only 1 fold[attribute]
            train_features, train_labels, train_mask, train_split = load_data(train,label_name,trial,type,num_feats,1,attribute)
            test_features, test_labels, test_mask, test_split = load_data(test,label_name,trial,type,num_feats,1,attribute)

            with open("data/{}/labels/{}_encoder.pkl".format(test, label_name), 'rb') as test_unpickler:
                test_encoder = pickle.load(test_unpickler)

            encoded_train_labels = [train_encoder[attribute][i] for i in train_labels]

            # reduce to top features
            train_features = select_features(train_features, train_labels, num_feats, False)
            test_features = select_features(test_features,'dont need',num_feats, train_features.columns)

            model = train_model(train_features, encoded_train_labels, model_type)
            predicted = predict(model, test_features, model_type)
            encoded_test_labels = [test_encoder[attribute][i] for i in test_labels]
            summary = evaluate_model(predicted, encoded_test_labels, model_type, dilutions, attribute, train_encoder[attribute])
            prec_recall = precision_recall_fscore_support(encoded_test_labels, predicted,average=None, labels=list(test_encoder[attribute].values()))
            prec_recall = np.transpose(prec_recall)
            prec_recall = pd.DataFrame(data=prec_recall, index = list(train_encoder[attribute].keys()),columns = ['Precision','Recall', 'F-Score','Supports'])

            decoder = {v:k for k,v in test_encoder[attribute].items()}
            decoded_predictions = [decoder[i] for i in predicted]
            predicted_df = pd.DataFrame(data=decoded_predictions, columns=[attribute], index=test_features.index)

    else:
        # for nested cross-validation (hyperparameter optimizations)
        if test in [train,validation]:
            raise Exception("You cannot use the testing set for training or validation")
        if validation == 'none' and test == 'none':
            #TODO call 3way split, nested cross validation on train
            pass
        elif validation == 'none':
            #TODO split train 80% train, 20% validation, used passed train
            pass
        elif test == 'none':
            #TODO split train 80% train, 20% testing, use passed validation
            pass
        else:
            # All datasets passed
            if train == validation:
                raise Exception("Same dataset used for training and validation, to split training set for validation, do not pass --validation")
            #TODO do no splitting, call for model execution
            pass

    # if the features are not included in a model, they need to be included somewhere
    return model, predicted_df, summary, prec_recall

def manual_call(arguments):
    # attr = getattr(arguments,arg)

    # which label to use
    label = getattr(arguments, 'label')
    # how many splits/folds in cross validation, e.g. 5-fold CV
    num_cv_splits = getattr(arguments, 'cv')

    trials = [i for i in range(getattr(arguments, 'trial'))]
    attribute_type = getattr(arguments, 'type')

    # First select what data needs to be setup
    datasets = []
    for dataset in ['train','test','validation']:
        if getattr(arguments, dataset) not in ['None','none']:
            datasets.append(dataset)

    # Make masks (selects which data is valid and which samples have attributes to be ignored)
    for dataset in datasets:
        label_df = pd.read_pickle("data/{}/labels/{}.df".format(dataset,label))
        if len(datasets) > 1:
            mask = make_mask(label_df, num_cv_splits)
        else:
            mask = make_mask(label_df, 1)
        mask.to_pickle("data/{}/features/masks/{}_{}.df".format(
            dataset,attribute_type,label))

    # Split data (only required when using one dataset (cross validating))
    if len(datasets) == 1:
        import glob
        samples = glob.glob("data/{}/wgs/raw/*.fasta".format(datasets[0]))
        for trial in trials:
            split = make_split(label_df, mask, num_cv_splits, samples)
            split.to_pickle("data/{}/masks/split{}_{}_{}_{}xCV.df".format(
                datasets[0],trial,attribute_type,label,num_cv_splits))

    # Build Model
    for trial in trials:
        model, predicted_df, summary, prec_recall = supervised_model.make_model(
            getattr(arguments,'model'),getattr(arguments,'train'),
            getattr(arguments,'test'),getattr(arguments,'validation'),
            getattr(arguments,'label'),getattr(arguments,'type'),
            getattr(arguments,'attribute'),getattr(arguments,'num_features'),
            getattr(arguments,'hyperparam'),getattr(arguments,'cv'),trial)

        out = "results/model={}_train={}_test={}_validate={}_feats={}_hyp={}_cvfolds={}_attribute={}_trial={}/".format(
        getattr(arguments,'model'), getattr(arguments,'train'), getattr(arguments,'test'),
        getattr(arguments,'validation'), getattr(arguments,'num_features'),
        getattr(arguments,'hyperparam'), getattr(arguments,'cv'),
        getattr(arguments,'attribute'), trial)

        joblib.dump(model, out+'/model.joblib')
        predicted_df.to_pickle(out+"/predictions.df")
        prec_recall.to_pickle(out+"/precision_recall_fscore_support.df")
        summary.to_pickle(out+"/summary.df")

        print("Model creation completed and saved. Results saved in {}".format(out))

if __name__ == "__main__":
    """
    If we know that the inputs are satisfied we can manually call this script.
    This is done using acheron, not by calling python acheron/workflows/supervised_model.py
    WARNING: This directly calls model creation, things like cluster support will be skipped!
    If using slurm, you need to wrap the call yourself.
    """
    raise Exception("To manually call this script, pass --manual to acheron build model. This allows you to still use the command line tool parser")
