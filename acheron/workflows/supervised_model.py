#!/usr/bin/env python

import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import precision_recall_fscore_support, mean_squared_error
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
    split_df = pd.DataFrame(data=np.zeros(label_matrix.shape),
        columns=label_matrix.columns, index = label_matrix.index, dtype='uint8')

    split_df = split_df[~split_df.index.duplicated(keep='first')]

    for col in label_matrix.columns:
        # which labels are valid in this specific column
        valid_labels = label_matrix[col].values[mask[col].values]
        # matching sample name for each i in valid_labels
        valid_samples = label_matrix.index.values[mask[col].values]

        if len(valid_samples) == 0:
            print("All samples in column "+col+" are invalid, skipping split")
            continue

        # in the event of duplicates, keep only the first seen instance
        processed = []

        # we also need to factor in that we only have the samples in /samples,
        # where a datasheet might have thousands of valid, but extra datapoints
        #seen_bool_mask = np.array([i in samples and i not in duplicates for i in valid_samples])
        seen_bool_mask = []

        for i in valid_samples:
            if i in processed:
                seen_bool_mask.append(False)
            else:
                processed.append(i)
                if i in samples:
                    seen_bool_mask.append(True)
                else:
                    seen_bool_mask.append(False)

        seen_bool_mask = np.array(seen_bool_mask)

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
        split = pd.read_pickle("data/{}/splits/split{}_{}_{}_{}xCV.df".format(dataset_name,trial,type,label_name,k))
    else:
        split = []

    features, labels = apply_mask(features, labels, mask)

    return features, labels, mask, split

def train_model(features, label, model_type, num_classes):
    """
    Converts feature and label matrices in a trained model
    Sci-kit models are at the end, as they share a fit method
    """

    # XGBoost
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

    # Artificial Neural Network
    elif model_type.upper() in ['ANN']:
        from keras.utils import np_utils, to_categorical
        from keras.layers.core import Dense, Dropout, Activation
        from keras.models import Sequential
        from keras.utils import np_utils, to_categorical
        from keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau

        cat_labels = to_categorical(label, num_classes)

        patience = 16
        early_stop = EarlyStopping(monitor='loss', patience=patience, verbose=1, min_delta=0.005, mode='auto')
        reduce_LR = ReduceLROnPlateau(monitor='loss', factor= 0.1, patience=(patience/2), verbose = 1, min_delta=0.005,mode = 'auto', cooldown=0, min_lr=0)

        num_feats = len(features.columns)
        model = Sequential()
        model.add(Dense(int(((num_feats+num_classes)/2)),activation='relu',input_dim=(num_feats)))
        model.add(Dropout(0.5))
        model.add(Dense(num_classes, kernel_initializer='uniform', activation='softmax'))

        if num_classes == 2:
            loss = "binary_crossentropy"
        else:
            loss = "poisson"

        model.compile(loss=loss, metrics=['accuracy'], optimizer='adam')
        model.fit(features.values, cat_labels, epochs=100, verbose=1, callbacks=[early_stop, reduce_LR])
        return model

    # Support Vector Classifier
    # https://scikit-learn.org/stable/modules/svm.html#classification
    if model_type.upper() in ['SVC']:
        from sklearn import svm
        model = svm.SVC()

    # Support Vector Regressor
    # https://scikit-learn.org/stable/modules/svm.html#regression
    elif model_type.upper() in ['SVR']:
        from sklearn import svm
        model = svm.SVR()

    # Stochastic Gradient Descent Classifier
    # https://scikit-learn.org/stable/modules/sgd.html#classification
    elif model_type.upper() in ['SGDC']:
        from sklearn.linear_model import SGDClassifier
        model = SGDClassifier(loss="hinge", penalty="l2", max_iter=25)

    # Perceptron
    # https://scikit-learn.org/stable/modules/linear_model.html#perceptron
    elif model_type.upper() in ['PERC']:
        from sklearn.linear_model import SGDClassifier
        model = SGDClassifier(loss="perceptron", eta0=1, learning_rate="constant", penalty=None)

    # Passive Aggressive Algorithms
    # https://scikit-learn.org/stable/modules/linear_model.html#passive-aggressive-algorithms
    # https://www.geeksforgeeks.org/passive-aggressive-classifiers
    elif model_type.upper() in ['PAC']:
        from sklearn.linear_model import PassiveAggressiveClassifier
        model = PassiveAggressiveClassifier(max_iter=100)

    # Nearest Neighbours Classifier
    # https://scikit-learn.org/stable/modules/neighbors.html#nearest-neighbors-classification
    elif model_type.upper() in ['NNC']:
        from sklearn.neighbors import KNeighborsClassifier
        model = KNeighborsClassifier(n_neighbors=3)

    # Nearest Neighbours Regressor
    # https://scikit-learn.org/stable/modules/neighbors.html#nearest-neighbors-regression
    elif model_type.upper() in ['NNR']:
        from sklearn.neighbors import KNeighborsRegressor
        model = KNeighborsRegressor(n_neighbors=3)

    # Gaussian Naive Bayes
    # https://scikit-learn.org/stable/modules/naive_bayes.html#gaussian-naive-bayes
    elif model_type.upper() in ['GNB']:
        from sklearn.naive_bayes import GaussianNB
        model = GaussianNB()

    # Multinomial Naive Bayes
    # https://scikit-learn.org/stable/modules/naive_bayes.html#multinomial-naive-bayes
    elif model_type.upper() in ['MNB']:
        from sklearn.naive_bayes import MultinomialNB
        model = MultinomialNB()

    # Categorical Naive Bayes
    # https://scikit-learn.org/stable/modules/naive_bayes.html#categorical-naive-bayes
    elif model_type.upper() in ['CNB']:
        from sklearn.naive_bayes import CategoricalNB
        model = CategoricalNB()

    # Decision Tree Classifier
    # https://scikit-learn.org/stable/modules/tree.html#classification
    elif model_type.upper() in ['DTC']:
        from sklearn import tree
        model = tree.DecisionTreeClassifier()

    # Decision Tree Regressor
    # https://scikit-learn.org/stable/modules/tree.html#regression
    elif model_type.upper() in ['DTR']:
        from sklearn import tree
        model = tree.DecisionTreeRegressor()

    # AdaBoost
    # https://scikit-learn.org/stable/modules/ensemble.html#adaboost
    elif model_type.upper() in ['ADA']:
        from sklearn.ensemble import AdaBoostClassifier
        model = AdaBoostClassifier()

    # Gradient Boosted Decision Trees
    # https://scikit-learn.org/stable/modules/ensemble.html#gradient-tree-boosting
    elif model_type.upper() in ['GBDT']:
        from sklearn.ensemble import GradientBoostingClassifier
        model = GradientBoostingClassifier()

    # Multi-layer Perceptron Classifier
    # https://scikit-learn.org/stable/modules/neural_networks_supervised.html#multi-layer-perceptron
    elif model_type.upper() in ['MLPC']:
        from sklearn.neural_network import MLPClassifier
        model = MLPClassifier()

    else:
        raise Exception("model type {} not defined".format(model_type))

    model.fit(features.values, label)
    return model

def train_hyper_model(x_train, y_train, x_val, y_val, model_type, num_classes):
    """
    Trains a hyperparameter optimized model
    """
    from hyperopt import hp, fmin, tpe, STATUS_OK, STATUS_FAIL, Trials
    from acheron.workflows import hyp
    from hyperas import optim

    # https://towardsdatascience.com/an-example-of-hyperparameter-optimization-on-xgboost-lightgbm-and-catboost-using-hyperopt-12bc41a271e
    # Search Space Subject to Change!!
    if model_type.upper() in ['XGB','XGBOOST']:
        params = {
                'learning_rate':    hp.choice('learning_rate',    np.arange(0.05, 0.31, 0.05)),
                'max_depth':        hp.choice('max_depth',        np.arange(1, 8, 1, dtype=int)),
                'min_child_weight': hp.choice('min_child_weight', np.arange(1, 8, 1, dtype=int)),
                'colsample_bytree': hp.choice('colsample_bytree', np.arange(0.3, 0.8, 0.1)),
                'subsample':        hp.uniform('subsample', 0.8, 1),
                'n_estimators':     100,
                }
        pass

    elif model_type.upper() in ['ANN','KERAS','TF','TENSORFLOW']:
        best_run, best_model = optim.minimize(
            model=hyp.create_model(x_train, y_train, x_test, y_test),
            data=(x_train, y_train, x_test, y_test),
            algo=tpe.suggest,
            max_evals=10,
            trials=Trials(),
            keep_temp=True)

        return best_model, best_run

    else:
        raise Exception("model type {} not defined".format(model_type))

    """
    # Minimal Cost-Complexity Pruning
    # https://scikit-learn.org/stable/modules/tree.html#minimal-cost-complexity-pruning
    elif model_type.upper() in ['MCCP']:
        from sklearn import tree
        model = tree.DecisionTreeClassifier()
        path = model.cost_complexity_pruning_path(features.values, labels)
        ccp_alphas, impurities = path.ccp_alphas, path.impurities
        models = []
        for ccp_alpha in ccp_alphas:
            clf = DecisionTreeClassifier(ccp_alpha=ccp_alpha)
            clf.fit(features.values, labels)
            clfs.append(clf)
        # now use validation set to see which model did best, use that alpha to train final model
    """

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

    data = [len(predicted),direct_accuracy]
    columns = ["Supports", "Accuracy"]

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

def mean_summaries(summaries):
    """
    Takes a list of model summaries
    averages them appropriately, relevant to number of supports
    """
    try:
        indx = summaries[0].index[0]
    except:
        indx = 0

    mean_df = pd.DataFrame(columns=summaries[0].columns, index=[indx])

    total_supports = 0
    proportion = {}

    for summary in summaries:
        num_sups = summary['Supports'][0]
        total_supports += num_sups
        for col in summary.columns:
            if col != "Supports":
                if col in proportion.keys():
                    proportion[col] += num_sups*summary[col][0]
                else:
                    proportion[col] = num_sups*summary[col][0]

    mean_df.loc[indx,'Supports'] = total_supports
    for k,v in proportion.items():
        mean_df.loc[indx, k] = v/total_supports

    return mean_df

def mean_prec_recall(prec_recall_dfs):
    """
    Takes a list of precision_recall_fscore_support dataframes
    Returns the mean df based on proportion of supports
    """
    indeces = prec_recall_dfs[0].index
    done_rows = []

    for indx in indeces:
        rows = [i[i.index == indx] for i in prec_recall_dfs]
        done_rows.append(mean_summaries(rows))

    return pd.concat(done_rows)

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

def split_data(features, labels, split, attribute, hyp, fold, num_splits, BLOCK_TRAIN=False):
    """
    Splits into train, testing and validation (if asked for)
    based on predetermined spliting computed beforehand
    No determination of which samples belong in which split happens at this point
    """

    # if BLOCK_TRAIN = true, folds are batched into blocks, so 10-fold would be 2 test, 2 val, 6 train
    # if left as false, its 1 fold test, 1 fold val, and the rest would be training, regardless to number of folds

    split = split[attribute]

    # ensure splits are in the same order as the features
    split = split[split.index.isin(features.index)]
    split = split.reindex(features.index)

    split_priority = [(i+int(fold)-1)%num_splits for i in range(num_splits)]
    """
    split_priority above defines which folds become validation or testing
    if fold 3 for a 5-fold cross-validation is given, split_priority will look like:
        [2, 3, 4, 0, 1]
    where fold 2 becomes the testing set, and 3 becomes validation (if requested, else becomes training)
    if fold 7 for 10-fold cross-validation is given, split_priority will look like:
        [6, 7, 8, 9, 0, 1, 2, 3, 4, 5]
    """

    if BLOCK_TRAIN:
        # If using 5 fold cross val, testing will be 1 block or 1 of 5
        # if using 10 fold, 20% will therefore be 2 blocks
        # 7 fold with then be 1 block for test 1 block val, 5 blocks training
        if num_splits < 5:
            folds_per_block = 1
        else:
            folds_per_block = num_splits//5
    else:
        folds_per_block = 1

    # which block of fold is to be used for testing
    test_blocks = split_priority[:folds_per_block]

    # Split for cross-validation, e.g. 80%-20% for 5 folds
    if hyp == False:
        test_mask = [i in test_blocks for i in split]
        train_mask = [i not in test_blocks for i in split]

        x_train = features[train_mask]
        y_train = labels[train_mask]

        x_test = features[test_mask]
        y_test = labels[test_mask]

        return x_train, y_train, x_test, y_test

    # Split for nested cross-validation, e.g. 60%-20%-20%, for 5 folds
    else:
        # which blocks to use for validation
        val_blocks = split_priority[folds_per_block:folds_per_block*2]

        test_mask = [i in test_blocks for i in split]
        val_mask = [i in val_blocks for i in split]
        train_mask = [i not in test_blocks and i not in val_blocks for i in split]

        x_train = features[train_mask]
        y_train = labels[train_mask]

        x_test = features[test_mask]
        y_test = labels[test_mask]

        x_val = features[val_mask]
        y_val = labels[val_mask]

        return x_train, y_train, x_val, y_val, x_test, y_test

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

    num_classes = len(train_encoder[attribute].keys())

    if 'MIC' in label_name or 'SIR' in label_name:
        dilutions = [1]
    else:
        dilutions = []

    # As models are created, we can save them here, average after
    final_models = []
    final_features = []
    final_labels = []


    if not do_hyp:
        """
        This block for k-fold cross validation (NOT nested)
        And direct dataset-to-dataset predictions
        """
        if validation not in ['none','None']:
            raise Exception("Validation set can only be used when doing hyperparameter optimizations")
        if test in ['none','None']:
            # k-fold cv split training set into train and test
            # x is data, y is labels, splits are which samples belong in each split for each cv
            # dataset_name, label_name, trial, type, num_feats, k
            features, labels, mask, split = load_data(train,label_name,trial,type,num_feats,cv_folds, attribute)
            for fold in range(cv_folds):
                # reduce to top features,
                # also only return the features that have been selected to be part of each set
                x_train, y_train, x_test, y_test = split_data(features, labels, split, attribute, do_hyp, fold, cv_folds)

                # select features based purely on training data, testing data cant be seen by feature selection
                x_train = select_features(x_train, y_train, num_feats, False)
                x_test = select_features(x_test,'dont need',num_feats, x_train.columns)

                encoded_y_train = [train_encoder[attribute][i] for i in y_train]

                model = train_model(x_train, encoded_y_train, model_type, num_classes)
                final_models.append(model)
                final_features.append(x_test)
                final_labels.append(y_test)

        else:
            # call with supplied train and test, no cross validation
            # splits in this case will be a matrix of 1's, as there is only 1 fold[attribute]
            x_train, y_train, train_mask, train_split = load_data(train,label_name,trial,type,num_feats,1,attribute)
            x_test, y_test, test_mask, test_split = load_data(test,label_name,trial,type,num_feats,1,attribute)

            with open("data/{}/labels/{}_encoder.pkl".format(test, label_name), 'rb') as test_unpickler:
                test_encoder = pickle.load(test_unpickler)

            encoded_y_train = [train_encoder[attribute][i] for i in y_train]

            # reduce to top features
            x_train = select_features(x_train, y_train, num_feats, False)
            x_test = select_features(x_test,'dont need',num_feats, y_train.columns)

            model = train_model(x_train, encoded_y_train, model_type, num_classes)
            final_models.append(model)
            final_features.append(x_test)
            final_labels.append(y_test)


    else:
        """
        This block for nested cross-validation (hyperparameter optimizations)

        There are 4 different ways this test could be shaped, the last 9 lines of code are all the same.
        If you want to write a new test, make sure the parameters of the repeating 9 lines are satisfied.

        The reason the 9 lines are outside of the loop is because different tests have the data fractioned
        differently. To avoid flagging memory limits on cluster managers, we need internally split datasets
        to be loaded as efficiently as possible and therefore different scope as other dataset.
        """

        if test in [train,validation] and test != 'none':
            raise Exception("You cannot use the testing set for training or validation")
        if validation == 'none' and test == 'none':
            """
            In this block, the training set is split into 60% training data,
            20% hyperparam option validation, 20% testing (for 5-fold cv)
            """
            x, y, mask, split = load_data(train,label_name,trial,type,num_feats,cv_folds,attribute)
            for fold in range(cv_folds):

                x_train, y_train, x_val, y_val, x_test, y_test = split_data(
                    x, y, split, attribute, do_hyp, fold, cv_folds)

                x_train = select_features(x_train, y_train, num_feats, False)
                x_val = select_features(x_val,'dont need',num_feats, x_train.columns)
                x_test = select_features(x_test,'dont need',num_feats, x_train.columns)

                encoded_y_train = [train_encoder[attribute][i] for i in y_train]
                encoded_y_val = [train_encoder[attribute][i] for i in y_val]

                model = train_hyper_model(x_train, encoded_y_train, x_val, encoded_y_val, model_type, num_classes)
                final_models.append(model)
                final_features.append(x_test)
                final_labels.append(y_test)

        elif validation == 'none':
            """
            In this block, the training set is split into 80% training data,
            20% hyperparam option validation, and the entirety of the dataset passed in as test
            is used to test each fold. For 5-fold cross-validation, each 20% of the training set
            gets used as a validation set, each time tested using 100% of the test set.
            """
            # split train into 80% train, 20% validation, used passed test for testing
            x, y, mask, split = load_data(train,label_name,trial,type,num_feats,cv_folds,attribute)
            x_test, y_test, test_mask, test_split = load_data(test,label_name,trial,type,num_feats,cv_folds,attribute)

            for fold in range(cv_folds):
                # call with do_hype=False so that it only splits 80%-20%
                x_train, y_train, x_val, y_val = split_data(
                    x, y, split, attribute, False, fold, cv_folds)

                x_train = select_features(x_train, y_train, num_feats, False)
                x_val = select_features(x_val,'dont need',num_feats, x_train.columns)
                x_test = select_features(x_test,'dont need',num_feats, x_train.columns)

                encoded_y_train = [train_encoder[attribute][i] for i in y_train]
                encoded_y_val = [train_encoder[attribute][i] for i in y_val]

                model = train_hyper_model(x_train, encoded_y_train, x_val, encoded_y_val, model_type, num_classes)
                final_models.append(model)
                final_features.append(x_test)
                final_labels.append(y_test)

        elif test == 'none':
            """
            In this block, the training set is split into 80% training data,
            20% testing data (for 5-fold cv). This is repeated 5 times, each time validating
            hyperparams using the passed validation set.
            """
            x, y, mask, split = load_data(train,label_name,trial,type,num_feats,cv_folds,attribute)
            x_val, y_val, val_mask, val_split = load_data(val,label_name,trial,type,num_feats,cv_folds,attribute)

            for fold in range(cv_folds):
                # call with do_hype=False so that it only splits 80%-20%
                x_train, y_train, x_test, y_test = split_data(
                    x, y, split, attribute, False, fold, cv_folds)

                x_train = select_features(x_train, y_train, num_feats, False)
                x_val = select_features(x_val,'dont need',num_feats, x_train.columns)
                x_test = select_features(x_test,'dont need',num_feats, x_train.columns)

                encoded_y_train = [train_encoder[attribute][i] for i in y_train]
                encoded_y_val = [train_encoder[attribute][i] for i in y_val]

                model = train_hyper_model(x_train, encoded_y_train, x_val, encoded_y_val, model_type, num_classes)
                final_models.append(model)
                final_features.append(x_test)
                final_labels.append(y_test)
        else:
            # All datasets passed
            if train == validation:
                raise Exception("Same dataset used for training and validation, to split training set for validation, do not pass --validation")
            """
            In this block, no splitting occurs.
            3 passed datasets, each used for its labeled purpose
            """
            x_train, y_train, train_mask, train_split = load_data(train,label_name,trial,type,num_feats,cv_folds,attribute)
            x_val, y_val, val_mask, val_split = load_data(val,label_name,trial,type,num_feats,cv_folds,attribute)
            x_test, y_test, test_mask, test_split = load_data(test,label_name,trial,type,num_feats,cv_folds,attribute)

            x_train = select_features(x_train, y_train, num_feats, False)
            x_val = select_features(x_val,'dont need',num_feats, x_train.columns)
            x_test = select_features(x_test,'dont need',num_feats, x_train.columns)

            encoded_y_train = [train_encoder[attribute][i] for i in y_train]
            encoded_y_val = [train_encoder[attribute][i] for i in y_val]

            model = train_hyper_model(x_train, encoded_y_train, x_val, encoded_y_val, model_type, num_classes)
            final_models.append(model)
            final_features.append(x_test)
            final_labels.append(y_test)

    trained_models = []
    predicted_dfs = []
    summaries = []
    prec_recall_dfs = []

    for model, x_test, y_test in zip(final_models, final_features, final_labels):

        predicted = predict(model, x_test, model_type)
        encoded_y_test = [train_encoder[attribute][i] for i in y_test]
        summary = evaluate_model(predicted, encoded_y_test, model_type, dilutions, attribute, train_encoder[attribute])
        prec_recall = precision_recall_fscore_support(encoded_y_test, predicted,average=None, labels=list(train_encoder[attribute].values()))
        prec_recall = np.transpose(prec_recall)
        prec_recall = pd.DataFrame(data=prec_recall, index = list(train_encoder[attribute].keys()),columns = ['Precision','Recall', 'F-Score','Supports'])

        decoder = {v:k for k,v in train_encoder[attribute].items()}
        decoded_predictions = [decoder[i] for i in predicted]
        predicted_df = pd.DataFrame(data=decoded_predictions, columns=[attribute], index=x_test.index)

        trained_models.append(model)
        predicted_dfs.append(predicted_df)
        summaries.append(summary)
        prec_recall_dfs.append(prec_recall)

    # Models stay in a list

    # stack predictions
    predicted_dfs = pd.concat(predicted_dfs)

    # mean summaries
    summaries = mean_summaries(summaries)

    # mean prec_recall_dfs
    prec_recall_dfs = mean_prec_recall(prec_recall_dfs)

    return trained_models, predicted_dfs, summaries, prec_recall_dfs

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
        models, predicted_df, summary, prec_recall = supervised_model.make_model(
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

        for fold_num, model in enumerate(models):
            joblib.dump(model, "{}/model{}.joblib".format(out,fold))
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
