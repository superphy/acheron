import numpy as np
import pandas as pd

from hyperopt import hp, fmin, tpe, STATUS_OK, STATUS_FAIL, Trials, space_eval

from acheron.workflows import supervised_model

def PHACs_next_top_model(trials):
    """
    input: hyperopt trials object
    returns the trained model with the lowest loss
    """
    trials = [i for i in trials if i['result']['status'] == STATUS_OK]
    best_index = np.argmin([i['result']['loss'] for i in trials])
    best_model = trials[best_index]['result']['model']

    return best_model


def load_hyp_data(test_id):
    """
    Loads saved data split for hyperparameter optimizations
    used only in train_hyper_model()
    if test_id isnt in scope, will have to move to new script. use sys call
    """
    save_path = "data/hyp_data/{}/".format(test_id)

    x_train = pd.read_pickle(save_path+'0'+".pkl")
    y_train = list(np.load(save_path+'1'+".npy"))
    x_val = pd.read_pickle(save_path+'2'+".pkl")
    y_val = list(np.load(save_path+'3'+".npy"))

    return x_train, y_train, x_val, y_val

def create_model(x_train, y_train, x_test, y_test):
    """
    This is for hyperas, which isnt used anymore, dont try
    to use this anywhere cause it isnt python
    """
    patience = {{choice([4,8,12,16])}}
    early_stop = EarlyStopping(monitor='loss', patience=patience, verbose=0, min_delta=0.005, mode='auto')
    reduce_LR = ReduceLROnPlateau(monitor='loss', factor= 0.1, patience=(patience/2), verbose = 0, min_delta=0.005,mode = 'auto', cooldown=0, min_lr=0)

    model = Sequential()

    # how many hidden layers are in our model
    num_layers = {{choice(['zero', 'one', 'two', 'three', 'four', 'five'])}}

    if(num_layers == 'zero'):
        model.add(Dense(num_classes,activation='softmax',input_dim=(x_train.shape[1])))
    else:
        # this isnt a for loop because each variable needs its own name to be independently trained
        if (num_layers in ['one','two','three','four','five']):
            model.add(Dense(int({{uniform(num_classes,x_train.shape[1])}}),activation='relu',input_dim=(x_train.shape[1])))
            model.add(Dropout({{uniform(0,1)}}))
        if (num_layers in ['two','three','four','five']):
            model.add(Dense(int({{uniform(num_classes,x_train.shape[1])}})))
            model.add(Dropout({{uniform(0,1)}}))
        if (num_layers in ['three','four','five']):
            model.add(Dense(int({{uniform(num_classes,x_train.shape[1])}})))
            model.add(Dropout({{uniform(0,1)}}))
        if (num_layers in ['four','five']):
            model.add(Dense(int({{uniform(num_classes,x_train.shape[1])}})))
            model.add(Dropout({{uniform(0,1)}}))
        if (num_layers == 'five'):
            model.add(Dense(int({{uniform(num_classes,x_train.shape[1])}})))
            model.add(Dropout({{uniform(0,1)}}))

        model.add(Dense(num_classes, kernel_initializer='uniform', activation='softmax'))

    model.compile(loss='poisson', metrics=['accuracy'], optimizer='adam')
    model.fit(x_train, y_train, epochs=100, verbose=0, batch_size=6000, callbacks=[early_stop, reduce_LR])

    score, acc = model.evaluate(x_test, y_test, verbose=0)
    return {'loss': -acc, 'status': STATUS_OK, 'model': model}

def xgboost_objective(args):
    import xgboost as xgb

    x_train = args['x_train']
    y_train = args['y_train']
    x_val = args['x_val']
    y_val = args['y_val']

    for i in ['x_train', 'y_train', 'x_val', 'y_val']:
        args.pop(i)

    xgb_matrix = xgb.DMatrix(x_train.values, y_train, feature_names=x_train.columns)
    booster = xgb.train(args, xgb_matrix)

    prediction = supervised_model.predict(booster, x_val, 'XGB')

    acc = np.sum([i==j for i,j in zip(prediction, y_val)])/len(prediction)

    return {'loss': -acc, 'status': STATUS_OK, 'model': booster}

def neural_engine(args):
    from keras.utils import np_utils, to_categorical
    from keras.layers.core import Dense, Dropout, Activation
    from keras.models import Sequential
    from keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau

    x_train = args['x_train'].values
    y_train = args['y_train']
    x_val = args['x_val'].values
    y_val = args['y_val']

    train_cat_labels = to_categorical(y_train, args['num_class'])
    val_cat_labels = to_categorical(y_val, args['num_class'])

    for i in ['x_train', 'y_train', 'x_val', 'y_val']:
        args.pop(i)

    patience = args['patience']

    early_stop = EarlyStopping(monitor='loss', patience=patience, verbose=0, min_delta=0.005, mode='auto')
    reduce_LR = ReduceLROnPlateau(monitor='loss', factor= 0.1, patience=(patience/2), verbose = 0, min_delta=0.005,mode = 'auto', cooldown=0, min_lr=0)

    model = Sequential()

    # This is if there are no hidden layers
    # declare dense connection from input layer of size features, to output layer of size num classes
    if args['num_layers'] == 0:
        model.add(Dense(args['num_class'],activation='relu',input_dim=(x_train.shape[1])))
        model.add(Dropout(args['dropout0']))

    # This is if there is at least one hidden layer
    else:
        for hidden_layer in range(args['num_layers']):
            num_neurons = int(args['layer{}'.format(hidden_layer+1)])
            dropout_rate = args['dropout{}'.format(hidden_layer+1)]

            # Declare dense connections from input of size features, to a hidden layer of size num_neurons
            if hidden_layer == 0:
                model.add(Dense(num_neurons,activation='relu',input_dim=(x_train.shape[1])))
                model.add(Dropout(dropout_rate))

            # For every additional hidden layer, connect a new hidden layer of size
            # num_neurons to the output of the previous hidden layer
            else:
                model.add(Dense(num_neurons))
                model.add(Dropout(dropout_rate))

        # Now connect the end of the last hidden layer to an output layer of
        # size num_classes. We use dropout0 as the dropout prior to output
        model.add(Dense(args['num_class'], kernel_initializer='uniform', activation='softmax'))
        model.add(Dropout(args['dropout0']))

    model.compile(loss='poisson', metrics=['accuracy'], optimizer='adam')
    model.fit(x_train, train_cat_labels, epochs=100, verbose=0, batch_size=6000, callbacks=[early_stop, reduce_LR])

    score, acc = model.evaluate(x_val, val_cat_labels, verbose=0)
    return {'loss': -acc, 'status': STATUS_OK, 'model': model}


def get_best(x_train, y_train, x_val, y_val, model_type, num_classes):
    """
    Trains model using x/y_train, validating hyperparams using x/y_val

    Returns:
    best_model - the model with the best performance, according to x/y_val
    best_params - the parameters used to train this model
    """
    trials = Trials()

    """
    XGBoost
    """
    if model_type.upper() in ['XGB','XGBOOST']:
        if num_classes == 2:
            objective = 'binary:logistic'
        else:
            objective = 'multi:softmax'

        search_params = {
            'learning_rate':    hp.choice('learning_rate',    np.arange(0.05, 0.31, 0.05)),
            'max_depth':        hp.choice('max_depth',        np.arange(1, 8, 1, dtype=int)),
            'min_child_weight': hp.choice('min_child_weight', np.arange(1, 8, 1, dtype=int)),
            'colsample_bytree': hp.choice('colsample_bytree', np.arange(0.3, 0.8, 0.1)),
            'subsample':        hp.uniform('subsample', 0.8, 1),
            'objective':        objective,
            'num_class':        num_classes
            }

        objective_function = xgboost_objective

    elif model_type.upper() in ['ANN','KERAS','TF','TENSORFLOW']:
        num_feats = x_train.shape[1]
        search_params = {
            'patience':         hp.choice('patience',   [4,8,12,16]),
            'num_layers':       hp.choice('num_layers', np.arange(6)),
            'num_class':        num_classes
        }
        # up to 5 hidden layers, each with a dropout
        for layer in range(6):
            if layer == 0:
                search_params['dropout0'] = hp.uniform('dropout0', 0, 1)
            else:
                search_params['layer'+str(layer)] = hp.uniform('layer'+str(layer),num_classes, num_feats)
                search_params['dropout'+str(layer)] = hp.uniform('dropout'+str(layer), 0, 1)

        objective_function = neural_engine

    else:
        raise Exception("model type {} not defined".format(model_type))

    search_params['x_train'] = x_train
    search_params['y_train'] = y_train
    search_params['x_val'] = x_val
    search_params['y_val'] = y_val

    best_index = fmin(
        fn=objective_function, space=search_params,
        algo=tpe.suggest, max_evals=10, trials=trials)

    best_params = space_eval(search_params, best_index)

    best_model = PHACs_next_top_model(trials)

    return best_model, best_params
