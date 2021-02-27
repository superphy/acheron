import numpy as np
from numpy.random import seed
import pandas as pd
from pandas import DataFrame
import sys
import pickle
from decimal import Decimal
import os

import tensorflow
from tensorflow import set_random_seed

from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

from hyperopt import Trials, STATUS_OK, tpe
from keras.layers.convolutional import Conv1D
from keras.layers.core import Dense, Dropout, Activation
from keras.layers import Flatten, BatchNormalization
from keras.models import Sequential, load_model
from keras.utils import np_utils, to_categorical
from keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau

from keras import backend as K
K.set_session(K.tf.Session(config=K.tf.ConfigProto(intra_op_parallelism_threads=16, inter_op_parallelism_threads=16)))

from sklearn import metrics
from sklearn.externals import joblib
from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit
from sklearn.feature_selection import SelectKBest, f_classif

def create_model(x_train, y_train, x_test, y_test):
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
