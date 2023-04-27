
""" Runs antigenicity prediction with euANTIGEN """

import os
import logging
import os
import numpy as np
import tensorflow
from tensorflow import keras
from keras.models import load_model
from tensorflow.keras.models import load_model
from tensorflow.python.ops.numpy_ops import np_config
from utils import *
from Protein import *


def euntigen(list_of_proteins, working_dir, NERVE_dir) -> list:

    logging.basicConfig(filename = os.path.join(working_dir, 'logfile.log'),
                        filemode = 'a',
                        level = logging.DEBUG,
                        force = True)

    model_dir = os.path.join(NERVE_dir, 'models/')
    std_devs = np.load(os.path.join(model_dir, 'std_devs_antigen.npy'))
    means = np.load(os.path.join(model_dir, 'means_antigen.npy'))
    projection_matrix = np.load(os.path.join(model_dir, 'projection_matrix_antigen.npy'))

    np_config.enable_numpy_behavior()
    model = load_model(os.path.join(model_dir, 'antigen.h5'), compile=False)


    config_antigen = model.get_config()
    model = tensorflow.keras.Model.from_config(config_antigen)


    for i, p in enumerate(list_of_proteins):
        data = p.standardize(means, std_devs, projection_matrix)
        data = np.array(data)[None, ...]
        p.p_antigen = float(model.predict(data, verbose=0))

    return list_of_proteins



