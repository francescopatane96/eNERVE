
""" Runs adhesin or adhesin-like probability prediction with euSPAAN (modified version) """

import tensorflow
from tensorflow.keras.models import load_model
from tensorflow.python.ops.numpy_ops import np_config
from Utils import *
from Protein import *

def euspaan(list_of_proteins, working_dir, NERVE_dir) -> list:
    """Run adhesin and adhesin-like predictor"""

    logging.basicConfig(filename = os.path.join(working_dir, 'logfile.log'),
                        filemode = 'a',
                        level = logging.DEBUG,
                        force = True)

    model_dir = os.path.join(NERVE_dir, 'models/')
    std_devs = np.load(os.path.join(model_dir, 'std_devs_adh.npy'))
    means = np.load(os.path.join(model_dir, 'means_adh.npy'))
    projection_matrix = np.load(os.path.join(model_dir, 'projection_matrix_adh.npy'))
    #np_config.enable_numpy_behavior()
    model = load_model(os.path.join(model_dir, 'adhesin.h5'), compile=False)

    #config_adhesin = model.get_config()
    #model = tensorflow.keras.Model.from_config(config_adhesin)

    for p in list_of_proteins:
        data = p.standardize(means, std_devs, projection_matrix)
        data = np.array(data)[None, ...]
        p.p_ad = float(model.predict(data, verbose=0))

    return list_of_proteins

