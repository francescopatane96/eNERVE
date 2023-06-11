
""" Runs subcellular localization prediction """

from Utils import *
from Protein import *
import joblib

def euloc(list_of_proteins, working_dir, NERVE_dir, loclimit) -> list:
    """Subcellular localization module"""
    
    logging.basicConfig(filename = os.path.join(working_dir, 'logfile.log'),
                        filemode = 'a',
                        level = logging.DEBUG,
                        force = True)

    model_dir = os.path.join(NERVE_dir, 'models/')
    model = joblib.load(os.path.join(model_dir, 'loc.joblib'))


    for p in list_of_proteins:
        data = p.model_raw_data
        data = np.array(data)[None, ...]
        p.reliability_out = model.predict_proba(data)[0][1]
        p.localization = 'out' if p.reliability_out >= loclimit else 'in'


    return list_of_proteins


