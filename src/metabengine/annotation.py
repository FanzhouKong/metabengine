# Author: Hauxu Yu

# A module to annotate metabolites based on their MS/MS spectra

# Import modules
import os
from ms_entropy import read_one_spectrum
import pickle

def laod_msms_db(path):
    """
    A function to load the MS/MS database in MSP format or pickle format.

    Parameters
    ----------
    path : str
        The path to the MS/MS database in MSP format.    
    """

    # get extension of path
    ext = os.path.splitext(path)[1]

    if ext.lower() == '.msp':
        db =[]
        for a in read_one_spectrum('D:/Database/NIST23/nist23_negative.MSP'):
            db.append(a)
        return db
    
    elif ext.lower() == '.pickle':
        db = pickle.load(open(path, 'rb'))


def annotate_features(feature_list, db_path):
    """
    A function to annotate features based on their MS/MS spectra and a MS/MS database.

    Parameters
    ----------
    feature_list : list
        A list of features.
    db_path : str
        The path to the MS/MS database in MSP or pickle format.   
    """

    # load the MS/MS database
    db = laod_msms_db(db_path)

    # annotate features
    for f in feature_list:
        f.annotate(db)
