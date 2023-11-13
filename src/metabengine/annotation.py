# Author: Hauxu Yu

# A module to annotate metabolites based on their MS/MS spectra

# Import modules
import os
from ms_entropy import read_one_spectrum, FlashEntropySearch
import pickle
import numpy as np

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
        entropy_search = FlashEntropySearch()
        entropy_search.build_index(db)

        return entropy_search
    
    elif ext.lower() == '.pickle':
        entropy_search = pickle.load(open(path, 'rb'))


def annotate_features(feature_list, db_path, params):
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
    entropy_search = laod_msms_db(db_path)

    # annotate features
    for f in feature_list:
        entropy_similarity = entropy_search.identity_search(precursor_mz=f.mz, peaks=f.best_ms2.peaks, ms1_tolerance_in_da=params.mz_tol_ms1, ms2_tolerance_in_da=params.mz_tol_ms2)
        idx = np.argmax(entropy_similarity)
        if entropy_similarity[idx] > params.ms2_sim_tol:
            matched = entropy_search[np.argmax(entropy_similarity)]
            f.annotation = matched.name
            f.similarity = entropy_similarity[idx]
            # f.matched_peak_number = len(matched.peaks)
            f.smiles = matched.smiles
            f.inchikey = matched.inchikey


def has_chlorine(iso):
    pass

def has_bromine(iso):
    pass
