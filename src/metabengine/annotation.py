# Author: Hauxu Yu

# A module to annotate metabolites based on their MS/MS spectra

# Import modules
import os
from ms_entropy import read_one_spectrum, FlashEntropySearch
import pickle
import numpy as np
import json

def load_msms_db(path):
    """
    A function to load the MS/MS database in MSP format or pickle format.

    Parameters
    ----------
    path : str
        The path to the MS/MS database in MSP format.    
    """

    print("Loading MS/MS database...")
    # get extension of path
    ext = os.path.splitext(path)[1]

    if ext.lower() == '.msp':
        db =[]
        for a in read_one_spectrum(path):
            db.append(a)
        entropy_search = FlashEntropySearch()
        entropy_search.build_index(db)
        print("MS/MS database loaded.")
        return entropy_search
    
    elif ext.lower() == '.pkl':
        entropy_search = pickle.load(open(path, 'rb'))
        print("MS/MS database loaded.")
        return entropy_search
    
    elif ext.lower() == '.json':
        db = json.load(open(path, 'r'))
        entropy_search = FlashEntropySearch()
        entropy_search.build_index(db)
        print("MS/MS database loaded.")
        return entropy_search


def annotate_features(feature_list, params):
    """
    A function to annotate features based on their MS/MS spectra and a MS/MS database.

    Parameters
    ----------
    feature_list : list
        A list of features.
    params : Params object
        The parameters for the workflow.
    """

    # load the MS/MS database
    entropy_search = load_msms_db(params.msms_library)

    for f in feature_list:

        f.annotation = None
        f.similarity = None
        f.matched_peak_number = None
        f.smiles = None
        f.inchikey = None
        f.matched_precursor_mz = None
        f.matched_peaks = None
        f.formula = None
        
        if f.best_ms2 is not None:
            peaks = entropy_search.clean_spectrum_for_search(f.mz, f.best_ms2.peaks)
            entropy_similarity, matched_peaks_number = entropy_search.identity_search(precursor_mz=f.mz, peaks=peaks, ms1_tolerance_in_da=params.mz_tol_ms1, 
                                                                                      ms2_tolerance_in_da=params.mz_tol_ms2, output_matched_peak_number=True)
            
            idx = np.argmax(entropy_similarity)
            if entropy_similarity[idx] > params.ms2_sim_tol:
                matched = entropy_search[np.argmax(entropy_similarity)]
                matched = {k.lower():v for k,v in matched.items()}
                f.annotation = matched['name']
                f.similarity = entropy_similarity[idx]
                f.matched_peak_number = matched_peaks_number[idx]
                f.smiles = matched['smiles'] if 'smiles' in matched else None
                f.inchikey = matched['inchikey'] if 'inchikey' in matched else None
                f.matched_precursor_mz = matched['precursor_mz']
                f.matched_peaks = matched['peaks']
                f.formula = matched['formula'] if 'formula' in matched else None


def annotate_rois(d):
    """
    A function to annotate rois based on their MS/MS spectra and a MS/MS database.

    Parameters
    ----------
    d : MSData object
        MS data.
    """

    # load the MS/MS database
    entropy_search = load_msms_db(d.params.msms_library)
    

    for f in d.rois:
        f.annotation = None
        f.similarity = None
        f.matched_peak_number = None
        f.smiles = None
        f.inchikey = None
        f.matched_precursor_mz = None
        f.matched_peaks = None
        f.formula = None

        if f.best_ms2 is not None:
            peaks = entropy_search.clean_spectrum_for_search(f.mz, f.best_ms2.peaks)
            entropy_similarity, matched_peaks_number = entropy_search.identity_search(precursor_mz=f.mz, peaks=peaks, ms1_tolerance_in_da=d.params.mz_tol_ms1, 
                                                                                      ms2_tolerance_in_da=d.params.mz_tol_ms2, output_matched_peak_number=True)
            
            idx = np.argmax(entropy_similarity)
            if entropy_similarity[idx] > d.params.ms2_sim_tol:
                matched = entropy_search[np.argmax(entropy_similarity)]
                matched = {k.lower():v for k,v in matched.items()}
                f.annotation = matched['name']
                f.similarity = entropy_similarity[idx]
                f.matched_peak_number = matched_peaks_number[idx]
                f.smiles = matched['smiles'] if 'smiles' in matched else None
                f.inchikey = matched['inchikey'] if 'inchikey' in matched else None
                f.matched_precursor_mz = matched['precursor_mz']
                f.matched_peaks = matched['peaks']
                f.formula = matched['formula'] if 'formula' in matched else None


def has_chlorine(iso):
    # to be constructed
    pass

def has_bromine(iso):
    # to be constructed
    pass