# Author: Hauxu Yu

# A module to support the data processing workflow by LC-Binbase

# Workflow
# 1. Load internal standard library with no retention time and MS/MS spectra
# 2. Get all method blank data
# 3. Run peak picking of method blank data
# 4. Targeted search for internal standards in method blank data
# 5. For observed internal standards, get retention time and MS/MS spectra, and add to the library (project-specific)
# 6. Generate aligned feature table for method blank data
# 7. Get all pooled quality control sample data
# 8. Run peak picking of pooled QC and correct retention time using internal standards
# 8. Align pooled QCs
# 9. Get all sample data
# 10. Run peak picking of sample data and correct retention time using internal standards
# 11. Align sample data

# Import modules
import json
import os
from . import read_raw_file_to_obj
from .params import Params
import numpy as np

def lcb_workflow(data_dir, sample_type, ion_mode, create_sub_folders=False):

    # create three folders for method blank, pooled QC, and sample data
    if create_sub_folders==True:
        _create_folder(data_dir)
    
    # Load internal standard library
    db = load_internal_standard_library(sample_type, ion_mode)

    # Get all method blank data
    mb_file_names = get_method_blank_data(data_dir)
    select_istd_from_mb(mb_file_names, db)







def load_internal_standard_library(sample_type, ion_mode):
    """
    Load the pre-defined internal standard library for the given sample type and ion mode.

    Parameters
    ----------
    sample_type : str
        Sample type, "lipidomics" or "metabolomics".
    ion_mode : str
        Ion mode, "positive" or "negative".    
    """

    if sample_type == "lipidomics":
        if ion_mode == "positive":
            fn = "lipid_istd_pos.json"
        elif ion_mode == "negative":
            fn = "lipid_istd_neg.json"
        
    elif sample_type == "metabolomics":
        if ion_mode == "positive":
            fn = "metabolite_istd_pos.json"
        elif ion_mode == "negative":
            fn = "metabolite_istd_neg.json"
    data_path = os.path.join(os.path.dirname(__file__), 'data', fn)

    with open(data_path, 'r') as f:
        return json.load(f)


def get_method_blank_data(data_dir):
    """
    Get all method blank data.

    Parameters
    ----------
    data_dir : str
        Path to the directory containing all data files.
    """

    subpath = os.path.join(data_dir, "method_blank")

    # Get all files in the directory ending with ".mzML" or ".mzXML"
    files = [f for f in os.listdir(subpath) if f.endswith(".mzML") or f.endswith(".mzXML")]


    if len(files) == 0:
        raise Exception("No data files found in the method blank folder.")
    
    files = [os.path.join(subpath, f) for f in files]
    
    return files
    

def _create_folder(data_dir):
    """
    Create three folders for method blank, pooled QC, and sample data.
    """

    os.makedirs(os.path.join(data_dir, "method_blank"))
    os.makedirs(os.path.join(data_dir, "pooled_qc"))
    os.makedirs(os.path.join(data_dir, "sample"))


def select_istd_from_mb(mb_file_names, db)
    """
    Select internal standards from method blank data as chemical anchors.

    The internal standards should be selected based on the following criteria:
    1. Detected in all blank files
    2. No other isobaric compounds with intensity higher than 20% of the highest intensity
    3. Have MS/MS spectra if the data acquisiton mode is data-dependent acquisition (DDA)

    Parameters
    ----------
    mb_file_names : list
        List of file names of method blank data.
    db : dict
        Internal standard library.
    """

    matched_multi_files = []
    target_mzs = []

    for istd in db:
        target_mzs.append(istd['common_adducts_mz'][istd['common_adducts'].index(istd["preferred_adduct"])])

    for fn in mb_file_names:
        d = read_raw_file_to_obj(fn)
        params = Params()
        params.estimate_params(d, estimate_mz_tol=False, estimate_cycle_time=False, estimate_int_tol=True)
        int_tol = params.int_tol
        matched = [False] * len(target_mzs)
        
        for i, mz in enumerate(target_mzs):
            _, eic_int = d.get_eic_data(mz, 0.01)
            
            # the largest intensity is above threshold (i.e. distinguishable from noise)
            if np.max(eic_int) > params.int_tol:
                # if there is a local maximum that is higher than 20% of the highest intensity
                # change all values in eic_int that are lower than 20% of the highest intensity to 0
                eic_int[eic_int < np.max(eic_int) * 0.2] = 0
                if _single_max(eic_int):
                    matched[i] = True
        matched_multi_files.append(matched)

    # get the internal standards that are detected in all blank files
    matched_multi_files = np.array(matched_multi_files)
    matched_multi_files = matched_multi_files.sum(axis=0)
    matched_multi_files = matched_multi_files == len(mb_file_names)
    target_mzs = np.array(target_mzs)
    target_mzs = target_mzs[matched_multi_files]

    # get the internal standards that have MS/MS spectra
    target_mzs = target_mzs[db['msms'] == True]
    
                        


def _single_max(a):
    """
    Return True if there is only one local maximum in the array.

    Parameters
    ----------
    a : numpy.ndarray
        Input array.    
    """
    counter = 0
    for i in range(1, len(a)-1):
        if a[i] == 0 and a[i+1] > 0:
            counter += 1
        if counter > 1:
            return False
    return True

                



    