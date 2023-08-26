# Author: Hauxu Yu

# A module to support the data processing workflow by LC-Binbase

# Workflow
# 1. Load internal standard library with no retention time and MS/MS spectra
# 2. Get all method blank data
# 3. Targeted search for internal standards in method blank data
# 4. For each file (Pooled QC to real samples to method blanks), run
#    Peak picking, find internal standards, run retention time correction, run m/z correction, run alignment

# Import modules
import json
import os
from . import read_raw_file_to_obj, feat_detection
from .params import Params
import numpy as np

def lcb_workflow(data_dir, sample_type, ion_mode, create_sub_folders=False, output_single_file=False):

    # create three folders for method blank, pooled QC, and sample data
    if create_sub_folders==True:
        _create_folder(data_dir)
        _ = input("Please put method blank data in the method_blank folder, pooled QC data in the pooled_qc folder, and sample data in the sample folder. Press Enter to continue...")
    
    # Load internal standard library
    db = load_internal_standard_library(sample_type, ion_mode)

    # Get all method blank data
    mb_file_names = get_data_name(data_dir, "method_blank")
    # Select internal standards from method blank data as chemical anchors
    istd_selected = select_istd_from_mb(mb_file_names, db)

    # Process each file from pooled QC to real samples to method blanks
    qc_file_names = get_data_name(data_dir, "pooled_qc")
    sample_file_names = get_data_name(data_dir, "sample")

    full_file_names = qc_file_names + sample_file_names + mb_file_names

    # create a list for aligned features
    feature_list = []

    for file_name in full_file_names:
        # feature detection
        d = feat_detection(file_name, output_single_file=output_single_file)

        # find the selected anchors (internal standards) in the file by m/z, rt, intensity, and MS/MS spectra (if available)
        matched_mz, matched_rt = find_itsd_from_rois(d, istd_selected)

        





        print('Running alignment on: ', file_name)
        # alignement(feature_list, d)
    
    # choose the best MS2 for each feature
    for feat in feature_list:
        feat.choose_best_ms2()
    
    return feature_list









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


def get_data_name(data_dir, sub_folder=None):
    """
    Get all method blank data.

    Parameters
    ----------
    data_dir : str
        Path to the directory containing all data files.
    sub_folder : str
        Name of the sub folder containing method blank data.
    """

    if sub_folder:
        subpath = os.path.join(data_dir, sub_folder)

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


def select_istd_from_mb(mb_file_names, db, mz_tol=0.01, rt_tol=1.0, match_ms2=False):
    """
    Select internal standards from method blank data as chemical anchors.

    The internal standards should be selected based on the following criteria:
    1. Detected in all blank files by matching the m/z, retention time, and MS/MS spectra (if available)
    2. No isobaric interference within the retention time window (>20% intensity)

    Parameters
    ----------
    mb_file_names : list
        List of file names of method blank data.
    db : dict
        Internal standard library.
    match_ms2 : bool
        Whether to match MS/MS spectra. Default is False.
    
    Returns
    ----------
    A list of internal standards that are detected in all blank files.
    An extra key will be added to record the average intensity of the internal standard in the blank files.
    """

    matched_multi_files = []
    target_mzs = []
    # target_rts = []

    target_mzs = [i['preferred_mz'] for i in db]
    target_rts = [i['rt'] for i in db]
    target_rt_ranges = []
    for i in range(len(db)):
        if target_rts[i] == target_rts[i]:
            target_rt_ranges.append([target_rts[i]-rt_tol, target_rts[i]+rt_tol])
        else:
            target_rt_ranges.append([0, np.inf])
    
    failed = []

    for fn in mb_file_names:

        # read raw file
        d = read_raw_file_to_obj(fn)

        # estimate intensity tolerance
        params = Params()
        params.estimate_params(d, estimate_mz_tol=False, estimate_cycle_time=False, estimate_int_tol=True)

        # record whether the internal standard is detected in the file
        matched = [False] * len(target_mzs)
        
        for i in range(len(target_mzs)):

            if i in failed:
                continue

            _, eic_int = d.get_eic_data(target_mz=target_mzs[i], mz_tol=0.01, rt_range=target_rt_ranges[i])
            
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
    id = np.where(matched_multi_files.sum(axis=0) == len(mb_file_names))[0]
    return [db[i] for i in id] 
                        

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


def find_itsd_from_rois(d, istds, mz_tol=0.01, rt_tol=1.0, dp_tol=0.7, match_ms2=False):
    """
    Find internal standards from ROIs for retention time and m/z correction.

    Parameters
    ----------------------------------------------------------
    d: MSData object
        Raw data object after feature detection.
    istds: list
        List of internal standards.
    mz_tol: float
        m/z tolerance.
    rt_tol: float
        Retention time tolerance.
    dp_tol: float
        Dot product tolerance.
    match_ms2: bool
        Whether to match MS/MS spectra. Default is False.

    Returns
    ----------------------------------------------------------
    Two lists:
    1. A list of m/z values of matched internal standards.
    2. A list of retention times of matched internal standards.
    """

    # ROIs have been sorted by m/z

    for istd in istds:
        mz = istd['preferred_mz']
        rt = istd['rt']

        matched_idx = np.where(np.logical_and(np.abs(d.rois_mz_seq-mz) < mz_tol, np.abs(d.rois_rt_seq-rt) < rt_tol))[0]
    
           



    