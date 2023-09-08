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
from scipy.interpolate import interp1d, splrep, BSpline
from sklearn.linear_model import LinearRegression
from .alignment import alignement
import pandas as pd


def lcb_workflow(data_dir, sample_type, ion_mode, create_sub_folders=False, output_single_file=False, istd_interval=0.3, 
                 validation=False, params=None, cut_roi=True, mz_method="smooth_spline", rt_method="smooth_spline", match_ms2=False):
    """
    A function to run the data processing workflow by LC-Binbase.

    Parameters
    ----------
    data_dir : str
        Path to the directory containing all data files.
    sample_type : str
        Sample type, "lipidomics" or "metabolomics".
    ion_mode : str
        Ion mode, "positive" or "negative".
    create_sub_folders : bool
        Whether to create three folders for method blank, pooled QC, and sample data. Default is False.
    output_single_file : bool
        Whether to output a single file for each data file. Default is False.
    istd_interval : float
        Retention time interval for selecting internal standards. Default is 0.3.
    validation : bool
        Whether to run validation. Default is False.    
    """

    # create three folders for method blank, pooled QC, and sample data
    if create_sub_folders==True:
        _create_folder(data_dir)
        _ = input("Please put method blank data in the method_blank folder, pooled QC data in the pooled_qc folder, and sample data in the sample folder. Press any key to continue...")
    
    # Load internal standard library
    db = load_internal_standard_library(sample_type, ion_mode)

    print(len(db), "internal standards loaded.")

    # Get all method blank data
    mb_file_names = get_data_name(data_dir, "method_blank")
    # Select internal standards from method blank data as chemical anchors
    istd_selected = select_istd_from_mb(mb_file_names, db)

    # sort the selected internal standards by retention time
    istd_selected = sorted(istd_selected, key=lambda x: x['rt'])

    # split the selected internal standards to training set and validation set
    istd_training = [istd_selected[0]]
    istd_validation = []

    rt_tmp = 0.0
    for i in range(1, len(istd_selected)-1):
        if istd_selected[i]['rt'] - rt_tmp > istd_interval:
            istd_training.append(istd_selected[i])
            rt_tmp = istd_selected[i]['rt']
        else:
            istd_validation.append(istd_selected[i])
    
    istd_training.append(istd_selected[-1])
    
    print("number of internal standards in training set: ", len(istd_training))
    print("number of internal standards in validation set: ", len(istd_validation))

    # Process each file from pooled QC to real samples to method blanks
    qc_file_names = get_data_name(data_dir, "pooled_qc")
    sample_file_names = get_data_name(data_dir, "sample")

    full_file_names = qc_file_names + sample_file_names + mb_file_names

    # create a list for aligned features
    feature_list = []

    # create lists to store the retention time and m/z of the internal standards in the valuation set if run validation
    if validation:
        mz_val_before_correct = []
        rt_val_before_correct = []
        mz_val_after_correct = []
        rt_val_after_correct = []
    
    matched_mzs_all = []
    matched_rts_all = []
        
    for idx_fn, file_name in enumerate(full_file_names):
        print("Processing file", idx_fn+1, "of", len(full_file_names), ":", file_name)
        # feature detection
        d = feat_detection(file_name, params=params, cut_roi=cut_roi, output_single_file=output_single_file)

        # find the selected anchors (internal standards) in the file by m/z, rt, intensity, and MS/MS spectra (if available)
        matched_mzs, matched_rts, _ = find_itsd_from_rois(d, istd_training)

        # if run validation, find the corrected m/z and retention time for the internal standards in the validation set
        if validation:
            matched_mzs_val, matched_rts_val, matched_roi_idx = find_itsd_from_rois(d, istd_validation)
            if None in matched_roi_idx:
                k = np.where(matched_roi_idx == None)
                return k, istd_validation
            else:
                matched_mzs_val = d.rois_mz_seq[matched_roi_idx]
                matched_rts_val = d.rois_rt_seq[matched_roi_idx]

            mz_val_before_correct.append(matched_mzs_val)
            rt_val_before_correct.append(matched_rts_val)

        matched_mzs_all.append(matched_mzs)
        matched_rts_all.append(matched_rts)
        # correct retention time
        correct_retention_time(d, matched_rts, istd_training, method=rt_method)

        # correct m/z
        correct_mz(d, matched_mzs, istd_training, method=mz_method)

        # if run validation, find the corrected m/z and retention time for the internal standards in the validation set
        if validation:
            matched_mzs_val = d.rois_mz_seq[matched_roi_idx]
            matched_rts_val = d.rois_rt_seq[matched_roi_idx]
            mz_val_after_correct.append(matched_mzs_val)
            rt_val_after_correct.append(matched_rts_val)
            
    #     # alignment
    #     print('Running alignment on: ', file_name)
    #     alignement(feature_list, d)

    # # choose the best MS2 for each feature
    # for feat in feature_list:
    #     feat.choose_best_ms2()
    
    # output the validation results if run validation
    if validation:
        before_correct_result = validate_istd(istd_validation, mz_val_before_correct, rt_val_before_correct)
        after_correct_result = validate_istd(istd_validation, mz_val_after_correct, rt_val_after_correct)

        # create a pandas dataframe to store the validation results
        validation_result = pd.DataFrame({
            "name": [i['name'] for i in istd_validation],
            "blank_int": [i['blank_int'] for i in istd_validation],
            "mz_MAE_before_correct": before_correct_result['mz_MAE'],
            "rt_MAE_before_correct": before_correct_result['rt_MAE'],
            "mz_MAXAE_before_correct": before_correct_result['mz_MAXAE'],
            "rt_MAXAE_before_correct": before_correct_result['rt_MAXAE'],
            "mz_MAE_after_correct": after_correct_result['mz_MAE'],
            "rt_MAE_after_correct": after_correct_result['rt_MAE'],
            "mz_MAXAE_after_correct": after_correct_result['mz_MAXAE'],
            "rt_MAXAE_after_correct": after_correct_result['rt_MAXAE']
        })

        # save the validation results to a csv file
        validation_result.to_csv(os.path.join(data_dir, "validation_result.csv"), index=False)
    
    return feature_list, istd_training, istd_validation, mz_val_before_correct, rt_val_before_correct, mz_val_after_correct, rt_val_after_correct, before_correct_result, after_correct_result, matched_mzs_all, matched_rts_all


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

    if sample_type.lower() == "lipidomics":
        if ion_mode.lower() == "positive":
            fn = "lipid_istd_pos.json"
        elif ion_mode.lower() == "negative":
            fn = "lipid_istd_neg.json"
        
    elif sample_type.lower() == "metabolomics":
        if ion_mode.lower() == "positive":
            fn = "metabolite_istd_pos.json"
        elif ion_mode.lower() == "negative":
            fn = "metabolite_istd_neg.json"
    data_path = os.path.join(os.path.dirname(__file__), 'data', fn)

    with open(data_path, 'r') as f:
        tmp = json.load(f)

    db = [i for i in tmp if i['rt'] == i['rt']]

    return db

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


    if sub_folder=="method_blank" and len(files) == 0:
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
    Detected in all blank files by matching the m/z, retention time, and MS/MS spectra (if available)

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

    target_mzs = [i['preferred_mz'] for i in db]
    target_rts = [i['rt'] for i in db]
    target_rt_ranges = []
    for i in range(len(db)):
        if target_rts[i] == target_rts[i]:
            target_rt_ranges.append([target_rts[i]-rt_tol, target_rts[i]+rt_tol])
        else:
            target_rt_ranges.append([0, np.inf])
    
    int_multi_files = []

    for fn in mb_file_names:

        # read raw file
        d = read_raw_file_to_obj(fn)

        # estimate intensity tolerance
        params = Params()
        params.estimate_params(d, estimate_mz_tol=False, estimate_int_tol=True)

        # record whether the internal standard is detected in the file
        matched = [False] * len(target_mzs)
        int_single_file = [0.0] * len(target_mzs)
        
        for i in range(len(target_mzs)):
            # skip if the internal standard retention time was not recorded
            if target_rts[i] != target_rts[i]:
                continue

            eic_int= d.get_eic_data(target_mz=target_mzs[i], mz_tol=0.01, rt_range=target_rt_ranges[i])[1]
            
            # the largest intensity is above threshold (i.e. distinguishable from noise)
            if np.max(eic_int) > params.int_tol:
                # if there is a local maximum that is higher than 50% of the highest intensity
                # change all values in eic_int that are lower than 50% of the highest intensity to 0
                eic_int[eic_int < np.max(eic_int) * 0.5] = 0
                if _single_max(eic_int):
                    matched[i] = True
                    int_single_file[i] = np.max(eic_int)

        matched_multi_files.append(matched)
        int_multi_files.append(int_single_file)

    # get the internal standards that are detected in all blank files
    matched_multi_files = np.array(matched_multi_files)
    id = np.where(matched_multi_files.sum(axis=0) == len(mb_file_names))[0]

    for i in id:
        db[i]['blank_int'] = np.mean([int_multi_files[j][i] for j in range(len(mb_file_names))])

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


def find_itsd_from_rois(d, istds, mz_tol=0.012, rt_tol=0.2, dp_tol=0.7, match_ms2=False):
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

    matched_mzs = []
    matched_rts = []
    matched_roi_idx = []

    # ROIs have been sorted by m/z

    for istd in istds:
        mz = istd['preferred_mz']
        rt = istd['rt']

        matched_idx = np.where(np.logical_and(np.abs(d.rois_mz_seq-mz) < 0.012, np.abs(d.rois_rt_seq-rt) < rt_tol))[0]

        if len(matched_idx) == 0:
            matched_mzs.append(None)
            matched_rts.append(None)
            matched_roi_idx.append(None)
            continue
        
        # if match by MS/MS spectra
        if match_ms2:
            pass
        
        if len(matched_idx) == 1:
            matched_mzs.append(d.rois_mz_seq[matched_idx[0]])
            matched_rts.append(d.rois_rt_seq[matched_idx[0]])
            matched_roi_idx.append(matched_idx[0])
            continue   
        
        # if there are multiple ROIs matched, choose the one with reasonable intensity (within 2 fold change) and closest retention time
        if len(matched_idx) > 1:
            roi_int = np.array([d.rois[i].peak_height for i in matched_idx])
            roi_rt = np.array([d.rois_rt_seq[i] for i in matched_idx])
            idx = np.where(np.logical_and(roi_int/istd['blank_int']>0.5, roi_int/istd['blank_int']<2))[0]
            # if there is no ROI with reasonable intensity, choose the one with the cloest retention time
            if len(idx) == 0:
                matched_mzs.append(d.rois_mz_seq[matched_idx[np.argmin(np.abs(roi_rt-rt))]])
                matched_rts.append(d.rois_rt_seq[matched_idx[np.argmin(np.abs(roi_rt-rt))]])
                matched_roi_idx.append(matched_idx[np.argmin(np.abs(roi_rt-rt))])
            # if there is only one ROI with reasonable intensity, choose it
            if len(idx) == 1:
                matched_mzs.append(d.rois_mz_seq[matched_idx[idx[0]]])
                matched_rts.append(d.rois_rt_seq[matched_idx[idx[0]]])
                matched_roi_idx.append(matched_idx[idx[0]])
            # if there are multiple ROIs with reasonable intensity, choose the one with the cloest retention time
            if len(idx) > 1:
                matched_idx = matched_idx[idx]
                roi_rt = np.array([d.rois_rt_seq[i] for i in matched_idx])
                try:
                    matched_mzs.append(d.rois_mz_seq[matched_idx[np.argmin(np.abs(roi_rt-rt))]])
                    matched_rts.append(d.rois_rt_seq[matched_idx[np.argmin(np.abs(roi_rt-rt))]])
                    matched_roi_idx.append(matched_idx[np.argmin(np.abs(roi_rt-rt))])
                except:
                    print(matched_idx, idx)
    
    matched_mzs = np.array(matched_mzs, dtype=np.float64)
    matched_rts = np.array(matched_rts, dtype=np.float64)
    matched_roi_idx = np.array(matched_roi_idx)

    return matched_mzs, matched_rts, matched_roi_idx


def correct_retention_time(d, matched_rts, istd_training, method="smooth_spline"):
    """
    Correct retention time for each ROI using the matched internal standards.
    
    Parameters
    ----------------------------------------------------------
    d: MSData object
        Raw data object after feature detection.
    matched_rts: list
        List of retention times of matched internal standards.
    istd_selected: list
        List of selected internal standards.
    method: str
        "smooth_spline": smooth spline interpolation
        "linear_interp": linear interpolation
        "linear_regression": linear regression
        "polynomial_regression": polynomial regression
        
    Utilities
    ----------------------------------------------------------
    Change the retention time of each ROI in d
    """

    # data retention time
    x = matched_rts
    # reference retention time
    y = np.array([i['rt'] for i in istd_training], dtype=np.float64)

    if method=="smooth_spline":
        try:
            tck = splrep(x, y, s=len(x)*0.005)
            d.rois_rt_seq = BSpline(*tck)(d.rois_rt_seq)
        except:
            print(x, y)

    if method=="linear_interp":
        try:
            f = interp1d(x, y)
            d.rois_rt_seq = f(d.rois_rt_seq)
        except:
            print(x, y)

    if method=="linear_regression":
        reg = LinearRegression().fit(x.reshape(-1, 1), y)
        d.rois_rt_seq = reg.predict(d.rois_rt_seq.reshape(-1, 1))

    if method=="polynomial_regression":
        reg = LinearRegression().fit(np.vstack([x, x**2, x**3]).T, y)
        d.rois_rt_seq = reg.predict(np.vstack([d.rois_rt_seq, d.rois_rt_seq**2, d.rois_rt_seq**3]).T)

    # Change the retention time of each ROI in d
    for i in range(len(d.rois)):
        d.rois[i].rt = d.rois_rt_seq[i]     


def correct_mz(d, matched_mzs, istd_training, method="smooth_spline"):
    """
    Correct m/z for each ROI using the matched internal standards.

    Parameters
    ----------------------------------------------------------
    d: MSData object
        Raw data object after feature detection.
    matched_mzs: list
        List of m/z values of matched internal standards.
    istd_selected: list
        List of selected internal standards.
    method: str 
        "smooth_spline": smooth spline interpolation
        "linear_interp": linear interpolation
        "linear_regression": linear regression
        "polynomial_regression": polynomial regression

    Utilities
    ----------------------------------------------------------
    Change the m/z of each ROI in d
    """

    # data m/z
    x = matched_mzs
    # reference m/z
    y = np.array([i['preferred_mz'] for i in istd_training], dtype=np.float64)

    # sort x and y by x
    idx = np.argsort(x)
    x = x[idx]
    y = y[idx]

    if method=="smooth_spline":
        try:
            tck = splrep(x, y, s=len(x)*0.0001)
            d.rois_mz_seq = BSpline(*tck)(d.rois_mz_seq)
        except:
            print(x, y)
    
    if method=="linear_interp":
        try:
            f = interp1d(x, y)
            d.rois_mz_seq = f(d.rois_mz_seq)
        except:
            print(x, y)

    if method=="linear_regression":
        reg = LinearRegression().fit(x.reshape(-1, 1), y)
        d.rois_mz_seq = reg.predict(d.rois_mz_seq.reshape(-1, 1))

    if method=="polynomial_regression":
        reg = LinearRegression().fit(np.vstack([x, x**2, x**3]).T, y)
        d.rois_mz_seq = reg.predict(np.vstack([d.rois_mz_seq, d.rois_mz_seq**2, d.rois_mz_seq**3]).T)

    # Change the m/z of each ROI in d
    for i in range(len(d.rois)):
        d.rois[i].mz = d.rois_mz_seq[i]


def validate_istd(istd_validation, mz_val, rt_val):
    """
    A function to validate the internal standards and generate a report.
    Evaluation criteria:
    1. median absolute error (MAE) of m/z
    2. median absolute error (MAE) of retention time
    3. maximum absolute error (MAXAE) of m/z
    4. maximum absolute error (MAXAE) of retention time

    Parameters
    ----------
    istd_validation : list
        List of internal standards in the validation set.
    mz_val : list
        List of corrected m/z values of the internal standards in the validation set.
    rt_val : list
        List of corrected retention times of the internal standards in the validation set.
    """

    target_mzs = np.array([i['preferred_mz'] for i in istd_validation])
    target_rts = np.array([i['rt'] for i in istd_validation])

    mz_MAE = []
    rt_MAE = []
    mz_MAXAE = []
    rt_MAXAE = []

    # to make one row for each internal standard
    mz_val = np.array(mz_val).T
    rt_val = np.array(rt_val).T

    for i in range(len(target_mzs)):
        mz_MAE.append(np.median(np.abs(mz_val[i] - target_mzs[i])))
        rt_MAE.append(np.median(np.abs(rt_val[i] - target_rts[i])))
        mz_MAXAE.append(np.max(np.abs(mz_val[i] - target_mzs[i])))
        rt_MAXAE.append(np.max(np.abs(rt_val[i] - target_rts[i])))

    result = {
        "mz_MAE": mz_MAE,
        "rt_MAE": rt_MAE,
        "mz_MAXAE": mz_MAXAE,
        "rt_MAXAE": rt_MAXAE
    }

    return result

    