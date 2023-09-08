# Author: Hauxu Yu

# A module to summarize the main data processing modules

# Import modules
from . import raw_data_utils as raw
from .params import Params
from .ann_feat_quality import predict_quality
from .feature_grouping import annotate_isotope
from .alignment import alignement
import pickle

def feat_detection(file_name, params=None, cut_roi=True, pred_quality_NN=False, anno_isotope=False, output_single_file=False, discard_short_roi=False):
    """
    A function to detect features from a raw file.

    Parameters
    ----------
    file_name : str
        The file name of the raw file.
    pred_quality_NN : bool
        Whether to predict the quality of the features. The default is False.
    anno_isotope : bool
        Whether to annotate the isotopes. The default is False.
    output_single_file : bool
        Whether to output the feature report to a single file. The default is False.
    path : str
        The path to the output file. The default is None.    
    """

    # create a MSData object1
    d = raw.MSData()

    if params is None:
        estimate_params = True
        # create a Params object
        params = Params()
    else:
        params = params
        estimate_params = False

    d.read_raw_data(file_name, params)  # read raw data

    if estimate_params:
        params.estimate_params(d, estimate_mz_tol=False, estimate_int_tol=True) # estimate intensity tolerance for noise filtering
    
    d.drop_ion_by_int(params)   # drop ions by intensity

    if estimate_params:
        params.estimate_params(d, estimate_mz_tol=True, estimate_int_tol=False) # estimate m/z tolerance for peak picking
    d.params = params   # assign the params object to the MSData object

    d.find_rois(params) # find ROIs

    if cut_roi:
        d.cut_rois(params)  # cut ROIs

    # sort ROI by m/z, find roi quality by length, find the best MS2
    d.process_rois(params, discard_short_roi=discard_short_roi)

    if pred_quality_NN:
        predict_quality(d)
    if anno_isotope:
        annotate_isotope(d, params)

    if output_single_file:
        d.output_roi_report(output_single_file)

    return d


def process_files(file_names, params=None, output_single_file=False, discard_short_roi=False, output_aligned_file=None):

    feature_list = []

    for file_name in file_names:
        d = feat_detection(file_name, params=params, output_single_file=output_single_file, discard_short_roi=discard_short_roi)
        print('Running alignment on: ', file_name)
        alignement(feature_list, d)
    
    # choose the best MS2 for each feature
    for feat in feature_list:
        feat.choose_best_ms2()
    
    try:
        if output_aligned_file:
            with open(output_aligned_file, 'wb') as f:
                pickle.dump(feature_list, f)
    except:
        print('Unable to save the aligned feature list.')

    return feature_list


def read_raw_file_to_obj(file_name, estimate_param=False):
    """
    Read a raw file and return a MSData object.

    Parameters
    ----------
    file_name : str
        The file name of the raw file.
    estimate_param : bool
        Whether to estimate the parameters of the MSData object. The default is False.
    """

    d = raw.MSData()
    params = Params()
    d.read_raw_data(file_name, params)

    if estimate_param:
        params.estimate_params(d)
        d.params = params
    
    return d