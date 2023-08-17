# Author: Hauxu Yu

# A module to summarize the main data processing modules

# Import modules
from . import raw_data_utils as raw
from .params import Params, Cross_file_params
from .ann_feat_quality import predict_quality
from .feature_grouping import annotate_isotope
from .alignment import alignement

def feat_detection(file_name, pred_quality_NN=False, anno_isotope=False, output_single_file=False, path=None):

    d = raw.MSData()
    params = Params()
    d.read_raw_data(file_name, params)
    params.estimate_params(d, estimate_mz_tol=False, estimate_int_tol=True)
    d.drop_ion_by_int(params)
    params.estimate_params(d, estimate_mz_tol=True, estimate_int_tol=False)
    d.params = params
    d.find_rois(params)
    d.cut_rois(params)

    # sort ROI by m/z, find roi quality by length, find the best MS2
    d.process_rois(params)

    if pred_quality_NN:
        predict_quality(d)
    if anno_isotope:
        annotate_isotope(d, params)

    if output_single_file:
        d.output_roi_report(path)

    return d


def process_files(file_names):

    feature_list = []

    cross_file_params = Cross_file_params()

    for file_name in file_names:
        d = feat_detection(file_name)
        print('Running alignment on: ', file_name)
        alignement(feature_list, d, cross_file_params)
    
    # choose the best MS2 for each feature
    for feat in feature_list:
        feat.choose_best_ms2()
    
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