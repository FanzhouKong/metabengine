# Author: Hauxu Yu

# A module to summarize the main data processing modules

# Import modules
from . import raw_data_utils as raw
from .params import Params, Cross_file_params
from .ann_feat_quality import predict_quality
from .feature_grouping import annotate_isotope
from .alignment import alignement

def feat_detection(file_name, pred_quality_NN=False, annotate_isotope=False, output_single_file=False, path=None):

    d = raw.MSData()
    params = Params()
    d.read_raw_data(file_name, params)
    params.estimate_params(d)
    print(params.int_tol)
    d.drop_ion_by_int(params)
    d.find_rois(params)
    d.cut_rois(params)
    d.process_rois(params)

    if pred_quality_NN:
        predict_quality(d)
    if annotate_isotope:
        annotate_isotope(d, params)

    if output_single_file:
        d.output_roi_report(path)

    return d


def process_files(file_names):

    feature_list = []

    cross_file_params = Cross_file_params()

    for file_name in file_names:
        d = feat_detection(file_name)
        alignement(feature_list, d, cross_file_params)
    
    return feature_list
