# Author: Hauxu Yu

# A module to summarize the main data processing modules

# Import modules
from . import raw_data_utils as raw
from .params import Params
from .ann_feat_quality import predict_quality
from .feature_grouping import annotate_isotope, annotate_adduct, annotate_in_source_fragment
from .alignment import alignement, sum_aligned_features, output_aligned_features
import pickle
import os
from keras.models import load_model
from .annotation import annotate_features, annotate_rois
import time


def feature_detection(file_name, params, annotation=False):
    """
    Feature detection from a raw LC-MS file (.mzML or .mzXML).

    Parameters
    ----------
    file_name : str
        File name of the raw file.
    parameters : Params object
        The parameters for the workflow.
    """

    # create a MSData object1
    d = raw.MSData()

    d.read_raw_data(file_name, params)  # read raw data

    d.drop_ion_by_int()

    d.find_rois() # find ROIs

    if d.params.cut_roi:
        d.cut_rois()  # cut ROIs

    # sort ROI by m/z, find roi quality by length, find the best MS2
    d.process_rois()

    # predict feature quality
    if d.params.ann_model is None:
        data_path_ann = os.path.join(os.path.dirname(__file__), 'model', "peak_quality_NN.keras")
        d.params.ann_model = load_model(data_path_ann)

    predict_quality(d)

    print("Number of regular ROIs: " + str(len(d.rois)))

    # annotate isotopes, adducts, and in-source fragments
    annotate_isotope(d)

    annotate_in_source_fragment(d)

    annotate_adduct(d)

    if annotation and d.params.msms_library is not None:
        annotate_rois(d)

    # output single file
    if d.params.output_single_file:
        d.output_single_file()

    return d


def process_files(file_names, params):
    """
    A function to process multiple raw files.

    Parameters
    ----------
    file_names : list
        A list of file names of the raw files in .mzML or .mzXML format.
    params : Params object
        The parameters for the workflow.
    """

    feature_list = []

    for file_name in file_names:
        d = feature_detection(file_name, params)
        alignement(feature_list, d)
        print("-----------------------------------")
    
    sum_aligned_features(feature_list)

    # annotation
    if params.msms_library is not None:
        annotate_features(feature_list, params)

    if params.output_aligned_file:
        output_aligned_features(feature_list, file_names, params.project_dir)

    return feature_list


def read_raw_file_to_obj(file_name, params=None):
    """
    Read a raw file to a MSData object.

    Parameters
    ----------
    file_name : str
        The file name of the raw file.
    """

    d = raw.MSData()
    if params is None:
        params = Params()
    d.read_raw_data(file_name, params)
    
    return d


def untargeted_workflow(parameters):
    """
    A function for the untargeted metabolomics workflow.

    Parameters
    ----------
    parameters : Params object
        The parameters for the workflow.
    """

    # Check the folder for creating the project
    if not os.path.exists(parameters.project_dir):
        raise ValueError("The project directory does not exist.")
    
    # Check if raw files exist

