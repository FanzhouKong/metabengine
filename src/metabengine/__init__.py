# Author: Hauxu Yu

# A module to summarize the main data processing modules

# Import modules
from . import raw_data_utils as raw
from .params import Params
import pandas as pd
import numpy as np
from .ann_feat_quality import predict_quality
from .feature_grouping import annotate_isotope, annotate_adduct, annotate_in_source_fragment
from .alignment import alignement, sum_aligned_features, output_aligned_features
import os
from keras.models import load_model
from .annotation import annotate_features, annotate_rois
import pickle

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

    annotated = False
    if annotation and d.params.msms_library is not None:
        annotate_rois(d)
        annotated = True

    # output single file
    if d.params.output_single_file:
        d.output_single_file(annotated)

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
        print("Processing file: " + file_name)
        d = feature_detection(file_name, params)
        d.rois = [roi for roi in d.rois if roi.length >= 5 or roi.best_ms2 is not None]
        print("Number of regular ROIs after discarding short ROIs: " + str(len(d.rois)))
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
    
    # Check if the three file folders have been created
    if not os.path.exists(os.path.join(parameters.project_dir, "qc")):
        os.makedirs(os.path.join(parameters.project_dir, "qc"))
    if not os.path.exists(os.path.join(parameters.project_dir, "blank")):
        os.makedirs(os.path.join(parameters.project_dir, "blank"))
    if not os.path.exists(os.path.join(parameters.project_dir, "sample")):
        os.makedirs(os.path.join(parameters.project_dir, "sample"))
    
    # Move files to the three folders by their names if there is any mzML or mzXML files in the project folder
    file_names = os.listdir(parameters.project_dir)
    file_names = [file_name for file_name in file_names if file_name.endswith(".mzML") or file_name.endswith(".mzXML")]
    if len(file_names) > 0:
        for i in file_names:
            if "qc" in i.lower() or 'pool' in i.lower():
                os.rename(os.path.join(parameters.project_dir, i), os.path.join(parameters.project_dir, "qc", i))
            elif "blank" in i.lower():
                os.rename(os.path.join(parameters.project_dir, i), os.path.join(parameters.project_dir, "blank", i))
            else:
                os.rename(os.path.join(parameters.project_dir, i), os.path.join(parameters.project_dir, "sample", i))
    
    # list the absolute paths of the files in the three folders
    qc_file_names = os.listdir(os.path.join(parameters.project_dir, "qc"))
    qc_file_names = [os.path.join(parameters.project_dir, "qc", file_name) for file_name in qc_file_names]
    blank_file_names = os.listdir(os.path.join(parameters.project_dir, "blank"))
    blank_file_names = [os.path.join(parameters.project_dir, "blank", file_name) for file_name in blank_file_names]
    sample_file_names = os.listdir(os.path.join(parameters.project_dir, "sample"))
    sample_file_names = [os.path.join(parameters.project_dir, "sample", file_name) for file_name in sample_file_names]

    file_names = qc_file_names + sample_file_names + blank_file_names

    # process files
    feature_list = process_files(file_names, parameters)
    
    # output feature list to a pickle file
    with open(os.path.join(parameters.project_dir, "project.pickle"), "wb") as f:
        pickle.dump(feature_list, f)


def load_project(project_dir):
    """
    Load a project from a project directory.

    Parameters
    ----------
    project_dir : str
        The project directory.
    """

    # load the project
    with open(os.path.join(project_dir, "project.pickle"), "rb") as f:
        feature_list = pickle.load(f)
    
    return feature_list


def targeted_workflow(parameters):
    """
    A workflow for targeted metabolomics research. The function tries to
    find the targeted features in the samples by m/z and retention time.

    Parameters
    ----------
    parameters : Params object
        The parameters for the workflow.    
    """

    # check if rhe targeted file list is there
    if parameters.targeted_list is None or not os.path.exists(parameters.targeted_list):
        raise ValueError("The targeted file list does not exist.")
    
    df = pd.read_csv(parameters.targeted_list)
    mz_seq = np.array(df["mz"])
    rt_seq = np.array(df["rt"])

    # get file names
    file_names = os.listdir(parameters.project_dir)
    file_names = [file_name for file_name in file_names if file_name.endswith(".mzML") or file_name.endswith(".mzXML")]
    file_names = [os.path.join(parameters.project_dir, file_name) for file_name in file_names]

    for file_name in file_names:
        d = feature_detection(file_name, parameters)
        rois = d.find_roi_by_mzrt(mz_seq, rt_seq)

        # get the ROI with the highest intensity
        roi = rois[0]
        for i in rois:
            if i.peak_height > roi.peak_height:
                roi = i

        # define the output
        output = {}

        
        