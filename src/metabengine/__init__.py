# Author: Hauxu Yu

# A module to summarize the main data processing modules

# Import modules
from . import raw_data_utils as raw
from .params import Params
from .ann_feat_quality import predict_quality
from .feature_grouping import annotate_isotope, annotate_adduct, annotate_in_source_fragment
from .alignment import alignement, summarize_aligned_features, output_aligned_features
import os
from keras.models import load_model
from .annotation import annotate_features, annotate_rois
import pickle
import pandas as pd


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

    # create a MSData object
    d = raw.MSData()

    # read raw data
    d.read_raw_data(file_name, params)

    # drop ions by intensity (defined in params.int_tol)
    d.drop_ion_by_int()

    # detect region of interests (ROIs)
    d.find_rois()

    # cut ROIs by MS2 spectra
    if d.params.cut_roi:
        d.cut_rois()

    # sort ROI by m/z, find roi quality by length, find the best MS2
    d.process_rois()

    # predict feature quality. If the model is not loaded, load the model
    if d.params.ann_model is None:
        data_path_ann = os.path.join(os.path.dirname(__file__), 'model', "peak_quality_NN.keras")
        d.params.ann_model = load_model(data_path_ann)
    predict_quality(d)

    print("Number of extracted ROIs: " + str(len(d.rois)))

    # annotate isotopes, adducts, and in-source fragments
    annotate_isotope(d)
    annotate_in_source_fragment(d)
    annotate_adduct(d)

    # annotate MS2 spectra
    if annotation and d.params.msms_library is not None:
        annotate_rois(d)

    # output single file to a csv file
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

    # generate a list to store all the features
    feature_list = []

    # load the ANN model for peak quality prediction
    data_path_ann = os.path.join(os.path.dirname(__file__), 'model', "peak_quality_NN.keras")
    params.ann_model = load_model(data_path_ann)

    # process each file
    for file_name in file_names:
        print("Processing file: " + os.path.basename(file_name))
        # feature detection
        d = feature_detection(file_name, params)

        # remove the features with scan number < 5 and no MS2 from the feature alignment
        d.rois = [roi for roi in d.rois if roi.length >= 5 or roi.best_ms2 is not None]
        print("Number of ROIs for alignment: " + str(len(d.rois)))

        # feature alignment
        alignement(feature_list, d)
        print("-----------------------------------")
    
    # summarize aligned features
    summarize_aligned_features(feature_list)

    # annotation
    if params.msms_library is not None:
        annotate_features(feature_list, params)

    # output aligned features to a csv file
    if params.output_aligned_file:
        output_file_names = [os.path.basename(file_name) for file_name in file_names]
        output_file_names = [os.path.splitext(file_name)[0] for file_name in output_file_names]
        output_aligned_features(feature_list, file_names, params.project_dir, normalization=params.run_normalization)

    return feature_list


def read_raw_file_to_obj(file_name, params=None):
    """
    Read a raw file to a MSData object.
    It's a useful function for data visualization or brief data analysis.

    Parameters
    ----------
    file_name : str
        The file name of the raw file.
    """

    # create a MSData object
    d = raw.MSData()

    # read raw data
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

    file_names, sample_groups = _untargeted_workflow_preparation(parameters)
    parameters.sample_groups = sample_groups

    # process files
    feature_list = process_files(file_names, parameters)
    
    # output feature list to a pickle file
    with open(parameters.project_dir + "mbe_project.mbe", "wb") as f:
        pickle.dump(feature_list, f)

    return feature_list


def load_project(project_file):
    """
    Load a project from a project directory.

    Parameters
    ----------
    project_file : str
        The path to the project file with .mbe format.
    """

    # load the project
    with open(project_file, "rb") as f:
        feature_list = pickle.load(f)
    
    return feature_list


def _untargeted_workflow_preparation(parameters):
    
    # Check the folder for creating the project
    if not os.path.exists(parameters.project_dir):
        raise ValueError("The project directory does not exist.")
    
    if parameters.project_dir[-1] != "/":
        parameters.project_dir += "/"
    
    sample_dir = parameters.project_dir + "sample/"
    single_file_dir = parameters.project_dir + "single_file_output/"
    
    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)
    if not os.path.exists(single_file_dir):
        os.makedirs(single_file_dir)
    
    # Move files to the sample folder if not moved
    file_names = os.listdir(parameters.project_dir)
    file_names = [file_name for file_name in file_names if file_name.endswith(".mzML") or file_name.endswith(".mzXML")]
    if len(file_names) > 0:
        for i in file_names:
            os.rename(os.path.join(parameters.project_dir, i), os.path.join(parameters.project_dir, "sample", i))
    
    # Check the raw files are loaded
    file_names = os.listdir(sample_dir)
    if len(file_names) == 0:
        raise ValueError("No raw files are found in the project directory.")
    # Get raw file extension
    file_ext = os.path.splitext(file_names[0])[1]
    
    # Sort the file names to process the data in the order of QC, sample, and blank
    # Check if the sample table is available
    if not os.path.exists(os.path.join(parameters.project_dir, "sample_table.csv")):
        print("No sample table is found in the project directory. All samples will be treated as regular samples.")
    else:
        sample_table = pd.read_csv(os.path.join(parameters.project_dir, "sample_table.csv"))
    
    # Sort the sample table by the order of sample, QC, and blank
    sample_groups = []
    for i in range(len(sample_table)):
        sample_groups.append((sample_table.iloc[i, 0], sample_table.iloc[i, 1]))
    qc_names = [i[0] for i in sample_groups if i[1].lower() == "qc"]
    blank_names = [i[0] for i in sample_groups if i[1].lower() == "blank"]
    sample_names = [i[0] for i in sample_groups if i[1].lower() != "qc" and i[1].lower() != "blank"]
    file_names = sample_names + qc_names + blank_names
    file_names = [sample_dir + name + file_ext for name in file_names]

    return file_names, sample_groups


def _bin_detection(file_name, params):
    """
    Bin detection from a raw LC-MS file (.mzML or .mzXML).
    A bin has to have MS2 spectra.

    Parameters
    ----------
    file_name : str
        File name of the raw file.
    parameters : Params object
        The parameters for the workflow.
    """

    # create a MSData object
    d = raw.MSData()

    # read raw data
    d.read_raw_data(file_name, params)

    # drop ions by intensity (defined in params.int_tol)
    d.drop_ion_by_int()

    # detect region of interests (ROIs)
    d.find_rois()

    # cut ROIs by MS2 spectra
    if d.params.cut_roi:
        d.cut_rois()

    # sort ROI by m/z, find roi quality by length, find the best MS2
    d._process_rois_bin_generation()
    # predict feature quality. If the model is not loaded, load the model
    if d.params.ann_model is None:
        data_path_ann = os.path.join(os.path.dirname(__file__), 'model', "peak_quality_NN.keras")
        d.params.ann_model = load_model(data_path_ann)
    predict_quality(d)

    print("Number of extracted ROIs: " + str(len(d.rois)))

    # annotate isotopes, adducts, and in-source fragments
    annotate_isotope(d)
    annotate_in_source_fragment(d)
    annotate_adduct(d)

    d._discard_isotopes()

    # output single file to a csv file
    d.output_single_file()

    return d


def process_files_for_bin_generation(file_names, params):
    """
    A function to process multiple raw files.

    Parameters
    ----------
    file_names : list
        A list of file names of the raw files in .mzML or .mzXML format.
    params : Params object
        The parameters for the workflow.
    """

    # load the ANN model for peak quality prediction
    data_path_ann = os.path.join(os.path.dirname(__file__), 'model', "peak_quality_NN.keras")
    params.ann_model = load_model(data_path_ann)

    # process each file
    for file_name in file_names:
        print("Processing file: " + os.path.basename(file_name))
        # feature detection
        d = _bin_detection(file_name, params)

        output_file_name = os.path.basename(file_name)
        output_file_name = os.path.splitext(output_file_name)[0]
        with open(params.project_dir + "processed_data/" + output_file_name + ".pkl", "wb") as f:
            pickle.dump(d, f)