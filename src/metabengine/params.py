# Author: Hauxu Yu

# A module to define and estimate the parameters

# Import
import numpy as np

# Define a class to store the parameters
class Params:
    """
    A class to store the parameters for individual files.
    """

    def __init__(self):
        """
        Function to initiate Params.
        ----------------------------------------------------------
        """

        # Need to be specified by the user
        self.project_dir = None   # Project directory, character string
        self.targeted = False     # Whether the workflow is targeted, default is False
        self.targeted_list = None

        self.rt_range = [0.0, 1000.0]   # RT range in minutes, list of two numbers
        self.mode = "dda"         # Acquisition mode, "dda", "dia", or "full_scan"
        self.ms2_sim_tol = 0.8    # MS2 similarity tolerance
        self.ion_mode = "positive"   # Ionization mode, "positive" or "negative"

        # Parameters for feature detection
        self.mz_tol_ms1 = 0.01    # m/z tolerance for MS1, default is 0.01
        self.mz_tol_ms2 = 0.015   # m/z tolerance for MS2, default is 0.015
        self.int_tol = 1000       # Intensity tolerance, default is 10000 for Orbitrap and 1000 for other instruments
        self.roi_gap = 2          # Gap within a feature, default is 2 (i.e. 2 consecutive scans without signal)
        self.min_ion_num = 10     # Minimum scan number a feature, default is 10
        self.ann_model = None     # Model for feature annotation, keras.src.engine.sequential.Sequential object
        self.cut_roi = True       # Whether to cut ROI, default is True
        self.discard_short_roi = True   # Whether to discard short ROI with length < 5 and without MS2, default is Tr

        # Parameters for feature alignment
        self.align_mz_tol = 0.01        # m/z tolerance for MS1, default is 0.01
        self.align_rt_tol = 0.2         # RT tolerance, default is 0.2

        # Parameters for feature annotation
        self.msms_library = None   # MS/MS library in MSP format, character string
        self.ppr = 0.8             # Peak peak correlation threshold, default is 0.8

        # Parameters for output
        self.output_single_file = False   # Whether to output a single file for each raw file, default is False
        self.output_aligned_file = True   # Output aligned file path, character string


    def show_params_info(self):
        """
        Function to print the parameters.
        ----------------------------------------------------------
        """

        print("m/z tolerance (MS1):", self.mz_tol_ms1)
        print("m/z tolerance (MS2):", self.mz_tol_ms2)
        print("Intensity tolerance:", self.int_tol)
        print("ROI gap:", self.roi_gap)
        print("MS2 similarity tolerance:", self.ms2_sim_tol)
        print("Acquisition mode:", self.mode)
        print("Retention time range:", self.rt_range)
        print("Project directory:", self.proj_dir)