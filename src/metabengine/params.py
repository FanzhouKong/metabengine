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

        self.rt_range = [0.0, 1000.0]   # RT range in minutes, list of two numbers
        self.mode = "dda"         # Acquisition mode, "dda", "dia", or "full_scan"
        self.ms2_sim_tol = 0.8    # MS2 similarity tolerance
        self.ion_mode = "positive"   # Ionization mode, "positive" or "negative"

        # Parameters for feature detection
        self.mz_tol_ms1 = 0.01    # m/z tolerance for MS1, default is 0.01
        self.mz_tol_ms2 = 0.015   # m/z tolerance for MS2, default is 0.015
        self.int_tol = 1000       # Intensity tolerance, default is 30000 for Orbitrap and 1000 for other instruments
        self.roi_gap = 2          # Gap within a feature, default is 2 (i.e. 2 consecutive scans without signal)
        self.min_ion_num = 10     # Minimum scan number a feature, default is 10
        self.cut_roi = True       # Whether to cut ROI, default is True

        # Parameters for feature alignment
        self.align_mz_tol = 0.01        # m/z tolerance for MS1, default is 0.01
        self.align_rt_tol = 0.2         # RT tolerance, default is 0.2
        self.discard_short_roi = True   # Whether to discard short ROI with length < 5 and without MS2 from feature alignment, default is True
        self.quant_method = "peak_height"   # Quantification method, "peak_height", "peak_area", or "top_average"

        # Parameters for feature annotation
        self.msms_library = None   # MS/MS library in MSP format, character string
        self.ppr = 0.7             # Peak peak correlation threshold, default is 0.7

        # Parameters for normalization
        self.run_normalization = False   # Whether to normalize the data, default is False
        self.normalization_method = "pqn"   # Normalization method, "pqn". See module normalization.py for details

        # Parameters for output
        self.output_single_file = False   # Whether to output a single file for each raw file, default is False
        self.output_aligned_file = True   # Output aligned file path, character string

        # Statistical analysis
        self.run_statistical_analysis = True   # Whether to perform statistical analysis, default is True

        # Other parameters (only change if necessary)
        self.ann_model = None     # ANN model for peak quality prediction, default is None


    def __str__(self):
        """
        Print the parameters.
        ----------------------------------------------------------
        """

        print("Project directory: " + self.project_dir)
        print("RT range: " + str(self.rt_range))
        print("Acquisition mode: " + self.mode)
        print("MS2 similarity tolerance: " + str(self.ms2_sim_tol))
        print("Ionization mode: " + self.ion_mode)
        print("m/z tolerance for MS1: " + str(self.mz_tol_ms1))
        print("m/z tolerance for MS2: " + str(self.mz_tol_ms2))
        print("Intensity tolerance: " + str(self.int_tol))
        print("Gap within a feature: " + str(self.roi_gap))
        print("Minimum scan number a feature: " + str(self.min_ion_num))
        print("Whether to cut ROI: " + str(self.cut_roi))
        print("Whether to discard short ROI: " + str(self.discard_short_roi))
        print("ANN model for peak quality prediction: " + str(self.ann_model))
        print("m/z tolerance for MS1: " + str(self.align_mz_tol))
        print("RT tolerance: " + str(self.align_rt_tol))
        print("MS/MS library in MSP format: " + str(self.msms_library))
        print("Peak peak correlation threshold: " + str(self.ppr))
        print("Whether to output a single file for each raw file: " + str(self.output_single_file))
        print("Output aligned file path: " + str(self.output_aligned_file))