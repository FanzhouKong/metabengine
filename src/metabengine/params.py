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

        self.rt_range = [0.0, 60.0]   # RT range in minutes, list of two numbers
        self.mode = "dda"         # Acquisition mode, "dda", "dia", or "full_scan"
        self.ms2_sim_tol = 0.7    # MS2 similarity tolerance
        self.ion_mode = "pos"     # Ionization mode, "pos" or "neg"

        self.output_single_file_path = None   # Output single file path, character string

        # Parameters for feature detection
        self.mz_tol_ms1 = 0.01    # m/z tolerance for MS1, default is 0.01
        self.mz_tol_ms2 = 0.015   # m/z tolerance for MS2, default is 0.015
        self.int_tol = 1000       # Intensity tolerance, default is 10000 for Orbitrap and 1000 for other instruments
        self.roi_gap = 2          # Gap within a feature, default is 2 (i.e. 2 consecutive scans without signal)
        self.min_ion_num = 5      # Minimum scan number a feature, default is 5

        # Parameters for feature alignment
        self.align_mz_tol_ms1 = 0.01  # m/z tolerance for MS1, default is 0.01
        self.align_rt_tol = 0.1       # RT tolerance, default is 0.1

        # Parameters for feature annotation
        self.msms_library = None   # MS/MS library in MSP format, character string


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