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
        self.rt_range = np.array([0, np.inf], dtype=float)   # in minute
        self.proj_dir = None    # Project directory, character string
        self.ms2_sim_tol = 0.7  # MS2 similarity tolerance
        self.mode = "dda"   # Acquisition mode, "dda", "dia", or "full_scan"
        self.ion_mode = "pos"   # Ionization mode, "pos", "neg" or "mixed"

        # Will be estimated by the program
        self.mz_tol_ms1 = 0.008
        self.mz_tol_ms2 = 0.015
        self.cycle_time = 0.0
        self.int_tol = 0

        # Constant
        self.roi_gap = 2
        self.min_ion_num = 5

        # Alignment parameters
        self.align_mz_tol_ms1 = 0.01
        self.align_mz_tol_ms2 = 0.02
        self.align_rt_tol = 0.3


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