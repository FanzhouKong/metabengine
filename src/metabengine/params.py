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
        self.ms2_sim_tol = 0.75  # MS2 similarity tolerance
        self.mode = "dda"   # Acquisition mode, "dda", "dia", or "full_scan"
        self.ion_mode = "pos"   # Ionization mode, "pos", "neg" or "mixed"

        # Will be estimated by the program
        self.mz_tol_ms1 = 0.01
        self.mz_tol_ms2 = 0.02
        self.cycle_time = 0.0
        self.int_tol = 0

        # Constant
        self.roi_gap = 2
        self.min_ion_num = 5
    

    def estimate_params(self, d, estimate_mz_tol=True, estimate_cycle_time=True, estimate_int_tol=True):
        """
        Function to estimate the parameters from MS data.
        The parameters to be estimated include:
            m/z tolerance: self.mz_tol_ms1, self.mz_tol_ms2
            cycle time: self.cycle_time
            intensity tolerance: self.int_tol

        Parameters
        ----------------------------------------------------------
        d: MSData object
            An MSData object that contains the MS data.
        """

        mz_ms1_diff = 100.0

        int_ms1 = np.array([], dtype=int)

        for i in d.ms1_idx:
            scan = d.scans[i]
            temp = np.diff(np.sort(scan.mz_seq))

            if estimate_mz_tol:
                if len(temp) > 0:
                    if np.min(temp) < mz_ms1_diff:
                        mz_ms1_diff = np.min(temp)
            
            if estimate_int_tol:
                int_ms1 = np.append(int_ms1, np.min(scan.int_seq))
        
        if estimate_mz_tol:
            # Estimate the m/z tolerance
            self.mz_tol_ms1 = mz_ms1_diff * 0.99
            self.mz_tol_ms2 = mz_ms1_diff * 1.99

        if estimate_cycle_time:
            # Estimate the cycle time
            self.cycle_time = np.median(np.diff(d.ms1_rt_seq))

        if estimate_int_tol:
            # Estimate the intensity tolerance
            self.int_tol = int(np.mean(int_ms1) + 10*np.std(int_ms1))


    def print_params(self):
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
    

class Cross_file_params:
    """
    A class to store the parameters for multiple files.
    """

    def __init__(self):

        self.mz_tol_ms1 = 0.01
        self.mz_tol_ms2 = 0.02

        self.rt_tol = 0.3
