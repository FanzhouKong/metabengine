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

        Method for estimating the parameters:
        1. m/z tolerance:
            1) Find the m/z that has the highest intensity in all MS1 scans
            2) Calculate its standard deviation
            3) Use the standard deviation * 5 as the m/z tolerance
        2. intensity tolerance:
            1) Find the last 5% of the lowest intensity in all MS1 scans
            2) Calculate the mean and standard deviation of the intensity
            3) Use the mean + 3 * standard deviation as the intensity tolerance

        Parameters
        ----------------------------------------------------------
        d: MSData object
            An MSData object that contains the MS data.
        """

        wanted_mz = 0.0
        wanted_rt = 0.0
        highest_int = 0.0

        int_ms1 = np.array([], dtype=int)

        if estimate_mz_tol:
            for i in d.ms1_idx:
                scan = d.scans[i]
                mz_highest_int = scan.mz_seq[np.argmax(scan.int_seq)]
                int_highest_int = np.max(scan.int_seq)
                if int_highest_int > highest_int:
                    wanted_mz = mz_highest_int
                    wanted_rt = scan.rt
                    highest_int = int_highest_int
            mz_seq = d.get_eic_data(mz=wanted_mz, rt=wanted_rt, mz_tol=0.005, rt_tol=1.0)[2]
            mz_seq = mz_seq[mz_seq > 0]
            if np.std(mz_seq) * 5 < 0.002:
                self.mz_tol_ms1 = 0.002
                self.mz_tol_ms2 = 0.004
            else:
                self.mz_tol_ms1 = np.std(mz_seq) * 5
                self.mz_tol_ms2 = 2 * self.mz_tol_ms1
            
        if estimate_int_tol:
            for i in d.ms1_idx:
                scan = d.scans[i]
                # get the last 5% of the lowest intensity
                int_ms1 = np.append(int_ms1, np.sort(scan.int_seq)[:(int(len(scan.int_seq)*0.05)+1)])
            # Estimate the intensity tolerance
            self.int_tol = np.mean(int_ms1) + 3 * np.std(int_ms1)


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
