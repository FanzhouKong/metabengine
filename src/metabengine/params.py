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

        wanted_mz_seq = []
        wanted_mz = 0.0
        highest_int = 0.0

        int_ms1 = np.array([], dtype=int)

        for i in d.ms1_idx:
            scan = d.scans[i]
            if estimate_mz_tol:
                highest_mz_temp = scan.mz_seq[np.argmax(scan.int_seq)]
                highest_int_temp = np.max(scan.int_seq)
                if abs(highest_mz_temp - wanted_mz) > 0.1:
                    if highest_int_temp > highest_int:
                        wanted_mz_seq = [highest_mz_temp]
                        wanted_mz = highest_mz_temp
                        highest_int = highest_int_temp
                else:
                    wanted_mz_seq.append(highest_mz_temp)
                    if highest_int_temp > highest_int:
                        wanted_mz = highest_mz_temp
                        highest_int = highest_int_temp
            
            if estimate_int_tol:
                # get the last 5% of the lowest intensity
                int_ms1 = np.append(int_ms1, np.sort(scan.int_seq)[:(int(len(scan.int_seq)*0.05)+1)])
        
        if estimate_mz_tol:
            # Estimate the m/z tolerance
            wanted_mz_seq = np.array(wanted_mz_seq)
            self.mz_tol_ms1 = np.std(wanted_mz_seq) * 5
            self.mz_tol_ms2 = 2 * self.mz_tol_ms1

        if estimate_int_tol:
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
