# Author: Hauxu Yu

# A module to handle the metadata of the dataset

import os

class SampleTable:
    """
    A class to define the properties of samples.    
    """

    def __init__(self):
        """
        Define the properties of samples.        
        """

        self.type = None    # "sample", "qc", "blank", "others"
        self.names = None    # sample name, with extension of mzXML or mzML


    def load_sample_name(self, dir):
        """
        Load the sample name in the directory with extension of mzXML or mzML.

        Parameters
        ----------
        dir : str
            The directory of the raw MS data.   
        """

        self.names = os.listdir(dir)
        self.names = [n + dir for n in self.names if n.endswith(".mzXML") or n.endswith(".mzML")]
    

    def determine_sample_type(self):
        """
        Determine the sample type based on sample names.       
        """

        if self.names is None:
            return None
        
        for n in self.names:
            if n.lower().startswith("QC"):
                self.type = "qc"
            elif n.lower().startswith("Blank"):
                self.type = "blank"
            elif n.lower().startswith("Sample"):
                self.type = "sample"
            else:
                self.type = "others"


class AnalyticalMethod:
    """
    A class to define the analytical method for, e.g., LC-MS/MS.
    """

    def __init__(self):
        # LC method
        self.column_name = None    # column type
        self.column_length = None    # column length
        self.column_inner_diameter = None    # column inner diameter
        self.column_particle_size = None    # column particle size
        self.column_temperature = None    # column temperature
        self.column_flow_rate = None    # column flow rate
        self.column_gradient = None    # column gradient
        self.column_injection_volume = None    # column injection volume

        # MS method
        self.ms_name = None    # MS name
        self.ms_ion_mode = None    # MS ionization mode, "pos", "neg", or "mixed"
        self.ms_scan_type = None    # MS scan type, "full_scan", "dda", or "dia"
        self.ms1_scan_range = None    # MS1 scan range
        self.ms2_scan_range = None    # MS2 scan range
    

    def load_method_from_raw(self, dir):
        """
        Load the analytical method from the raw MS data.

        Parameters
        ----------
        dir : str
            The directory of the raw MS data.   
        """

        pass