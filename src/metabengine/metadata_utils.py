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