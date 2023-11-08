# Author: Hauxu Yu

# A module to annotate metabolites based on their MS/MS spectra

# Import modules
import os


def laod_msms_db(path):
    """
    A function to load the MS/MS database in MSP format.

    Parameters
    ----------
    path : str
        The path to the MS/MS database in MSP format.    
    """

    # get extension of path
    ext = os.path.splitext(path)[1]

    if ext.lower() == '.msp':
        with open(path, 'r') as f:
            lines = f.readlines()












def msp_parser(path):
    """
    A function to parse the msp file into a list of dictionary.

    Parameters
    ----------
    path : str
        The path to the msp file.   
    """

    





