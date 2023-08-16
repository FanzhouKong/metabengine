# Author: Hauxu Yu

# A module to provide utilities for internal standards
from pyteomics.mass import calculate_mass
from chemparse import parse_formula
import json


class InternalStandards:
    """
    A class defines a internal standard.
    """

    def __init__(self):
        """
        Define the initial attributes of internal standard.
        
        Parameters
        ----------------------------------------------------------        
        """

        self.name = None
        self.formula = None

        self.inchi = None
        self.inchikey = None
        self.smiles = None

        self.common_adducts = []
        self.common_adducts_mz = []
        self.preferred_adduct = None
        
        self.see_positive_mode = False
        self.see_negative_mode = False

        self.retention_time = None

        self.commercial_source = None


    def define_internal_std(self, **kwargs):
        """
        Define an internal standard by providing the required information.

        Parameters
        ----------------------------------------------------------
        kwargs: dict
            A dictionary contains the required information of internal standard.        
        """

        for key, value in kwargs.items():
            # if key is an attribute of internal standard, assign the value to the attribute
            if hasattr(self, key):
                setattr(self, key, value)