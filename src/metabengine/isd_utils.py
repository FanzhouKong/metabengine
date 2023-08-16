# Author: Hauxu Yu

# A module to provide utilities for internal standards
from pyteomics.mass import calculate_mass
from chemparse import parse_formula
import re

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


    def calculate_mz_for_adducts(self):
        """
        Calculate the m/z values for internal standards based on the provided formula and common adducts.
        """

        mz = 0.0

        # if the formula and preferred adduct are provided, calculate the m/z value
        if self.formula and len(self.common_adducts) > 0:
            parsed = parse_formula(self.formula)

            if "D" in parsed.keys():
                parsed["H[2]"] = parsed.pop("D")
            if "T" in parsed.keys():
                parsed["H[3]"] = parsed.pop("T")
            
            mz = calculate_mass(parsed)
        else:
            return None
        
        # calculate the m/z value for the preferred adduct
        for adduct in self.common_adducts:
            
            # get charge state
            charge_state = adduct.split("]")[1]
            if len(charge_state) == 1:
                charge_state += "1"
            else:
                charge_state = charge_state[-1] + charge_state[:-1]
            charge_state = int(charge_state)

            # get formation of adduct: number of molecule, added items, removed items
            adduct_formation = re.findall(r'\[(.*?)\]', adduct)[0]
            molecule_number = adduct_formation.split("M")[0]
            if len(molecule_number) == 0:
                molecule_number = 1
            else:
                molecule_number = int(molecule_number)

            adduct_formation = adduct_formation.split("M")[1]

            if len(adduct_formation) == 0:
                mz_final = mz
            else:
                mz_final = mz + _mass_offset[adduct_formation]

            mz_final = (mz_final - charge_state*ELECTRON_MASS) / abs(charge_state)

            self.common_adducts_mz.append(mz_final)

