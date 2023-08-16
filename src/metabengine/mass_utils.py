# Author: Hauxu Yu

# A module to support 
from pyteomics.mass import calculate_mass
from chemparse import parse_formula
import re

def calculate_mz_for_adducts(formula, adduct):
    """
    Calculate the m/z values for chemicals based on formula and adduct.

    Parameters
    ----------
    formula : str
        The formula of the chemical.
    adduct : str
        The adduct of the chemical.
    """

    mz = 0.0

    # if the formula and preferred adduct are provided, calculate the m/z value
    try:
        parsed = parse_formula(formula)
    except:
        return None

    if "D" in parsed.keys():
        parsed["H[2]"] = parsed.pop("D")
    if "T" in parsed.keys():
        parsed["H[3]"] = parsed.pop("T")
    
    mz = calculate_mass(parsed)
      
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
    
    mz = mz * molecule_number

    adduct_formation = adduct_formation.split("M")[1]

    if len(adduct_formation) == 0:
        mz = mz
    else:
        mz = mz + _mass_offset[adduct_formation]

    return (mz - charge_state*ELECTRON_MASS) / abs(charge_state)



_mass_offset = {
    "+H": 1.007825,
    "+Na": 22.989769,
    "+K": 38.963707,
    "+NH4": 18.034374,
    "+H-H2O": -17.002740,
    "+Cl": 34.968853,
    "-H": -1.007825,
    '-H-H2O': -19.018390,
    "+HCOO": 44.997654,
    '+CH3COO': 59.013304,
}

ELECTRON_MASS = 0.000549