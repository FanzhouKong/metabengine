from pyteomics import mass
from ms_entropy import calculate_entropy_similarity
import numpy as np


def cal_ion_mass(formula, adduct, charge):
    """
    A function to calculate the ion mass of a given formula, adduct and charge.

    Parameters
    ----------------------------------------------------------
    formula: str
        The chemical formula of the ion.
    adduct: str
        Adduct of the ion.
    charge: int
        Charge of the ion. 
        Use signs for specifying ion modes: +1 for positive mode and -1 for negative mode.

    Returns
    ----------------------------------------------------------
    ion_mass: float
        The ion mass of the given formula, adduct and charge.
    """

    # if there is a D in the formula, and not followed by a lowercase letter, replace it with H[2]
    if 'D' in formula and not formula[formula.find('D') + 1].islower():
        formula = formula.replace('D', 'H[2]')

    # calculate the ion mass
    final_formula = formula + adduct
    ion_mass = (mass.calculate_mass(formula=final_formula) - charge * _ELECTRON_MASS) / abs(charge)
    return ion_mass

_ELECTRON_MASS = 0.00054858


def ms2_grouping(ms2_list, precursor_mz_tol=0.01, similarity_tol=0.8):
    """
    A function to group ms2 by precursor m/z and entropy similarity. (TODO)

    Parameters
    ----------------------------------------------------------
    ms2_list: list
        A list of MS2 spectra. Each MS2 spectrum is a dictionary with
        keys: 'precursor_mz', 'peaks'. For example:
        {'precursor_mz': 100.1012, 'peaks': [[69.7894, 100], [100.1012, 200]]}
    precursor_mz_tol: float
        The tolerance for precursor m/z grouping. Default is 0.01.
    similarity_tol: float
        The threshold for entropy similarity. Default is 0.8.
    """

    return None

    # sort ms2 by precursor m/z
    ms2_list.sort(key=lambda x: x['precursor_mz'])
    mz_seq = np.array([ms2['precursor_mz'] for ms2 in ms2_list])

    # group ms2 by precursor m/z
    ms2_groups = []

    for idx, ms2 in enumerate(ms2_list):
        matched = np.where(abs(mz_seq - ms2['precursor_mz']) < precursor_mz_tol)[0]
        matched = matched[matched > idx]
        if len(matched) == 0:
            ms2_groups.append([ms2])
            continue
        else:
            cur_group = [ms2]
            for i in matched:
                if calculate_entropy_similarity(ms2['peaks'], ms2_list[i]['peaks']) > similarity_tol:
                    cur_group.append(ms2_list[i])
            ms2_groups.append(cur_group)
    

