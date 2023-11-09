# Author: Hauxu Yu

# A module to group metabolic features from unique compounds
# 1. annotate isotopes
# 2. annotate adducts
# 3. annotate in-source fragments

# Import modules
import numpy as np


def annotate_isotope(d, params):
    """
    Function to annotate isotopes in the MS data.
    
    Parameters
    ----------------------------------------------------------
    d: MSData object
        An MSData object that contains the detected rois to be grouped.
    """

    # rank the rois (d.rois) in each file by m/z
    d.rois = d.rois.sort(key=lambda x: x.mz)
    mz_seq = np.array([roi.mz for roi in d.rois])
    rt_seq = np.array([roi.rt for roi in d.rois])

    # check rois
    for idx, r in enumerate(d.rois):
        mz = r.mz + 1.003355
        rt = r.rt

        mz_seq = d.scans[r.scan_number]
        
        while True:
            if abs(mz - m) < params.mz_tol_ms1:
                break


    # one compound generates different ions at the same time
    # here, we group the rois by their apexes' scan numbers with a tolerance of 2 scans (i.e., they should have apexes differing by 2 scans at most)

    # generate a vector to record whether a roi is grouped
    group_bool = np.zeros(len(d.rois), dtype=bool)
    
    for idx1, roi1 in enumerate(d.rois):
        # if the roi has been grouped, skip
        if group_bool[idx1]:
            continue
        
        # mark the roi as grouped
        group_bool[idx1] = True

        found_child = []

        # find the rois that belong to the same feature
        for idx2, roi2 in enumerate(d.rois[idx1+1:]):
            # if the roi has been grouped, skip
            if group_bool[idx1+1+idx2]:
                continue
            
            # if the apex the second roi is larger than the first roi by more than 2 scans, skip
            if roi2.scan_number - roi1.scan_number > 2:
                break

            known_mz_diff_output = known_mz_diff(roi1.mz, roi2.mz, params)
            
            # if the m/z differnce between the two rois is unknown, skip
            if not known_mz_diff_output:
                continue

            # if the peak-peak correlation between the two rois is less than 0.9, skip
            if peak_peak_correlation(roi1, roi2) < 0.9:
                continue

            # if pass all the above tests

            roi2.isotope["isotope"] = True
            roi2.isotope["isotope_number"] = known_mz_diff_output[1]
            roi2.isotope["isotope_element"] = known_mz_diff_output[2]
            roi2.isotope["parent_roi_id"] = roi1.id
            found_child.append(roi2.id)
        
        # if the roi has a son, mark the roi as the parent
        if len(found_child) > 0:
            roi1.isotope["isotope"] = True
            roi1.isotope["isotope_number"] = 0
            roi1.isotope["child_roi_id"] = found_child



def known_mz_diff(mz1, mz2, params):
    """
    A function to check whether the m/z difference between two rois is known.
    The second m/z should be larger than the first m/z.

    Parameters
    ----------------------------------------------------------
    mz1: float
        The m/z value of the first roi. (mz2 > mz1)
    mz2: float
        The m/z value of the second roi. (mz2 > mz1)
    params: Params object
        A Params object that contains the parameters.
    """

    mz_diff = mz2 - mz1

    # check isotopes
    fds = np.round(mz_diff / _isotopic_mass_diffence['mass_diffs'])
    diff = abs(_isotopic_mass_diffence['mass_diffs']*fds - mz_diff)

    idx_min_diff = np.argmin(diff)
    min_diff = diff[idx_min_diff]

    if min_diff < params.mz_tol_ms1 and fds[idx_min_diff] < 5:
        return ['isotope', fds[idx_min_diff], _isotopic_mass_diffence['elements'][idx_min_diff]]

    return False


def peak_peak_correlation(roi1, roi2):
    """
    A function to find the peak-peak correlation between two rois.

    Parameters
    ----------------------------------------------------------
    roi1: ROI object
        An ROI object.
    roi2: ROI object
        An ROI object.
    
    Returns
    ----------------------------------------------------------
    pp_cor: float
        The peak-peak correlation between the two rois.
    """

    # find the common scans in the two rois
    common_scans = np.intersect1d(roi1.scan_idx_seq, roi2.scan_idx_seq)

    if len(common_scans) < 2:
        return 1.0

    # find the intensities of the common scans in the two rois
    int1 = roi1.int_seq[np.isin(roi1.scan_idx_seq, common_scans)]
    int2 = roi2.int_seq[np.isin(roi2.scan_idx_seq, common_scans)]

    # calculate the correlation
    pp_cor = np.corrcoef(int1, int2)[0, 1]

    return pp_cor


_isotopic_mass_diffence = {
    'elements': ['H', 'C', 'N', 'O', 'S', 'Cl'],
    'mass_diffs': np.array([1.006277, 1.003355, 0.997035, 2.004246, 1.995796, 1.99705])
}

# _adduct_mass_diffence_neg = {
#     '-H': -1.007276,
#     '-H-H2O': -19.01839,
#     '+Na-2H': 20.974666,
#     '+Cl': 34.969402,
#     '+K-2H': 36.948606,
#     "+HCOO": 44.998201,
#     '+CH3COO': 59.013851,
#     '+Br': 78.918885,
#     '+CF3COO': 112.985586,
#     '2M-H': -1.007276,
#     '3M-H': -1.007276,
# }

# _adduct_mass_diffence_pos = {
#     '+H': 1.007276,
#     '+H-H2O': 19.01839,
#     '+Na': 22.989218,
#     '+K': 38.963158,
#     '+NH4': 18.033823,
#     '+CH3OH+H': 33.033489,
#     '+ACN+H': 42.033823,
#     '+2ACN+H': 83.060370,
#     '+ACN+Na': 64.015765,
#     '+Li': 7.016003,
#     "+Ag": 106.905093,
#     '+2Na-H': 44.97116,
#     '+2K-H': 76.919039,
#     '+IPA+H': 61.06534,
# }
