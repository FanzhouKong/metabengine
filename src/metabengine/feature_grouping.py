# Author: Hauxu Yu

# A module to group metabolic features from unique compounds
# 1. annotate isotopes
# 2. annotate adducts
# 3. annotate in-source fragments

# Import modules
import numpy as np


def annotate_isotope(d):
    """
    Function to annotate isotopes in the MS data.
    
    Parameters
    ----------------------------------------------------------
    d: MSData object
        An MSData object that contains the detected rois to be grouped.
    """

    # rank the rois (d.rois) in each file by m/z
    d.rois.sort(key=lambda x: x.mz)
    mz_seq = np.array([roi.mz for roi in d.rois])
    scan_seq = np.array([roi.scan_number for roi in d.rois])

    labeled_roi = np.zeros(len(d.rois), dtype=bool)

    for idx, r in enumerate(d.rois):
        
        if labeled_roi[idx]:
            idx += 1
            continue

        r = d.rois[idx]
        r.isotope_int_seq = [r.peak_height]
        r.isotope_mz_seq = [r.mz]

        # go to that scan and determine the charge state
        isotopes, _, charge_state = _find_iso_from_scan(d.scans[r.scan_number], r.mz)
        r.charge_state = charge_state

        # find roi using isotope list
        for j, iso in enumerate(isotopes):
            v = np.where(np.logical_and(np.abs(mz_seq - iso) < 0.005, np.abs(scan_seq - r.scan_number) <= 2))[0]

            if len(v) == 0:
                continue

            if len(v) > 1:
                # select the one with the lowest scan difference
                v = v[np.argmin(np.abs(scan_seq[v] - r.scan_number))]
            else:
                v = v[0]
            
            cor = peak_peak_correlation(r, d.rois[v])

            if cor > 0.9:
                labeled_roi[v] = True
                
                r.isotope_int_seq.append(d.rois[v].peak_height)
                r.isotope_mz_seq.append(d.rois[v].mz)

                d.rois[v].isotope_state = j+1


def annotate_in_source_fragment(d):
    """
    Function to annotate in-source fragments in the MS data.
    Only [M+O] (roi.isotope_state=0) will be considered in this function.
    Two criteria are used to annotate in-source fragments:
    1. The precursor m/z of the child is in the MS2 spectrum of the parent.
    2. Peak-peak correlation > 0.9
    
    Parameters
    ----------------------------------------------------------
    d: MSData object
        An MSData object that contains the detected rois to be grouped.
    params: Params object
        A Params object that contains the parameters.
    """

    # sort ROI by m/z from high to low
    d.rois.sort(key=lambda x: x.mz, reverse=True)

    # find in-source fragments
    for idx, r in enumerate(d.rois):

        if r.isotope_state != 0:
            continue

        mz_seq = r.best_ms2.peaks[:, 0]

        for m 







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


def _find_iso_from_scan(scan, mz):
    """
    Find the charge state of a m/z value based on a scan.  
    """

    isotopes = []
    distribution = []
    mass_diff = scan.mz_seq - mz
    charge_state = 1
    for idx, md in enumerate(mass_diff):
        if md < 0.04:
            continue
        if md > 10:
            break

        tmp = md/(1.003355/2)
        a = round(tmp)
        if abs(tmp-a) < 0.012:
            isotopes.append(scan.mz_seq[idx])
            distribution.append(scan.int_seq[idx])
            if a%2 == 1:
                charge_state = 2

    return isotopes, distribution, charge_state



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
