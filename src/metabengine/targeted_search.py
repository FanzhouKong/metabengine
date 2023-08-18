# Author: Hauxu Yu

# A module for searching features in a targeted manner
# Searching is based on m/z, retention time, MS/MS, and existance of defined adducts

# Import modules
from . import raw_data_utils as raw
from .params import Params
import numpy as np


def target_search(d, targets, mz_tol=0.01, rt_tol=0.3, method=None):
    """
    A function to search for a feature in a targeted manner.

    Parameters
    ----------
    d: MSData object
        The MSData object to search. The MSData object should contain the ROI after peak picking.
    targets: list of dict
        A dictionary containing the target feature information. For example:
        {
            'name': "Caffeine",
            'formula': "C8H10N4O2",
            'retention_time': 0.0,
            'common_adducts': ["[M+H]+", "[M+NH4]+"],
            "common_adducts_mz": [195.0886, 217.0706],
            'preferred_adduct': "[M+H]+",
            'retention_time': [{
                'method': '5_min',
                'rt': 0.0
            },],
            'blank_int_seq': np.array([100.0, 200.0]),
        }
    """

    # d.rois should be sorted by m/z

    for target in targets:
        mzs = np.array(target['common_adducts_mz'])
        rt = None

        # check if the retention time of the specified method has been specified.
        if method is not None and 'retention_time' in target.keys():
            for item in target['retention_time']:
                if item['method'] == method:
                    rt = item['rt']
                    break
        
        for mz in mzs:
            matched_idx = np.isclose(d.rois_mz_seq, mz, atol=mz_tol).nonzero()[0]
            if len(matched_idx) == 0:
                continue

            if rt is not None:
                matched_idx = np.array([idx for idx in matched_idx if abs(d.rois[idx].rt - rt) < rt_tol])
            
            if len(matched_idx) == 0:
                continue

            # if there are multiple matched rois, choose 



    






# def target_search(file_name, target_mz_seq, target_rt_seq=None, mz_tol=0.01, rt_tol=0.5):

#     d = raw.MSData()
#     params = Params()
#     d.read_raw_data(file_name, params)
#     params.estimate_params(d)
#     params.mz_tol_ms1 = mz_tol
#     d.drop_ion_by_int(params)
#     d.find_rois(params)
#     d.cut_rois(params)

#     # move all short rois to the regular roi list and clear the short roi list
#     for roi in d.rois_short:
#         d.rois.append(roi)
#     d.rois_short = []
    
#     if target_rt_seq is None:
#         target_rt_seq = [False] * len(target_mz_seq)
    
#     results = []

#     for i in range(len(target_mz_seq)):
#         matched_idx = []

#         for j, roi in enumerate(d.rois):
#             if abs(roi.mz - target_mz_seq[i]) < mz_tol:
#                 matched_idx.append(j)
        
#         if len(matched_idx) == 0:
#             results.append(None)
#             continue

#         ints = np.array([d.rois[idx].peak_height for idx in matched_idx])
#         results.append(d.rois[matched_idx[np.argmax(ints)]])

#         # # if retention time if not specified, search for the ROI with the highest intensity
#         # if target_rt_seq[i] is False:

#         #     ints = np.array([d.rois[idx].peak_height for idx in matched_idx])
#         #     results.append(d.rois[matched_idx[np.argmax(ints)]])
        
#         # # if retention time is specified, search for the ROI with the closest retention time
#         # else:
#         #     rts = np.array([d.rois[idx].rt for idx in matched_idx])
#         #     rt_diff = abs(rts - target_rt_seq[i])
#         #     results.append(d.rois[matched_idx[np.argmin(rt_diff)]])

#     return results