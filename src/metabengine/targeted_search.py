# Author: Hauxu Yu

# A module for targeted search

# Import modules
import raw_data_utils as raw
from .params import Params
import numpy as np


def target_search(file_name, target_mz_seq, target_rt_seq=None, mz_tol=0.01, rt_tol=0.5):

    d = raw.MSData()
    params = Params()
    d.read_raw_data(file_name, params)
    params.estimate_params(d)
    params.mz_tol_ms1 = mz_tol
    d.drop_ion_by_int(params)
    d.find_rois(params)
    d.cut_rois(params)

    # move all short rois to the regular roi list and clear the short roi list
    for roi in d.rois_short:
        d.rois.append(roi)
    d.rois_short = []
    
    if target_rt_seq is None:
        target_rt_seq = [False] * len(target_mz_seq)
    
    results = []

    for i in range(len(target_mz_seq)):
        matched_idx = []

        for j, roi in enumerate(d.rois):
            if abs(roi.mz - target_mz_seq[i]) < mz_tol:
                matched_idx.append(j)
        
        if len(matched_idx) == 0:
            results.append(None)
            continue

        ints = np.array([d.rois[idx].peak_height for idx in matched_idx])
        results.append(d.rois[matched_idx[np.argmax(ints)]])

        # # if retention time if not specified, search for the ROI with the highest intensity
        # if target_rt_seq[i] is False:

        #     ints = np.array([d.rois[idx].peak_height for idx in matched_idx])
        #     results.append(d.rois[matched_idx[np.argmax(ints)]])
        
        # # if retention time is specified, search for the ROI with the closest retention time
        # else:
        #     rts = np.array([d.rois[idx].rt for idx in matched_idx])
        #     rt_diff = abs(rts - target_rt_seq[i])
        #     results.append(d.rois[matched_idx[np.argmin(rt_diff)]])

    return results