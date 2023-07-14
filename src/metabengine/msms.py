# Author: Hauxu Yu

# A module for MS/MS data processing

# Import modules
import numpy as np


def spec_similarity(spec1, spec2, mz_tol=0.01, sim_method='cosine', search_type='exact', filter_low_int=True, exclude_precursor=True, remove_precursor_mz=1.6):
    '''
    A function to calculate the dot product of two mass spectra.

    Parameters
    ----------------------------------------------------------
    spec1: Dict
        A dictionary that contains the first MS/MS spectrum.
    spec2: Dict
        A dictionary that contains the second MS/MS spectrum.
    mz_tol: float
        The m/z tolerance for MS/MS spectra. Default is 0.01.
    sim_method: str
        The method to calculate the spectral similarity. Default is cosine.
        'cosine': cosine similarity or dot product.
        'entropy': entropy similarity.
    search_type: str
        The type of search. Default is exact.
        'exact': exact search.
        'open': open search.
        'neutral': neutral loss search.
        'hybrid': hybrid search.
    filter_low_int: bool
        Whether to filter out low intensity peaks. Default is True.
    exclude_precursor: bool
        Whether to exclude precursor peak in MS/MS for similarity calculation. Default is True.
    '''

    if len(spec1.prod_mz_seq) == 0 or len(spec2.prod_mz_seq) == 0:
        return "Empty spectrum detected."
    
    mz1 = spec1.prod_mz_seq
    mz2 = spec2.prod_mz_seq
    int1= spec1.prod_int_seq
    int2= spec2.prod_int_seq

    # remove precursor peak
    if exclude_precursor:
        if search_type == 'hybrid':
            smaller_mz = np.min([spec1.precs_mz, spec2.precs_mz])
            int1= int1[mz1 < smaller_mz-remove_precursor_mz]
            int2= int2[mz2 < smaller_mz-remove_precursor_mz]
            mz1 = mz1[mz1 < smaller_mz-remove_precursor_mz]
            mz2 = mz2[mz2 < smaller_mz-remove_precursor_mz]
        else:
            int1= int1[mz1 < spec1.precs_mz-remove_precursor_mz]
            int2= int2[mz2 < spec2.precs_mz-remove_precursor_mz]
            mz1 = mz1[mz1 < spec1.precs_mz-remove_precursor_mz]
            mz2 = mz2[mz2 < spec2.precs_mz-remove_precursor_mz]
    
    if len(mz1) == 0 or len(mz2) == 0:
        return "Empty spectrum after filtering."

    # filter out low intensity peaks
    if filter_low_int:
        mz1 = mz1[int1 > np.max(int1) * 0.01]
        mz2 = mz2[int2 > np.max(int2) * 0.01]
        int1= int1[int1 > np.max(int1) * 0.01]
        int2= int2[int2 > np.max(int2) * 0.01]

    # scale int1 and int2 respectively to make inner product = 1
    int1 = int1 / np.inner(int1, int1) ** 0.5
    int2 = int2 / np.inner(int2, int2) ** 0.5

    score = 0.0

    id1 = id2 = id1_ = id2_ = []

    # align two spectra based on search type
    if search_type == 'exact' or search_type == 'open':
        [id1, id2] = align_frag(mz1, mz2, mz_tol=mz_tol, mz_diff=0.0)
    elif search_type == 'neutral':
        [id1, id2] = align_frag(mz1, mz2, mz_tol=mz_tol, mz_diff=spec2.precs_mz-spec1.precs_mz)
    elif search_type == 'hybrid':
        [id1, id2] = align_frag(mz1, mz2, mz_tol=mz_tol, mz_diff=0.0)
        mz1_rest = mz1[np.isin(np.arange(len(mz1)), id1, invert=True)]
        mz2_rest = mz2[np.isin(np.arange(len(mz2)), id2, invert=True)]
        if len(mz1_rest) > 0 and len(mz2_rest) > 0:
            [id1_, id2_] = align_frag(mz1_rest, mz2_rest, mz_tol=mz_tol, mz_diff=spec2.precs_mz-spec1.precs_mz)
    
    # calculate spectral similarity
    if sim_method == 'cosine':
        if len(id1) > 0 and len(id2) > 0:
            score = np.inner(int1[id1], int2[id2])
        
        if len(id1_) > 0 and len(id2_) > 0:
            int1_rest= int1[np.isin(np.arange(len(mz1)), id1, invert=True)]
            int2_rest= int2[np.isin(np.arange(len(mz2)), id2, invert=True)]
            score += np.inner(int1_rest[id1_], int2_rest[id2_])
    
    elif sim_method == 'entropy':
        pass

    return score


def align_frag(mz1, mz2, mz_tol=0.01, mz_diff=0.0):
    '''
    A function to align two mass spectra.

    Parameters
    ----------------------------------------------------------
    mz1: np.array
        A numpy array that contains the m/z values of the first MS/MS spectrum.
    mz2: np.array
        A numpy array that contains the m/z values of the second MS/MS spectrum.
    mz_tol: float
        The m/z tolerance for MS/MS spectra. Default is 0.01.
    mz_diff: float
        The m/z difference for matching. Default is 0.0.

    Returns
    ----------------------------------------------------------
    [id1, id2]: list
        A list of indices of matched peaks.
    '''

    matched_pairs = []

    mz1 = mz1 + mz_diff

    mz_all = np.concatenate((mz1, mz2))
    id_all = np.concatenate((np.arange(len(mz1))+1, -np.arange(len(mz2))-1))

    # sort mz_all
    order = np.argsort(mz_all)
    mz_all = mz_all[order]
    id_all = id_all[order]
    
    mz_itvs = mz_all[1:] - mz_all[:-1]

    mz_itvs[id_all[1:] * id_all[:-1] > 0] = np.inf
    mz_itvs[mz_itvs > mz_tol] = np.inf

    while np.min(mz_itvs) < np.inf:
        p = np.argmin(mz_itvs)
        matched_pairs.append([id_all[p], id_all[p+1]])
        mz_itvs[p] = np.inf
        if p < len(mz_itvs) - 1:
            mz_itvs[p+1] = np.inf
        if p > 0:
            mz_itvs[p-1] = np.inf

    id1 = []
    id2 = []

    for pair in matched_pairs:
        if pair[0] > 0:
            id1.append(pair[0]-1)
            id2.append(-pair[1]-1)
        else:
            id1.append(pair[1]-1)
            id2.append(-pair[0]-1)
    
    return [id1, id2]