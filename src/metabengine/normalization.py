# Author: Hauxu Yu

# A module for sample normalization

import numpy as np

def find_normalization_factors(array, method='pqn'):
    """ 
    A function to find the normalization factors for a data frame.

    Parameters
    ----------
    array : numpy array
        The data to be normalized.
    method : str
        The method to find the normalization factors.
        'pqn': probabilistic quotient normalization.

    Returns
    -------
    numpy array
        Normalization factor.
    """

    # find the reference sample
    ref_idx = find_reference_sample(array)

    if method == 'pqn':
        array = array[array[:, ref_idx] != 0, :]
        # calculate the normalization factor
        return np.median(array / array[:, ref_idx][:, None], axis=0)
    

def sample_normalization_by_factors(array, v):
    """
    A function to normalize a data frame by a vector.

    Parameters
    ----------
    array : numpy array
        The data to be normalized.
    v : numpy array
        The normalization factor.
    """

    # change all zeros to ones
    v[v == 0] = 1

    return array / v


def find_reference_sample(array, method='number'):
    """
    A function to find the reference sample for normalization.
    Note, samples are in columns and features are in rows.

    Parameters
    ----------
    array : numpy array
        The data to be normalized.
    method : str
        The method to find the reference sample. 
        'number': the reference sample has the most detected features.
        'total_intensity': the reference sample has the highest total intensity.
        'median_intensity': the reference sample has the highest median intensity.

    Returns
    -------
    int
        The index of the reference sample.
    """

    if method == 'number':
        # find the reference sample with the most detected features
        return np.argmax(np.count_nonzero(array, axis=0))
    elif method == 'total_intensity':
        # find the reference sample with the highest total intensity
        return np.argmax(np.sum(array, axis=0))
    elif method == 'median_intensity':
        # find the reference sample with the highest median intensity
        return np.argmax(np.median(array, axis=0))
    

def normalize_feature_list(feature_list, method='pqn', blank_sample_idx=None):
    """
    A function to normalize samples using a feature list.

    Parameters
    ----------
    feature_list : list
        A list of features.
    method : str
        The method to find the normalization factors.
        'pqn': probabilistic quotient normalization.

    Returns
    -------
    numpy array
        Normalization factors.
    """

    # Two steps: 1. find the normalization factors using feature with good peak shape
    #            2. normalize all features using the normalization factors

    array_good = np.array([f.peak_height_seq for f in feature_list if f.quality == 'good'])
    v = find_normalization_factors(array_good, method=method)
    v[v == 0] = 1

    # don't normalize blank samples
    if blank_sample_idx is not None:
        v[blank_sample_idx] = 1

    for f in feature_list:
        for i in range(len(f.peak_height_seq)):
            f.peak_height_seq[i] /= v[i]
            f.peak_area_seq[i] /= v[i]
            f.top_average_seq[i] /= v[i]
    return v