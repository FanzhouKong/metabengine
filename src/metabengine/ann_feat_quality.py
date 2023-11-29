# Author: Hauxu Yu

# A module to load a trained model to predict the quality of features
# Prediction is based on peak shape

# Import modules
import numpy as np
from scipy.interpolate import interp1d


def predict_quality(d, model=None, threshold=0.5):
    """
    Function to predict the quality of a feature as an ROI.

    Parameters
    ----------------------------------------------------------
    d: MSData object
        An MSData object that contains the MS data.
    model: keras model
        A keras model that is trained to predict the quality of a feature.
    threshold: float
        A threshold to determine the quality of a feature.
    """

    if model is None:
        model = d.params.ann_model

    temp = np.array([peak_interpolation(roi.int_seq) for roi in d.rois])
    q = model.predict(temp, verbose=0)[:,0] > threshold

    for i in range(len(d.rois)):
        # if the roi quality is not good, then skip and don't overwrite
        if d.rois[i].quality == 'good' and q[i] == 0:
            d.rois[i].quality = 'bad peak shape'


def peak_interpolation(peak):
    '''
    A function to interpolate a peak to a vector of a given size.

    Parameters
    ----------------------------------------------------------
    peak: numpy array
        A numpy array that contains the peak to be interpolated.
    '''
    
    peak_interp_rule = interp1d(np.arange(len(peak)), peak, kind='linear')
    interp_seed = np.linspace(0, len(peak)-1, 64)
    peak_interp = peak_interp_rule(interp_seed)

    peak_interp = peak_interp / np.max(peak_interp)

    return peak_interp