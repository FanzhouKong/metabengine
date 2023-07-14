# Author: Hauxu Yu

# A module to load a trained model to predict the quality of features
# Prediction is based on peak shape

# Import modules
import numpy as np
from scipy.interpolate import interp1d
from keras.models import model_from_json


def predict_quality(d, threshold=0.5):
    """
    Function to predict the quality of a feature as an ROI.

    Parameters
    ----------------------------------------------------------
    d: MSData object
        An MSData object that contains the MS data.
    """

    model = model_from_json(open('model/model_architecture.json').read())
    model.load_weights('model/model_weights.h5')

    temp = np.array([peak_interpolation(roi.int_seq) for roi in d.rois])
    q = model.predict(temp, verbose=0)[:,0] > threshold

    for i in range(len(d.rois)):
        if q[i] == 1:
            d.rois[i].quality = 'good'
        else:
            d.rois[i].quality = 'bad peak shape'


def peak_interpolation(peak):
    '''
    A function to interpolate a peak to a vector of a given size.

    Parameters
    ----------------------------------------------------------
    peak: numpy array
        A numpy array that contains the peak to be interpolated.
    '''
    
    temp = np.array([0,0])
    temp = np.insert(temp, 1, peak)
    peak_interp_rule = interp1d(np.arange(len(temp)), temp, kind='linear')
    interp_seed = np.linspace(0, len(temp)-1, 64)
    peak_interp = peak_interp_rule(interp_seed)

    return(peak_interp)
