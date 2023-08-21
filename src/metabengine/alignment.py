# Author: Hauxu Yu

# A module to align metabolic features from different samples

# Import modules
import numpy as np
from tqdm import tqdm


def alignement(feature_list, d):
    """
    A function to correct the retention time (RT) of the MS data to be aligned.

    Parameters
    ----------------------------------------------------------
    feature_list: list
        A list of features to be aligned.
    d: MSData
        The MS data to be aligned.
    """

    d.rois_mz_seq = np.array([roi.mz for roi in d.rois])
    d.rois_rt_seq = np.array([roi.rt for roi in d.rois])

    # Initiate the aligned features using the first MS data if feature_list is empty
    if len(feature_list) == 0:
        for roi in d.rois:
            aligned_feature = AlignedFeature()
            aligned_feature.initiate_feat(roi)
            feature_list.append(aligned_feature)

        # sort features in feature list by peak height from high to low
        feature_list.sort(key=lambda x: x.height_seq[0], reverse=True)

        return None
    
    file_count = len(feature_list[0].mz_seq)

    # If feature_list is not empty, align the MS data to the existing features by targeted search
    for feat in tqdm(feature_list):
        roi = find_roi_from_data(feat, d)

        if roi is not None:
            feat.extend_feat(roi)
        else:
            feat.extend_gap()
    
    # For newly detected features, add them to the feature list
    for roi in d.rois:
        aligned_feature = AlignedFeature()
        aligned_feature.initiate_feat(roi)
        aligned_feature.extend_gap_front(file_count)
        feature_list.append(aligned_feature)
    
    feature_list.sort(key=lambda x: x.height_seq[0], reverse=True)


class AlignedFeature:
    """
    A class to model an aligned feature from different files.
    """

    def __init__(self):
        """
        Function to initiate Aligned.
        """

        self.mz = 0.0
        self.rt = 0.0
        self.mz_seq = np.array([])
        self.rt_seq = np.array([])
        self.height_seq = np.array([])
        self.area_seq = np.array([])
        self.average_height_seq = np.array([])
        self.ms2_seq = []

        self.best_ms2 = None
    

    def initiate_feat(self, roi):

        self.mz = roi.mz
        self.rt = roi.rt
        self.mz_seq = np.array([roi.mz])
        self.rt_seq = np.array([roi.rt])
        self.height_seq = np.array([roi.peak_height])
        self.area_seq = np.array([roi.peak_area])
        self.average_height_seq = np.array([roi.peak_height_by_ave])
        self.ms2_seq = [roi.best_ms2]


    def extend_feat(self, roi):
        """
        A function to extend the feature with a new ROI.

        Parameters
        ----------------------------------------------------------
        roi: ROI object
            The new ROI to be added.
        """
        
        self.mz_seq = np.append(self.mz_seq, roi.mz)
        self.rt_seq = np.append(self.rt_seq, roi.rt)
        self.height_seq = np.append(self.height_seq, roi.peak_height)
        self.area_seq = np.append(self.area_seq, roi.peak_area)
        self.average_height_seq = np.append(self.average_height_seq, roi.peak_height_by_ave)
        self.ms2_seq.append(roi.best_ms2)

        self.mz = np.nanmean(self.mz_seq)
        self.rt = np.nanmean(self.rt_seq)
    

    def extend_gap(self):
        """
        A function to extend the feature with a gap.
        """

        self.mz_seq = np.append(self.mz_seq, np.nan)
        self.rt_seq = np.append(self.rt_seq, np.nan)
        self.height_seq = np.append(self.height_seq, 0.0)
        self.area_seq = np.append(self.area_seq, 0.0)
        self.average_height_seq = np.append(self.average_height_seq, 0.0)
        self.ms2_seq.append(None)
    

    def extend_gap_front(self, file_count):
        """
        A function to add gaps with the number of file_count to the front of the feature.

        Parameters
        ----------------------------------------------------------
        file_count: int
            The number of files that have been processed.
        """

        self.mz_seq = np.concatenate((np.full(file_count, np.nan), self.mz_seq))
        self.rt_seq = np.concatenate((np.full(file_count, np.nan), self.rt_seq))
        self.height_seq = np.concatenate((np.full(file_count, 0.0), self.height_seq))
        self.area_seq = np.concatenate((np.full(file_count, 0.0), self.area_seq))
        self.average_height_seq = np.concatenate((np.full(file_count, 0.0), self.average_height_seq))
        self.ms2_seq = [None] * file_count + self.ms2_seq

    
    def choose_best_ms2(self):
        """
        A function to choose the best MS2 for the feature. 
        The best MS2 is the one with the highest summed intensity.
        """

        total_ints = []

        for ms2 in self.ms2_seq:
            if ms2 is not None:
                total_ints.append(np.sum(ms2.prod_int_seq))
            else:
                total_ints.append(0.0)

        self.best_ms2 = self.ms2_seq[np.argmax(total_ints)]


def find_roi_from_data(feat, d):
    """
    A function to find the feature from the MS data.

    Parameters
    ----------------------------------------------------------
    feat: AlignedFeature object
        The feature to be found.
    d: MSData
        The MS data to be searched.
    """

    idx = None

    idxes = np.where(np.logical_and(np.abs(d.rois_mz_seq-feat.mz) < d.params.align_mz_tol_ms1, np.abs(d.rois_rt_seq-feat.rt) < d.params.align_rt_tol))[0]
    if len(idxes) == 0:
        return None
    
    int_seq = np.array([d.rois[i].peak_height for i in idxes])
    
    idx = idxes[np.argmax(int_seq)]

    d.rois_mz_seq = np.delete(d.rois_mz_seq, idx)
    d.rois_rt_seq = np.delete(d.rois_rt_seq, idx)
    
    # remove the ROI from the MS data
    return d.rois.pop(idx)









# Back up
# def alignement(feature_list, d, cross_file_params):
#     """
#     A function to correct the retention time (RT) of the MS data to be aligned.

#     Parameters
#     ----------------------------------------------------------
#     feature_list: list
#         A list of features to be aligned.
#     d: MSData
#         The MS data to be aligned.
#     Cross_file_params: Cross_file_params object
#         The parameters for data alignment.
#     """

#     # Initiate the aligned features using the first MS data if feature_list is empty
#     if len(feature_list) == 0:
#         for roi in d.rois:
#             aligned_feature = AlignedFeature()
#             aligned_feature.initiate_feat(roi)
#             feature_list.append(aligned_feature)

#         # sort features in feature list by peak height from high to low
#         feature_list.sort(key=lambda x: x.height_seq[0], reverse=True)

#         return None
    
#     file_count = len(feature_list[0].mz_seq)

#     # If feature_list is not empty, align the MS data to the existing features by targeted search
#     for feat in tqdm(feature_list):
#         roi = find_roi_from_data(feat, d, cross_file_params)

#         if roi is not None:
#             feat.extend_feat(roi)
#         else:
#             feat.extend_gap()
    
#     # For newly detected features, add them to the feature list
#     for roi in d.rois:
#         aligned_feature = AlignedFeature()
#         aligned_feature.initiate_feat(roi)
#         aligned_feature.extend_gap_front(file_count)
#         feature_list.append(aligned_feature)
    
#     feature_list.sort(key=lambda x: x.height_seq[0], reverse=True)