# Author: Hauxu Yu

# A module to align metabolic features from different samples

# Import modules
import numpy as np
import pandas as pd
from .visualization import mirror_ms2_db

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

    # Initiate the aligned features using the first MS data if feature_list is empty
    if len(feature_list) == 0:
        for roi in d.rois:
            aligned_feature = AlignedFeature()
            aligned_feature.extend_feat(roi)
            feature_list.append(aligned_feature)

        # sort features in feature list by peak height from high to low
        feature_list.sort(key=lambda x: x.highest_roi_intensity, reverse=True)

    else:
        existed_file_num = len(feature_list[0].mz_seq)

        mz_seq = np.array([roi.mz for roi in d.rois])
        rt_seq = np.array([roi.rt for roi in d.rois])
        int_seq = np.array([roi.peak_height for roi in d.rois])
        labeled_roi = np.ones(len(d.rois), dtype=bool)

        for feat in feature_list:
            v = np.logical_and(np.abs(mz_seq - feat.mz) < d.params.align_mz_tol, np.abs(rt_seq - feat.rt) <= d.params.align_rt_tol)
            v = np.where(np.logical_and(v, labeled_roi))[0]

            if len(v) == 0:
                feat.extend_feat(roi=None)
                continue

            if len(v) > 1:
                # select the one with the highest intensity
                v = v[np.argmax(int_seq[v])]
            else:
                v = v[0]

            feat.extend_feat(roi=d.rois[v])
            labeled_roi[v] = False
        
        # For newly detected features, add them to the feature list
        for i in range(len(labeled_roi)):
            if labeled_roi[i]:
                aligned_feature = AlignedFeature()
                aligned_feature.extend_feat(roi=d.rois[i], front_zeros=[0.0] * existed_file_num)
                feature_list.append(aligned_feature)
        
        feature_list.sort(key=lambda x: x.highest_roi_intensity, reverse=True)


class AlignedFeature:
    """
    A class to model an aligned feature from different files.
    """

    def __init__(self):
        """
        Define the attributes of a aligned feature.
        """

        self.mz = 0.0
        self.rt = 0.0
        self.mz_seq = []
        self.rt_seq = []
        self.peak_height_seq = []
        self.peak_area_seq = []
        self.top_average_seq = []
        self.ms2_seq = []
        self.best_ms2 = None
        self.highest_roi = None
        self.highest_roi_intensity = 0.0

        # annotation by MS2 matching
        self.annotation = None
        self.similarity = None
        self.matched_peak_number = None
        self.smiles = None
        self.inchikey = None
        self.matched_precursor_mz = None
        self.matched_peaks = None

        # isotope, in-source fragment, and adduct information
        self.charge_state = 1
        self.is_isotope = False
        self.is_in_source_fragment = False
        self.is_adduct = False

    def extend_feat(self, roi, front_zeros=[]):
        """
        A function to extend the feature with a new ROI.

        Parameters
        ----------------------------------------------------------
        roi: ROI object
            The new ROI to be added.
        zeros: list
            A list of zeros to be added to the feature.
        """
        
        if len(self.mz_seq) == 0:
            set_init_mzrt = True
        else:
            set_init_mzrt = False

        if len(front_zeros) > 0:
            self.mz_seq.extend(front_zeros)
            self.rt_seq.extend(front_zeros)
            self.peak_height_seq.extend(front_zeros)
            self.peak_area_seq.extend(front_zeros)
            self.top_average_seq.extend(front_zeros)
            self.ms2_seq.extend([None] * len(front_zeros))

        if roi is not None:
            self.mz_seq.append(roi.mz)
            self.rt_seq.append(roi.rt)
            self.peak_height_seq.append(roi.peak_height)
            self.peak_area_seq.append(roi.peak_area)
            self.top_average_seq.append(roi.top_average)
            self.ms2_seq.append(roi.best_ms2)

            if roi.peak_height > self.highest_roi_intensity:
                self.highest_roi_intensity = roi.peak_height
                self.highest_roi = roi
            
            if set_init_mzrt:
                self.mz = roi.mz
                self.rt = roi.rt

        else:
            self.mz_seq.append(0.0)
            self.rt_seq.append(0.0)
            self.peak_height_seq.append(0.0)
            self.peak_area_seq.append(0.0)
            self.top_average_seq.append(0.0)
            self.ms2_seq.append(None)

    
    def choose_best_ms2(self):
        """
        A function to choose the best MS2 for the feature. 
        The best MS2 is the one with the highest summed intensity.
        """

        total_ints = []
        for ms2 in self.ms2_seq:
            if ms2 is not None:
                total_ints.append(np.sum(ms2.peaks[:,1]))
            else:
                total_ints.append(0.0)
        self.best_ms2 = self.ms2_seq[np.argmax(total_ints)]
    

    def show_feature_info(self):
        """
        A function to show the information of the feature.
        """

        print('m/z: ', self.mz)
        print('RT: ', self.rt)
        print('Area sequence: ', self.peak_area_seq)
        print('Height sequence: ', self.peak_height_seq)
    

    def plot_match_result(self, output=False):

        if self.matched_peaks is not None:
            mirror_ms2_db(self, output=output)
        else:
            print("No matched MS/MS spectrum found.")


def sum_aligned_features(feature_list):
    """
    A function to summarize the aligned features.

    Parameters
    ----------------------------------------------------------
    feature_list: list
        A list of aligned features.   
    """

    for f in feature_list:
        f.mz_seq = np.array(f.mz_seq)
        f.rt_seq = np.array(f.rt_seq)
        f.mz = np.mean(f.mz_seq[f.mz_seq > 0])
        f.rt = np.mean(f.rt_seq[f.rt_seq > 0])

        f.choose_best_ms2()

        f.charge_state = f.highest_roi.charge_state
        f.is_isotope = f.highest_roi.is_isotope
        f.is_in_source_fragment = f.highest_roi.in_source_fragment
        f.adduct_type = f.highest_roi.adduct_type


def output_aligned_features(feature_list, file_names, path, int_values="peak_area"):
    """
    A function to output the aligned features.

    Parameters
    ----------------------------------------------------------
    feature_list: list
        A list of aligned features.
    output_path: str
        The path to the output file.
    """

    result = []

    for idx, f in enumerate(feature_list):
        roi = f.highest_roi
        
        iso_dist = ""
        for i in range(len(roi.isotope_mz_seq)):
            iso_dist += str(np.round(roi.isotope_mz_seq[i], decimals=4)) + ";" + str(np.round(roi.isotope_int_seq[i], decimals=0)) + "|"
        iso_dist = iso_dist[:-1]

        ms2 = ""
        if f.best_ms2 is not None:
            for i in range(len(f.best_ms2.peaks)):
                ms2 += str(np.round(f.best_ms2.peaks[i, 0], decimals=4)) + ";" + str(np.round(f.best_ms2.peaks[i, 1], decimals=0)) + "|"
            ms2 = ms2[:-1]

        if int_values.lower()=="peak_area":
            int_seq = f.peak_area_seq
        elif int_values.lower()=="peak_height":
            int_seq = f.peak_height_seq
        elif int_values.lower()=="top_average":
            int_seq = f.top_average_seq

        temp = [idx+1, f.mz, f.rt, ms2, f.charge_state, f.is_isotope, iso_dist,
                f.is_in_source_fragment, f.adduct_type, f.annotation, f.similarity, 
                f.matched_peak_number, f.smiles, f.inchikey, roi.quality]
                
        temp.extend(int_seq)

        result.append(temp)

    # convert result to a pandas dataframe
    columns = ["id", "mz", "rt", "ms2", "charge_state", "is_isotope", "isotope_dist",
                "in_source_fragment", "adduct_type", "annotation", "similarity_score", 
                "matched_peak_number", "smiles", "inchikey", "quality"]
    columns.extend(file_names)
    df = pd.DataFrame(result, columns=columns)
    
    # save the dataframe to csv file
    path = path + "aligned_feature_table.csv"
    df.to_csv(path, index=False)