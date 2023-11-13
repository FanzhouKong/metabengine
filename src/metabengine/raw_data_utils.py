# Author: Hauxu Yu

# A module to read and process the raw MS data
# Classes are defined in order to handle the data

# Import modules
from pyteomics import mzml, mzxml
import numpy as np
import os
from . import peak_detect
import matplotlib.pyplot as plt
import pandas as pd


class MSData:
    """
    A class that models a single file (mzML or mzXML) and
    processes the raw data.

    Attributes
    ----------------------------------------------------------
    scans: list, a list of Scan objects
    ms1_rt_seq: numpy array, retention times of all MS1 scans
    bpc_int: numpy array, intensity of the BPC
    rois: list, a list of ROI objects
    rois_mz_seq: numpy array, m/z of all ROIs
    params: Params object, a Params object that contains the parameters
    """


    def __init__(self):
        """
        Function to initiate MSData.
        ----------------------------------------------------------
        """

        self.scans = []   # A list of MS scans
        self.ms1_rt_seq = []  # Retention times of all MS1 scans
        self.bpc_int = [] # Intensity of the BPC
        self.rois = []  # A list of ROIs
        self.params = None  # A Params object
        self.file_name = None  # File name of the raw data without extension


    def read_raw_data(self, file_name, params):
        """
        Function to read raw data to MS1 and MS2 (if available)
        (supported by pyteomics package).

        Parameters
        ----------------------------------------------------------
        file_name: str
            File name of raw MS data (mzML or mzXML).
        params: Params object
            A Params object that contains the parameters.
        """

        self.params = params

        if os.path.isfile(file_name):
            # get extension from file name
            ext = os.path.splitext(file_name)[1]

            self.file_name = os.path.splitext(os.path.basename(file_name))[0]

            if ext.lower() == ".mzml":
                with mzml.MzML(file_name) as reader:
                    self.extract_scan_mzml(reader)
            elif ext.lower() == ".mzxml":
                with mzxml.MzXML(file_name) as reader:
                    self.extract_scan_mzxml(reader)
            else:
                print("Unsupported raw data format. Raw data must be in mzML or mzXML.")
        else:
            print("File does not exist.")


    def extract_scan_mzml(self, spectra):
        """
        Function to extract all scans and convert them to Scan objects.

        Parameters
        ----------------------------------------------------------
        spectra: pyteomics object
            An iteratable object that contains all MS1 and MS2 scans.
        """

        idx = 0     # Scan number
        self.ms1_idx = []   # MS1 scan index
        self.ms2_idx = []   # MS2 scan index

        rt_unit = spectra[0]['scanList']['scan'][0]['scan start time'].unit_info

        # Iterate over all scans
        for spec in spectra:
            # Get the retention time and convert to minute
            try:
                rt = spec['scanList']['scan'][0]['scan start time']
            except:
                rt = spec['scanList']['scan'][0]['scan time']
            
            if rt_unit == 'second':
                rt = rt / 60

            # Check if the retention time is within the range
            if self.params.rt_range[0] < rt < self.params.rt_range[1]:
                if spec['ms level'] == 1:
                    temp_scan = Scan(level=1, scan=idx, rt=rt)
                    mz_array = spec['m/z array']
                    int_array = spec['intensity array']

                    temp_scan.add_info_by_level(mz_seq=mz_array, int_seq=int_array)
                    self.ms1_idx.append(idx)

                    # update base peak chromatogram
                    self.bpc_int.append(np.max(spec['intensity array']))
                    self.ms1_rt_seq.append(rt)

                elif spec['ms level'] == 2:
                    temp_scan = Scan(level=2, scan=idx, rt=rt)
                    precursor_mz = spec['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                    peaks = np.array([spec['m/z array'], spec['intensity array']], dtype=np.float64).T
                    temp_scan.add_info_by_level(precursor_mz=precursor_mz, peaks=peaks)
                    _clean_ms2(temp_scan, int_threshold=self.params.int_tol)
                    self.ms2_idx.append(idx)
                
                self.scans.append(temp_scan)
                idx += 1

        # print the number of extracted ms1 and ms2 scans
        print(str(len(self.ms1_idx)) + " MS1 and " + str(len(self.ms2_idx)) + " MS2 scans extracted.")


    def extract_scan_mzxml(self, spectra):
        """
        Function to extract all scans and convert them to Scan objects.

        Parameters
        ----------------------------------------------------------
        spectra: pyteomics object
            An iteratable object that contains all MS1 and MS2 scans.
        """

        idx = 0     # Scan number
        self.ms1_idx = []   # MS1 scan index
        self.ms2_idx = []   # MS2 scan index

        rt_unit = spectra[0]['scanList']['scan'][0]['scan start time'].unit_info

        # Iterate over all scans
        for spec in spectra:
            # Get the retention time and convert to minute
            rt = spec["retentionTime"]    # retention time of mzXML is in minute

            if rt_unit == 'second':
                rt = rt / 60

            # Check if the retention time is within the range
            if self.params.rt_range[0] < rt < self.params.rt_range[1]:
                if spec['msLevel'] == 1:
                    temp_scan = Scan(level=1, scan=idx, rt=rt)
                    mz_array = spec['m/z array']
                    int_array = spec['intensity array']

                    temp_scan.add_info_by_level(mz_seq=mz_array, int_seq=int_array)
                    self.ms1_idx.append(idx)

                    # update base peak chromatogram
                    self.bpc_int = np.append(self.bpc_int, np.max(spec['intensity array']))
                    self.ms1_rt_seq = np.append(self.ms1_rt_seq, rt)

                elif spec['msLevel'] == 2:
                    temp_scan = Scan(level=2, scan=idx, rt=rt)
                    precursor_mz = spec['precursorMz'][0]['precursorMz']
                    peaks = np.array([spec['m/z array'], spec['intensity array']], dtype=np.float64).T
                    temp_scan.add_info_by_level(precursor_mz=precursor_mz, peaks=peaks)
                    _clean_ms2(temp_scan, int_threshold=self.params.int_tol)
                    self.ms2_idx.append(idx)
                
                self.scans.append(temp_scan)
                idx += 1

        # print the number of extracted ms1 and ms2 scans
        print(str(len(self.ms1_idx)) + " MS1 and " + str(len(self.ms2_idx)) + " MS2 scans extracted.")

    
    def drop_ion_by_int(self):
        """
        Function to drop ions by intensity.

        Parameters
        ----------------------------------------------------------
        tol: int
            Intensity tolerance.
        """

        for idx in self.ms1_idx:
            self.scans[idx].mz_seq = self.scans[idx].mz_seq[self.scans[idx].int_seq > self.params.int_tol]
            self.scans[idx].int_seq = self.scans[idx].int_seq[self.scans[idx].int_seq > self.params.int_tol]


    def find_rois(self):
        """
        Function to find ROI in MS1 scans.

        Parameters
        ----------------------------------------------------------
        params: Params object
            A Params object that contains the parameters.
        """

        self.rois = peak_detect.roi_finder(self, self.params)
    

    def cut_rois(self, return_cut_rois=False):
        """
        Function to cut ROI into smaller pieces.
        """

        small_rois = []
        to_be_removed = []

        pairs = []
        for i, roi in enumerate(self.rois):
            
            positions = peak_detect.find_roi_cut(roi, self.params)
            
            if positions is not None:

                # append each item in a list to small_rois
                small_rois_tmp = peak_detect.roi_cutter(roi, positions)
                small_rois.extend(small_rois_tmp)
                to_be_removed.append(i)

                pairs.append([roi, small_rois_tmp])
        
        # remove the original rois
        for i in sorted(to_be_removed, reverse=True):
            del self.rois[i]

        # append the small rois to the original rois
        self.rois.extend(small_rois)

        if return_cut_rois:
            return pairs


    def process_rois(self):
        """
        Function to process ROIs.

        Parameters
        ----------------------------------------------------------
        params: Params object
            A Params object that contains the parameters.
        """

        # 1. sort ROI by m/z
        self.rois.sort(key=lambda x: x.mz)

        for idx, roi in enumerate(self.rois):
            roi.id = idx
            roi.sum_roi()

            roi.int_seq = np.array(roi.int_seq)
            
            # 1. find roi quality by length
            if roi.length >= self.params.min_ion_num:
                roi.quality = 'good'
            else:
                roi.quality = 'short'
            
            # 2. find the best MS2
            roi.find_best_ms2()
        
        if self.params.discard_short_roi:
            self.rois = [roi for roi in self.rois if roi.quality != 'short']
            print("Number of regular ROIs after discarding short ROIs: " + str(len(self.rois)))
    

    def plot_bpc(self, output=False):
        """
        Function to plot base peak chromatogram.

        Parameters
        ----------------------------------------------------------
        output: str
            Output file name. If not specified, the plot will be shown.
        """

        plt.figure(figsize=(10, 3))
        plt.rcParams['font.size'] = 14
        plt.rcParams['font.family'] = 'Arial'
        plt.plot(self.ms1_rt_seq, self.bpc_int, linewidth=1, color="black")
        plt.xlabel("Retention Time (min)", fontsize=18, fontname='Arial')
        plt.ylabel("Intensity", fontsize=18, fontname='Arial')
        plt.xticks(fontsize=14, fontname='Arial')
        plt.yticks(fontsize=14, fontname='Arial')

        if output:
            plt.savefig(output, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()
        
    
    def output_roi_report(self, path, detailed=True):
        """
        Function to generate a report for rois in csv format.'

        Parameters
        ----------------------------------------------------------
        path: str
            Output file path.
        detailed: bool, default False
            True: output detailed information for each ROI.
            False: output only basic information for each ROI.
        """

        result = []

        for roi in self.rois:
            iso_dist = ""
            for i in len(roi.isotope_mz_seq):
                iso_dist += str(np.round(roi.isotope_mz_seq[i], decimals=4)) + ";" + str(np.round(roi.isotope_int_seq[i], decimals=0)) + "|"
            iso_dist = iso_dist[:-1]

            ms2 = ""
            for i in len(roi.best_ms2.peaks):
                ms2 += str(np.round(roi.best_ms2.peaks[i, 0], decimals=4)) + ";" + str(np.round(roi.best_ms2.peaks[i, 1], decimals=0)) + "|"
            ms2 = ms2[:-1]

            temp = [roi.id+1, roi.mz, roi.rt, roi.length, roi.rt_seq[0],
                    roi.rt_seq[1], roi.peak_area, roi.peak_height,
                    roi.peak_height_by_ave, ms2,
                    roi.charge_state, roi.isotope_state, iso_dist,
                    roi.in_source_fragment, roi.isf_parent_roi_id, roi.isf_child_roi_id,
                    roi.adduct_type, roi.adduct_parent_roi_id, roi.adduct_child_roi_id,
                    roi.quality]
            result.append(temp)

        # convert result to a pandas dataframe
        df = pd.DataFrame(result, columns=["ID", "m/z", "RT", "Length", "RT_start",
                                           "RT_end", "Peak_area", "Peak_height",
                                           "Peak_height_by_ave", 
                                           "MS2", "Charge_state", "Isotope_state",
                                           "Isotope_distribution", "In_source_fragment",
                                           "ISF_parent_ID", "ISF_child_ID", "Adduct_type",
                                           "Adduct_parent_ID", "Adduct_child_ID", "Quality"])
        
        # save the dataframe to csv file
        path = path + self.file_name + ".csv"
        df.to_csv(path, index=False)
    

    def get_eic_data(self, target_mz, mz_tol=0.005, rt_range=[0, np.inf]):
        
        eic_rt = np.array([])
        eic_int = np.array([])
        eic_mz = np.array([])
        eic_scan_idx = np.array([])

        for i in self.ms1_idx:
            if self.scans[i].rt > rt_range[0]:
                mz_diff = np.abs(self.scans[i].mz_seq - target_mz)
                if np.min(mz_diff) < mz_tol:
                    eic_rt = np.append(eic_rt, self.scans[i].rt)
                    eic_int = np.append(eic_int, self.scans[i].int_seq[np.argmin(mz_diff)])
                    eic_mz = np.append(eic_mz, self.scans[i].mz_seq[np.argmin(mz_diff)])
                    eic_scan_idx = np.append(eic_scan_idx, i)
                else:
                    eic_rt = np.append(eic_rt, self.scans[i].rt)
                    eic_int = np.append(eic_int, 0.0)
                    eic_mz = np.append(eic_mz, 0.0)
                    eic_scan_idx = np.append(eic_scan_idx, i)

            if self.scans[i].rt > rt_range[1]:
                break
        
        return eic_rt, eic_int, eic_mz, eic_scan_idx


    def find_roi_by_mzrt(self, mz, rt, mz_tol=0.005, rt_tol=0.3):
        rois = []
        for roi in self.rois:
            if np.abs(roi.mz - mz) < mz_tol and np.abs(roi.rt - rt) < rt_tol:
                roi.show_roi_info()
                print("a total of " + str(len(roi.rt_seq)) + " scans")
                print("roi start: " + str(roi.rt_seq[0]) and "roi end: " + str(roi.rt_seq[-1]))
                print("------------------")
                rois.append(roi)

        return rois     
    

    def plot_eic(self, target_mz, mz_tol=0.005, rt_range=[0, np.inf], output=False):
        """
        Function to plot EIC of a target m/z.
        """

        # get the eic data
        eic_rt, eic_int, _, eic_scan_idx = self.get_eic_data(target_mz, mz_tol=mz_tol, rt_range=rt_range)

        plt.figure(figsize=(10, 3))
        plt.rcParams['font.size'] = 14
        plt.rcParams['font.family'] = 'Arial'
        plt.plot(eic_rt, eic_int, linewidth=1, color="black")
        plt.xlabel("Retention Time (min)", fontsize=18, fontname='Arial')
        plt.ylabel("Intensity", fontsize=18, fontname='Arial')
        plt.xticks(fontsize=14, fontname='Arial')
        plt.yticks(fontsize=14, fontname='Arial')

        if output:
            plt.savefig(output, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()
        
        return eic_rt[np.argmax(eic_int)], np.max(eic_int), eic_scan_idx[np.argmax(eic_int)]
        
    
    def find_ms2_by_mzrt(self, mz_target, rt_target, mz_tol=0.01, rt_tol=0.3, return_best=False):
        """
        Function to find MS2 scan by precursor m/z and retention time.

        Parameters
        ----------------------------------------------------------
        mz_target: float
            Precursor m/z.
        rt_target: float
            Retention time.
        mz_tol: float
            m/z tolerance.
        rt_tol: float
            Retention time tolerance.
        return_best: bool
            True: only return the best MS2 scan with the highest intensity.
            False: return all MS2 scans as a list.
        """

        matched_ms2 = []

        for id in self.ms2_idx:
            rt = self.scans[id].rt

            if rt < rt_target - rt_tol:
                continue
            
            mz = self.scans[id].precursor_mz
            
            if abs(mz - mz_target) < mz_tol and abs(rt - rt_target) < rt_tol:
                matched_ms2.append(self.scans[id])
        
            if rt > rt_target + rt_tol:
                break

        if return_best:
            if len(matched_ms2) > 1:
                total_ints = [np.sum(ms2.peaks[:,1]) for ms2 in matched_ms2]
                return matched_ms2[np.argmax(total_ints)]
            elif len(matched_ms2) == 1:
                return matched_ms2[0]
            else:
                return None
        else:
            return matched_ms2      


    def plot_roi(self, roi_idx, mz_tol=0.005, rt_range=[0, np.inf], rt_window=None, output=False):
        """
        Function to plot EIC of a target m/z.
        """

        if rt_window is not None:
            rt_range = [self.rois[roi_idx].rt - rt_window, self.rois[roi_idx].rt + rt_window]

        # get the eic data
        eic_rt, eic_int, _, eic_scan_idx = self.get_eic_data(self.rois[roi_idx].mz, mz_tol=mz_tol, rt_range=rt_range)
        idx_start = np.where(eic_scan_idx == self.rois[roi_idx].scan_idx_seq[0])[0][0]
        idx_end = np.where(eic_scan_idx == self.rois[roi_idx].scan_idx_seq[-1])[0][0] + 1

        plt.figure(figsize=(7, 3))
        plt.rcParams['font.size'] = 14
        plt.rcParams['font.family'] = 'Arial'
        plt.plot(eic_rt, eic_int, linewidth=1, color="black")
        plt.fill_between(eic_rt[idx_start:idx_end], eic_int[idx_start:idx_end], color="black", alpha=0.4)
        plt.xlabel("Retention Time (min)", fontsize=18, fontname='Arial')
        plt.ylabel("Intensity", fontsize=18, fontname='Arial')
        plt.xticks(fontsize=14, fontname='Arial')
        plt.yticks(fontsize=14, fontname='Arial')

        if output:
            plt.savefig(output, dpi=300, bbox_inches="tight")
            plt.close()
            return None
        else:
            plt.show()
            return eic_rt[np.argmax(eic_int)], np.max(eic_int), eic_scan_idx[np.argmax(eic_int)]


    def plot_all_rois(self, output_path, mz_tol=0.01, rt_range=[0, np.inf], rt_window=None, quality=None):
        """
        Function to plot EIC of all ROIs.
        """

        if output_path[-1] != "/":
            output_path += "/"

        for idx, roi in enumerate(self.rois):

            if quality and roi.quality != quality:
                continue

            if rt_window is not None:
                rt_range = [roi.rt_seq[0] - rt_window, roi.rt_seq[-1] + rt_window]

            # get the eic data
            eic_rt, eic_int, _, eic_scan_idx = self.get_eic_data(roi.mz, mz_tol=mz_tol, rt_range=rt_range)
            idx_start = np.where(eic_scan_idx == roi.scan_idx_seq[0])[0][0]
            idx_end = np.where(eic_scan_idx == roi.scan_idx_seq[-1])[0][0] + 1

            plt.figure(figsize=(9, 3))
            plt.rcParams['font.size'] = 14
            plt.rcParams['font.family'] = 'Arial'
            plt.plot(eic_rt, eic_int, linewidth=0.5, color="black")
            plt.fill_between(eic_rt[idx_start:idx_end], eic_int[idx_start:idx_end], color="black", alpha=0.2)
            plt.axvline(x = roi.rt, color = 'b', linestyle = '--', linewidth=1)
            plt.xlabel("Retention Time (min)", fontsize=18, fontname='Arial')
            plt.ylabel("Intensity", fontsize=18, fontname='Arial')
            plt.xticks(fontsize=14, fontname='Arial')
            plt.yticks(fontsize=14, fontname='Arial')
            plt.text(eic_rt[0], np.max(eic_int)*0.95, "m/z = {:.4f}".format(roi.mz), fontsize=12, fontname='Arial')
            plt.text(eic_rt[0] + (eic_rt[-1]-eic_rt[0])*0.2, np.max(eic_int)*0.95, "Quality = {}".format(roi.quality), fontsize=12, fontname='Arial', color="blue")
            plt.text(eic_rt[0] + (eic_rt[-1]-eic_rt[0])*0.6, np.max(eic_int)*0.95, self.file_name, fontsize=10, fontname='Arial', color="gray")

            file_name = output_path + "roi{}_".format(idx) + str(roi.mz.__round__(4)) + ".png"

            plt.savefig(file_name, dpi=300, bbox_inches="tight")
            plt.close()
    

    def sum_roi_quality(self):
        """
        Function to calculate the sum of all ROI qualities.
        """

        counter1 = 0
        counter2 = 0
        counter3 = 0

        for r in self.rois:
            if r.quality == "good":
                counter1 += 1
            elif r.quality == "short":
                counter2 += 1
            elif r.quality == "bad peak shape":
                counter3 += 1
        
        print("Number of good ROIs: " + str(counter1))
        print("Number of short ROIs: " + str(counter2))
        print("Number of bad peak shape ROIs: " + str(counter3))


class Scan:
    """
    A class that represents a MS scan.
    A MS1 spectrum has properties including:
        scan number, retention time, 
        m/z and intensities.
    A MS2 spectrum has properties including:
        scan number, retention time,
        precursor m/z, product m/z and intensities.
    """

    def __init__(self, level=None, scan=None, rt=None):
        """
        Function to initiate MS1Scan by precursor mz,
        retention time.

        Parameters
        ----------------------------------------------------------
        level: int
            Level of MS scan.
        scan: int
            Scan number.
        rt: float
            Retention time.
        """

        self.level = level
        self.scan = scan
        self.rt = rt

        # for MS1 scans:
        self.mz_seq = None
        self.int_seq = None

        # for MS2 scans:
        self.precursor_mz = None
        self.peaks = None
    

    def add_info_by_level(self, **kwargs):
        """
        Function to add scan information by level.
        """

        if self.level == 1:
            self.mz_seq = kwargs['mz_seq']
            self.int_seq = np.int64(kwargs['int_seq'])

        elif self.level == 2:
            self.precursor_mz = kwargs['precursor_mz']
            self.peaks = kwargs['peaks']


    def show_scan_info(self):
        """
        Function to print a scan's information.

        Parameters
        ----------------------------------------------------------
        scan: MS1Scan or MS2Scan object
            A MS1Scan or MS2Scan object.
        """

        print("Scan number: " + str(self.scan))
        print("Retention time: " + str(self.rt))

        if self.level == 1:
            print("m/z: " + str(np.around(self.mz_seq, decimals=4)))
            print("Intensity: " + str(np.around(self.int_seq, decimals=0)))

        elif self.level == 2:
            # keep 4 decimal places for m/z and 0 decimal place for intensity
            print("Precursor m/z: " + str(np.round(self.precursor_mz, decimals=4)))
            print(self.peaks)
    

    def plot_scan(self, mz_range=None):
        """
        Function to plot a scan.
        
        Parameters
        ----------------------------------------------------------
        """

        if self.level == 1:
            x = self.mz_seq
            y = self.int_seq
        elif self.level == 2:
            x = self.peaks[:, 0]
            y = self.peaks[:, 1]
        
        if mz_range is None:
            mz_range = [np.min(x)-10, np.max(x)+10]
        else:
            y = y[np.logical_and(x > mz_range[0], x < mz_range[1])]
            x = x[np.logical_and(x > mz_range[0], x < mz_range[1])]

        plt.figure(figsize=(10, 3))
        plt.rcParams['font.size'] = 14
        plt.rcParams['font.family'] = 'Arial'
        # plt.scatter(eic_rt, eic_int, color="black")
        plt.vlines(x = x, ymin = 0, ymax = y, color="black", linewidth=1.5)
        plt.hlines(y = 0, xmin = mz_range[0], xmax = mz_range[1], color="black", linewidth=1.5)
        plt.xlabel("m/z, Dalton", fontsize=18, fontname='Arial')
        plt.ylabel("Intensity", fontsize=18, fontname='Arial')
        plt.xticks(fontsize=14, fontname='Arial')
        plt.yticks(fontsize=14, fontname='Arial')
        plt.show()


def _clean_ms2(ms2, offset=1.5, int_threshold=1000):
    """
    A function to clean MS/MS by
    1. Drop ions with m/z > (precursor_mz - offset)   
    2. Drop ions with intensity < 1% of the base peak intensity
    3. Drop ions with intensity lower than threshold
    """
    
    if ms2.peaks.shape[0] > 0:
        ms2.peaks = ms2.peaks[ms2.peaks[:, 0] < ms2.precursor_mz - offset]
    if ms2.peaks.shape[0] > 0:
        ms2.peaks = ms2.peaks[ms2.peaks[:, 1] > 0.01 * np.max(ms2.peaks[:, 1])]
    if ms2.peaks.shape[0] > 0:
        ms2.peaks = ms2.peaks[ms2.peaks[:, 1] > int_threshold]