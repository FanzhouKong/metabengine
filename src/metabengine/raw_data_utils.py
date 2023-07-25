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
    """


    def __init__(self):
        """
        Function to initiate MSData.
        ----------------------------------------------------------
        """

        self.scans = []   # A list of MS scans
        self.ms1_rt_seq = np.array([])  # Retention times of all MS1 scans
        self.bpc_int = np.array([]) # Intensity of the BPC
        self.rois = []  # A list of ROIs
        self.params = None  # A Params object


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

        if os.path.isfile(file_name):
            # get extension from file name
            ext = os.path.splitext(file_name)[1]

            self.file_name = file_name

            if ext.lower() == ".mzml":
                with mzml.MzML(file_name) as reader:
                    self.extract_scan_mzml(reader, params=params)
            elif ext.lower() == ".mzxml":
                with mzxml.MzXML(file_name) as reader:
                    self.extract_scan_mzxml(reader, params=params)
            else:
                print("Unsupported raw data format. Raw data must be in mzML or mzXML.")
        else:
            print("File does not exist.")


    def extract_scan_mzml(self, spectra, params):
        """
        Function to extract all scans and convert them to Scan objects.

        Parameters
        ----------------------------------------------------------
        spectra: pyteomics object
            An iteratable object that contains all MS1 and MS2 scans.
        params: Params object
            A Params object that contains the parameters.
        """

        idx = 0     # Scan number
        self.ms1_idx = np.array([], dtype=int)   # MS1 scan index
        self.ms2_idx = np.array([], dtype=int)   # MS2 scan index

        # Iterate over all scans
        for spec in spectra:
            # Get the retention time and convert to minute
            try:
                rt = spec['scanList']['scan'][0]['scan start time'] / 60
            except:
                rt = spec['scanList']['scan'][0]['scan time'] / 60

            # Check if the retention time is within the range
            if params.rt_range[0] < rt < params.rt_range[1]:
                if spec['ms level'] == 1:
                    temp_scan = Scan(level=1, scan=idx, rt=rt)
                    mz_array = spec['m/z array']
                    int_array = spec['intensity array']

                    temp_scan.add_info_by_level(mz_seq=mz_array, int_seq=int_array)
                    self.ms1_idx = np.append(self.ms1_idx, idx)

                    # update base peak chromatogram
                    self.bpc_int = np.append(self.bpc_int, np.max(spec['intensity array']))
                    self.ms1_rt_seq = np.append(self.ms1_rt_seq, rt)

                elif spec['ms level'] == 2:
                    temp_scan = Scan(level=2, scan=idx, rt=rt)
                    precs_mz = spec['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                    prod_int_seq = spec['intensity array']
                    prod_mz_seq = spec['m/z array']
                    if len(prod_int_seq) > 0:
                        prod_mz_seq = prod_mz_seq[prod_int_seq > np.max(prod_int_seq)*0.01]
                        prod_int_seq = prod_int_seq[prod_int_seq > np.max(prod_int_seq)*0.01]
                    temp_scan.add_info_by_level(precs_mz=precs_mz, prod_mz_seq=prod_mz_seq, prod_int_seq=prod_int_seq)
                    self.ms2_idx = np.append(self.ms2_idx, idx)
                
                self.scans.append(temp_scan)
                idx += 1

        # print the number of extracted ms1 and ms2 scans
        print("Number of extracted MS1 scans: " + str(len(self.ms1_idx)))
        print("Number of extracted MS2 scans: " + str(len(self.ms2_idx)))


    def extract_scan_mzxml(self, spectra, params):
        """
        Function to extract all scans and convert them to Scan objects.

        Parameters
        ----------------------------------------------------------
        spectra: pyteomics object
            An iteratable object that contains all MS1 and MS2 scans.
        params: Params object
            A Params object that contains the parameters.
        """

        idx = 0     # Scan number
        self.ms1_idx = np.array([], dtype=int)   # MS1 scan index
        self.ms2_idx = np.array([], dtype=int)   # MS2 scan index

        # Iterate over all scans
        for spec in spectra:
            # Get the retention time and convert to minute
            rt = spec["retentionTime"]    # retention time of mzXML is in minute

            # Check if the retention time is within the range
            if params.rt_range[0] < rt < params.rt_range[1]:
                if spec['msLevel'] == 1:
                    temp_scan = Scan(level=1, scan=idx, rt=rt)
                    mz_array = spec['m/z array']
                    int_array = spec['intensity array']

                    temp_scan.add_info_by_level(mz_seq=mz_array, int_seq=int_array)
                    self.ms1_idx = np.append(self.ms1_idx, idx)

                    # update base peak chromatogram
                    self.bpc_int = np.append(self.bpc_int, np.max(spec['intensity array']))
                    self.ms1_rt_seq = np.append(self.ms1_rt_seq, rt)

                elif spec['msLevel'] == 2:
                    temp_scan = Scan(level=2, scan=idx, rt=rt)
                    precs_mz = spec['precursorMz'][0]['precursorMz']
                    prod_int_seq = spec['intensity array']
                    prod_mz_seq = spec['m/z array']
                    if len(prod_int_seq) > 0:
                        prod_mz_seq = prod_mz_seq[prod_int_seq > np.max(prod_int_seq)*0.01]
                        prod_int_seq = prod_int_seq[prod_int_seq > np.max(prod_int_seq)*0.01]
                    temp_scan.add_info_by_level(precs_mz=precs_mz, prod_mz_seq=prod_mz_seq, prod_int_seq=prod_int_seq)
                    self.ms2_idx = np.append(self.ms2_idx, idx)
                
                self.scans.append(temp_scan)
                idx += 1

        # print the number of extracted ms1 and ms2 scans
        print("Number of extracted MS1 scans: " + str(len(self.ms1_idx)))
        print("Number of extracted MS2 scans: " + str(len(self.ms2_idx)))

    
    def drop_ion_by_int(self, params):
        """
        Function to drop ions by intensity.

        Parameters
        ----------------------------------------------------------
        tol: int
            Intensity tolerance.
        """

        for idx in self.ms1_idx:
            self.scans[idx].mz_seq = self.scans[idx].mz_seq[self.scans[idx].int_seq > params.int_tol]
            self.scans[idx].int_seq = self.scans[idx].int_seq[self.scans[idx].int_seq > params.int_tol]
    

    def find_rois(self, params):
        """
        Function to find ROI in MS1 scans.

        Parameters
        ----------------------------------------------------------
        params: Params object
            A Params object that contains the parameters.
        """

        self.rois = peak_detect.roi_finder(self, params)

        print("Number of regular ROIs: " + str(len(self.rois)))
    

    def cut_rois(self, params):
        """
        Function to cut ROI into smaller pieces.

        Parameters
        ----------------------------------------------------------
        params: Params object
            A Params object that contains the parameters.
        """

        small_rois = []
        to_be_removed = []

        for i, roi in enumerate(self.rois):
            positions = peak_detect.find_roi_cut(roi, params)
            
            if positions is not None:

                # append each item in a list to small_rois
                small_rois.extend(peak_detect.roi_cutter(roi, positions))

                to_be_removed.append(i)
        
        # remove the original rois
        for i in sorted(to_be_removed, reverse=True):
            del self.rois[i]

        # append the small rois to the original rois
        self.rois.extend(small_rois)

        # print("Number of regular ROIs: " + str(len(self.rois)))
        # print("Number of short ROIs: " + str(len(self.rois_short)))

    def process_rois(self, params):
        """
        Function to process ROIs.

        Parameters
        ----------------------------------------------------------
        params: Params object
            A Params object that contains the parameters.
        """

        # 1. sort ROI by m/z
        self.rois.sort(key=lambda x: x.mz)

        for roi in self.rois:
            # 1. find roi quality by length
            if len(roi.mz_seq) >= params.min_ion_num:
                roi.quality = 'good'
            else:
                roi.quality = 'short'
            
            # 2. find the best MS2
            roi.find_best_ms2()


    def find_roi_quality_by_length(self, params):
        """
        Function to find the quality of each ROI by its length.

        Parameters
        ----------------------------------------------------------
        params: Params object
            A Params object that contains the parameters.
        """
        
        for roi in self.rois:
            if len(roi.mz_seq) >= params.min_ion_num:
                roi.quality = 'good'
            else:
                roi.quality = 'short'
    

    def sort_rois_by_rt(self):
        """
        Function to sort rois by rt.
        """

        self.rois.sort(key=lambda x: x.rt)
    

    def clean_rois_single_ms2(self):
        """
        Function to clean rois by selecting the best MS2 scan (highest total intensity).
        """

        for roi in self.rois:
            roi.keep_single_ms2()
    

    def plot_bpc(self, output=False):
        """
        Function to plot base peak chromatogram.
        """

        plt.figure(figsize=(10, 3))
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
        
    
    def output_roi_report(self, path):
        """
        Function to generate a report for rois in csv format.
        """

        result = []
        for roi in self.rois:
            roi.sum_roi()
            temp = np.array([roi.mz, roi.rt, roi.roi_mz_error(), roi.roi_length(), roi.peak_area, roi.peak_height, roi.peak_height_by_ave])
            result.append(temp)

        # convert the result to a numpy array
        result = np.array(result)

        # convert result to a pandas dataframe
        df = pd.DataFrame(result, columns=['mz', 'rt', 'mz_error', 'length', 'peak_area', 'peak_height', 'peak_height_by_ave'])

        # save the dataframe to csv file
        df.to_csv(path, index=False)
    

    def get_eic_data(self, mz, rt, mz_tol=0.005):
        rt_start = rt - 1.0
        rt_end = rt + 1.0
        
        eic_int_seq = np.array([])
        eic_rt_seq = np.array([])

        for scan in self.scans:
            if scan.rt > rt_start:
                if scan.level == 1:
                    eic_rt_seq = np.append(eic_rt_seq, scan.rt)
                    mz_diff = np.abs(scan.mz_seq - mz)
                    if np.min(mz_diff) < mz_tol:
                        eic_int_seq = np.append(eic_int_seq, scan.int_seq[np.argmin(mz_diff)])
                    else:
                        eic_int_seq = np.append(eic_int_seq, 0)
            if scan.rt > rt_end:
                break
        
        return [eic_int_seq, eic_rt_seq]


    def find_roi_by_mz(self, mz, mz_tol=0.005):
        rois = []
        for roi in self.rois:
            if np.abs(roi.mz - mz) < mz_tol:
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
        eic_int = np.array([])
        eic_rt = np.array([])
        for scan in self.scans:
            if scan.rt > rt_range[0]:
                if scan.level == 1:
                    mz_diff = np.abs(scan.mz_seq - target_mz)
                    if np.min(mz_diff) < mz_tol:
                        eic_rt = np.append(eic_rt, scan.rt)
                        eic_int = np.append(eic_int, scan.int_seq[np.argmin(mz_diff)])
                    else:
                        eic_rt = np.append(eic_rt, scan.rt)
                        eic_int = np.append(eic_int, 0)
            if scan.rt > rt_range[1]:
                break

        plt.figure(figsize=(10, 3))
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

    def __init__(self, level, scan, rt):
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
        self.precs_mz = None
        self.prod_mz_seq = None
        self.prod_int_seq = None
    

    def add_info_by_level(self, **kwargs):
        """
        Function to add scan information by level.
        """

        if self.level == 1:
            self.mz_seq = kwargs['mz_seq']
            self.int_seq = np.int64(kwargs['int_seq'])

        elif self.level == 2:
            self.precs_mz = kwargs['precs_mz']
            self.prod_mz_seq = kwargs['prod_mz_seq']
            self.prod_int_seq = np.int64(kwargs['prod_int_seq'])


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
            print("Precursor m/z: " + str(np.round(self.precs_mz, decimals=4)))
            print("Product m/z: " + str(np.around(self.prod_mz_seq, decimals=4)))
            print("Product intensity: " + str(np.around(self.prod_int_seq, decimals=0)))



def plot_eic(data, mz_seq, rt_seq, mz_tol=0.005):
    """
    A function to plot the EICs

    Parameters
    ----------
    data : MSData object
    """

    for i in range(len(mz_seq)):
        mz = mz_seq[i]
        rt = rt_seq[i]

        # get EIC data
        eic = data.get_eic_data(mz, rt, mz_tol=mz_tol)

        # plot EIC
        plt.plot(eic[1], eic[0], linewidth=1, color="black")
        plt.xlabel("Retention Time (min)", fontsize=18, fontname='Arial')
        plt.ylabel("Intensity", fontsize=18, fontname='Arial')
        plt.xticks(fontsize=14, fontname='Arial')
        plt.yticks(fontsize=14, fontname='Arial')
        # save to png
        plt.savefig("eic_" + str(mz) + "_" + str(rt) + ".png", dpi=300, bbox_inches="tight")
        plt.close()
