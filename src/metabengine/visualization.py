# Author: Hauxu Yu

# A module for data visualization.

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import random
import numpy as np
from . import read_raw_file_to_obj as rfo

def plot_bpcs(data_list=None, file_name_list=None, output=None, autocolor=False):
    """
    A function to plot the base peak chromatograms (overlapped) of a list of data.
    
    Parameters
    ----------
    data_list : list of MSData objects
        A list of data to be plotted.
    """

    if file_name_list is not None:
        data_list = []
        for file_name in file_name_list:
            data_list.append(rfo(file_name))

    if data_list is not None:
        if autocolor:
            color_list = _color_list
        else:
            color_list = ["black"] * len(data_list)

        plt.figure(figsize=(10, 4))
        plt.rcParams['font.size'] = 14
        plt.rcParams['font.family'] = 'Arial'

        for i, d in enumerate(data_list):
            plt.plot(d.ms1_rt_seq, d.bpc_int, color=color_list[i], linewidth=0.5)
            plt.fill_between(d.ms1_rt_seq, d.bpc_int, color=color_list[i], alpha=0.05)
            plt.xlabel("Retention Time (min)", fontsize=18, fontname='Arial')
            plt.ylabel("Intensity", fontsize=18, fontname='Arial')
            plt.xticks(fontsize=14, fontname='Arial')
            plt.yticks(fontsize=14, fontname='Arial')

        if output:
            plt.savefig(output, dpi=600, bbox_inches="tight")
            plt.close()
        else:
            plt.show()


def random_color_generator():
    # set seed
    color = random.choice(list(mcolors.CSS4_COLORS.keys()))
    return color


_color_list = ["red", "blue", "green", "orange", "purple", "brown", "pink", "gray", "olive", "cyan"]


def plot_roi(d, roi, mz_tol=0.01, rt_range=[0, np.inf], rt_window=None, output=False):
    """
    Function to plot EIC of a target m/z.
    """

    if rt_window is not None:
        rt_range = [roi.rt_seq[0] - rt_window, roi.rt_seq[-1] + rt_window]

    # get the eic data
    eic_rt, eic_int, _, eic_scan_idx = d.get_eic_data(roi.mz, mz_tol=mz_tol, rt_range=rt_range)
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
    plt.text(eic_rt[0] + (eic_rt[-1]-eic_rt[0])*0.6, np.max(eic_int)*0.95, d.file_name, fontsize=10, fontname='Arial', color="gray")

    if output:
        plt.savefig(output, dpi=600, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def plot_hist(arr, bins, x_label, y_label):

    plt.figure(figsize=(6, 3))
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Arial'
    plt.hist(arr, bins=bins, color='lightgrey', edgecolor='black', linewidth=0.5)
    plt.xlabel(x_label, fontsize=18, fontname='Arial')
    plt.ylabel(y_label, fontsize=18, fontname='Arial')
    plt.xticks(fontsize=14, fontname='Arial')
    plt.yticks(fontsize=14, fontname='Arial')

    plt.show()


def mirror_ms2(scan1, scan2, output=False):
    """
    Plot a mirror image of two MS2 spectra for comparison.
    """
    if scan1.level == 2 and scan2.level == 2:

        plt.figure(figsize=(10, 3))
        plt.rcParams['font.size'] = 14
        plt.rcParams['font.family'] = 'Arial'
        # plot precursor
        plt.vlines(x = scan1.precursor_mz, ymin = 0, ymax = 1, color="cornflowerblue", linewidth=1.5, linestyles='dashed')
        plt.vlines(x = scan2.precursor_mz, ymin = 0, ymax = -1, color="lightcoral", linewidth=1.5, linestyles='dashed')

        # plot fragment ions
        plt.vlines(x = scan1.peaks[:, 0], ymin = 0, ymax = scan1.peaks[:, 1] / np.max(scan1.peaks[:, 1]), color="blue", linewidth=1.5)
        plt.vlines(x = scan2.peaks[:, 0], ymin = 0, ymax = -scan2.peaks[:, 1] / np.max(scan2.peaks[:, 1]), color="red", linewidth=1.5)

        # plot zero line
        plt.hlines(y = 0, xmin = 0, xmax = max([scan1.precursor_mz, scan2.precursor_mz])*1.1, color="black", linewidth=1.5)
        plt.xlabel("m/z, Dalton", fontsize=18, fontname='Arial')
        plt.ylabel("Intensity", fontsize=18, fontname='Arial')
        plt.xticks(fontsize=14, fontname='Arial')
        plt.yticks(fontsize=14, fontname='Arial')

        if output:
            plt.savefig(output, dpi=600, bbox_inches="tight")
            plt.close()
        else:
            plt.show()