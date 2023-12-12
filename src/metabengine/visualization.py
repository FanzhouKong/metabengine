# Author: Hauxu Yu

# A module for data visualization.

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import random
import numpy as np
import networkx as nx
from .annotation import feature_to_feature_search

def plot_bpcs(data_list=None, output=None, autocolor=False):
    """
    A function to plot the base peak chromatograms (overlapped) of a list of data.
    
    Parameters
    ----------
    data_list : list of MSData objects
        A list of data to be plotted.
    """

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


def plot_roi(d, roi, mz_tol=0.01, rt_range=[0, np.inf], rt_window=None, output=False, break_scan=None, label_quality=True):
    """
    Function to plot EIC of a target m/z.
    """

    if rt_window is not None:
        rt_range = [roi.rt_seq[0] - rt_window, roi.rt_seq[-1] + rt_window]
    
    # get the eic data
    eic_rt, eic_int, _, eic_scan_idx = d.get_eic_data(roi.mz, mz_tol=mz_tol, rt_range=rt_range)
    idx_start = np.where(eic_scan_idx == roi.scan_idx_seq[0])[0][0]
    idx_end = np.where(eic_scan_idx == roi.scan_idx_seq[-1])[0][0] + 1

    if break_scan is not None:
        idx_middle = np.where(eic_scan_idx == break_scan)[0][0]

    plt.figure(figsize=(9, 3))
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Arial'
    plt.plot(eic_rt, eic_int, linewidth=0.5, color="black")

    if break_scan is not None:
        plt.fill_between(eic_rt[idx_start:(idx_middle+1)], eic_int[idx_start:(idx_middle+1)], color="blue", alpha=0.2)
        plt.fill_between(eic_rt[idx_middle:idx_end], eic_int[idx_middle:idx_end], color="red", alpha=0.2)
    else:
        plt.fill_between(eic_rt[idx_start:idx_end], eic_int[idx_start:idx_end], color="black", alpha=0.2)
    plt.axvline(x = roi.rt, color = 'b', linestyle = '--', linewidth=1)
    plt.xlabel("Retention Time (min)", fontsize=18, fontname='Arial')
    plt.ylabel("Intensity", fontsize=18, fontname='Arial')
    plt.xticks(fontsize=14, fontname='Arial')
    plt.yticks(fontsize=14, fontname='Arial')
    plt.text(eic_rt[0], np.max(eic_int)*0.95, "m/z = {:.4f}".format(roi.mz), fontsize=12, fontname='Arial')
    if label_quality:
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


def mirror_ms2_from_scans(scan1, scan2, output=False):
    """
    Plot a mirror image of two MS2 spectra for comparison.
    """

    if scan1.level == 2 and scan2.level == 2:
        mirror_ms2(precursor_mz1=scan1.precursor_mz, precursor_mz2=scan2.precursor_mz, peaks1=scan1.peaks, peaks2=scan2.peaks, output=output)


def mirror_ms2(precursor_mz1, precursor_mz2, peaks1, peaks2, output=False):

    plt.figure(figsize=(10, 3))
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Arial'
    # plot precursor
    plt.vlines(x = precursor_mz1, ymin = 0, ymax = 1, color="cornflowerblue", linewidth=1.5, linestyles='dashed')
    plt.vlines(x = precursor_mz2, ymin = 0, ymax = -1, color="lightcoral", linewidth=1.5, linestyles='dashed')

    # plot fragment ions
    plt.vlines(x = peaks1[:, 0], ymin = 0, ymax = peaks1[:, 1] / np.max(peaks1[:, 1]), color="blue", linewidth=1.5)
    plt.vlines(x = peaks2[:, 0], ymin = 0, ymax = -peaks2[:, 1] / np.max(peaks2[:, 1]), color="red", linewidth=1.5)

    # plot zero line
    plt.hlines(y = 0, xmin = 0, xmax = max([precursor_mz1, precursor_mz2])*1.1, color="black", linewidth=1.5)
    plt.xlabel("m/z, Dalton", fontsize=18, fontname='Arial')
    plt.ylabel("Intensity", fontsize=18, fontname='Arial')
    plt.xticks(fontsize=14, fontname='Arial')
    plt.yticks(fontsize=14, fontname='Arial')

    if output:
        plt.savefig(output, dpi=600, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def mirror_ms2_db(f, output=False):

    precursor_mz1 = f.mz
    precursor_mz2 = f.matched_precursor_mz
    peaks1 = f.best_ms2.peaks
    peaks2 = f.matched_peaks

    plt.figure(figsize=(10, 3))
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Arial'
    # plot precursor
    plt.vlines(x = precursor_mz1, ymin = 0, ymax = 1, color="cornflowerblue", linewidth=1.5, linestyles='dashed')
    plt.vlines(x = precursor_mz2, ymin = 0, ymax = -1, color="lightcoral", linewidth=1.5, linestyles='dashed')

    # plot fragment ions
    plt.vlines(x = peaks1[:, 0], ymin = 0, ymax = peaks1[:, 1] / np.max(peaks1[:, 1]), color="blue", linewidth=1.5)
    plt.vlines(x = peaks2[:, 0], ymin = 0, ymax = -peaks2[:, 1] / np.max(peaks2[:, 1]), color="red", linewidth=1.5)

    xmax = max([precursor_mz1, precursor_mz2])*1.2
    # plot zero line
    plt.hlines(y = 0, xmin = 0, xmax = xmax, color="black", linewidth=1.5)
    plt.xlabel("m/z, Dalton", fontsize=18, fontname='Arial')
    plt.ylabel("Intensity", fontsize=18, fontname='Arial')
    plt.xticks(fontsize=14, fontname='Arial')
    plt.yticks(fontsize=14, fontname='Arial')

    # note name and similarity score
    plt.text(xmax*0.9, 0.9, "Experiment", fontsize=12, fontname='Arial', color="grey")
    plt.text(xmax*0.9, -0.9, "Database", fontsize=12, fontname='Arial', color="grey")
    plt.text(0, 0.9, "similarity = {:.3f}".format(f.similarity), fontsize=12, fontname='Arial', color="blue")
    plt.text(0, -0.95, f.annotation.lower(), fontsize=12, fontname='Arial', color="black")

    if output:
        plt.savefig(output, dpi=600, bbox_inches="tight")
        plt.close()
    else:
        plt.show()



def plot_network(feature_list, annotation_type="hybrid_and_identity", feature_quality="all", output=False):
    """
    A function to plot a network graph.

    Parameters
    ----------
    feature_list : list of Feature objects
        A list of features to be plotted.
    annotation_type : str
        Type of annotation to be plotted. Default is "all".
        "all" - all the features with MS2 spectra.
        "hybrid_and_identity" - features with identity and hybrid annotation.
        "identity_only" - only features with identity annotation.
        "hybrid_only" - only features with hybrid annotation.
    feature_quality : str
        Quality of features to be plotted. Default is "all".
        "all" - all the features.
        "good" - only good features (quality=="good").
        "bad" - only bad features (quality=="bad peak shape").
    """

    # prepare feature list
    selected_features = _prepare_feature_list_for_network(feature_list, annotation_type, feature_quality)

    df = feature_to_feature_search(selected_features)

    hybrid_features = [f for f in selected_features if f.annotation_mode == "hybrid_search"]

    identity_search_names = [f.annotation for f in selected_features if f.annotation_mode == "identity_search"]

    if len(hybrid_features) > 0:
        for f in hybrid_features:
            if f.annotation in identity_search_names:
                df.loc[len(df)] = [f.network_name, f.annotation, f.similarity, f.id, "DB"]
            else:
                df.loc[len(df)] = [f.network_name, "DB_"+f.annotation, f.similarity, f.id, "DB"]

    # Create a new graph
    G = nx.Graph()

    # prepare nodes
    nodes = df["feature_name_1"].tolist() + df["feature_name_2"].tolist()
    nodes = list(set(nodes))

    # prepare edges
    edges = []
    for i in range(len(df)):
        edges.append((df["feature_name_1"][i], df["feature_name_2"][i]))

    # define node colors: identity - red, hybrid - "#8BABD3", database - black, unknown - gray
    node_color = []
    for n in nodes:
        if n.startswith("hybrid"):
            node_color.append("#8BABD3")
        elif n.startswith("unknown"):
            node_color.append("gray")
        elif n.startswith("DB"):
            node_color.append("black")
        else:
            node_color.append("red")

    # define edge colors as a gradient of similarity
    edge_color = _edge_color_gradient(df["similarity"])

    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    # make plot
    pos = nx.spring_layout(G, iterations=25)  # positions for all nodes
    nx.draw_networkx_edges(G, pos, edge_color=edge_color, width=1)
    nx.draw_networkx_nodes(G, pos, node_color=node_color, node_size=50, alpha=0.5, edgecolors="none")

    if output:
        plt.savefig(output, dpi=1000)
        plt.close()
    else:
        plt.show()


def _prepare_feature_list_for_network(feature_list, annotation_type="hybrid_and_identity", feature_quality="all"):
    """
    A function to prepare the feature list for plotting.
    
    Parameters
    ----------
    feature_list : list of Feature objects
        A list of features to be plotted.
    annotation_type : str
        Type of annotation to be plotted. Default is "all".
        "all" - all the features with MS2 spectra.
        "hybrid_and_identity" - features with identity and hybrid annotation.
        "identity_only" - only features with identity annotation.
        "hybrid_only" - only features with hybrid annotation.
    feature_quality : str
        Quality of features to be plotted. Default is "all".
        "all" - all the features.
        "good" - only good features (quality=="good").
        "bad" - only bad features (quality=="bad peak shape").

    Returns
    -------
    selected_features : list of Feature objects
        A list of features to be plotted.
    """
    
    selected_features = [f for f in feature_list if f.best_ms2 is not None]

    if annotation_type == "all":
        selected_features = feature_list
    elif annotation_type == "hybrid_and_identity":
        selected_features = [f for f in feature_list if f.annotation_mode in ["identity_search", "hybrid_search"]]
    elif annotation_type == "identity_only":
        selected_features = [f for f in feature_list if f.annotation_mode == "identity_search"]
    elif annotation_type == "hybrid_only":
        selected_features = [f for f in feature_list if f.annotation_mode == "hybrid_search"]
    else:
        raise ValueError("Invalid annotation_type: {}".format(annotation_type))
    
    if feature_quality == "all":
        pass
    elif feature_quality == "good":
        selected_features = [f for f in selected_features if f.quality == "good"]
    elif feature_quality == "bad":
        selected_features = [f for f in selected_features if f.quality == "bad peak shape"]
    else:
        raise ValueError("Invalid feature_quality: {}".format(feature_quality))
    
    for f in selected_features:
        if f.annotation_mode == "identity_search":
            f.network_name = "identity_{}".format(f.id) + "_" + f.annotation
        elif f.annotation_mode == "hybrid_search":
            f.network_name = "hybrid_{}".format(f.id)
        else:
            f.network_name = "unknown_{}".format(f.id)
    
    return selected_features


def _edge_color_gradient(similarity_array, color_1="lightgrey", color_2="black"):
    """
    A function to generate a list of edge colors based on the similarity scores.
    """

    colors = []

    similarity_array = np.array(similarity_array)

    similarity_array = (np.max(similarity_array) - similarity_array) / (np.max(similarity_array) - np.min(similarity_array))

    for s in similarity_array:
        colors.append(_color_gradient(s, color_1, color_2))

    return colors


def _color_gradient(s, color_1, color_2):
    """
    A function to generate a color based on the similarity score.
    """

    color_1 = mcolors.to_rgb(color_1)
    color_2 = mcolors.to_rgb(color_2)

    r = color_1[0] + s * (color_2[0] - color_1[0])
    g = color_1[1] + s * (color_2[1] - color_1[1])
    b = color_1[2] + s * (color_2[2] - color_1[2])

    return (r, g, b)

