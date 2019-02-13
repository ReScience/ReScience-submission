# To work with Python 2.7
from __future__ import absolute_import, division, print_function, unicode_literals

import os
import numpy as np
import string
import matplotlib.pyplot as plt
import multiprocessing as mp

from tqdm import tqdm

from tabak import scale_conductance, scale_tabak
from analysis import robustness, tabak_parallel

from uncertainpy.plotting.prettyplot.prettyplot import prettyBar, set_latex_font, set_style, spines_color

dt_default = 0.01
A_default = 3.1415927e-6

dts = [0.05, 0.005, 0.001]
areas = [3.1415927e-9, 3.1415927e-3, 1]
areas_burstiness = [3.1415927e-6, 3.1415927e-9, 3.1415927e-3]


burstiness_factor_reruns = 100

# Equivalent to ../article/figures
figure_folder = os.path.join(os.pardir, "article", "figures")
data_folder = "data"
output_file_area = "robustness_area.txt"
output_file_dt = "robustness_dt.txt"

# Plotting parameters
figure_width = 7.08
titlesize = 13
labelsize = 11
fontsize = 8
fontweight = "medium"
plot_label_weight = "bold"
figure_format = ".eps"
label_x = -0.08
label_y = 1.08
axis_grey = (0.6, 0.6, 0.6)


# Set the random seed to increase reproducability
np.random.seed(10)


# Set default options for plotting
params = {
    "xtick.color": axis_grey,
    "ytick.color": axis_grey,
    "axes.edgecolor": axis_grey,
    "xtick.bottom": True,
    "ytick.left": True,
    "axes.spines.bottom": True,
    "axes.spines.left": True,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1,
    "axes.labelsize": labelsize,
    "font.family": "serif",
}


def burstiness_factors(g_BK, A, dt):
    """
    Calculate the mean and standard deviation of the burstiness factor
    for a given area, timestep and g_BK.

    Parameters
    ----------
    g_BK : float
        The maximal conductance of BK channels, in S/cm^2.
    A : {int, float}
        Area to use in the model, in cm^2.
    dt : {int, float}
        Timestep used in the model, in ms.

    Returns
    ----------
    results : tuple
        A tuple of the results, on the form:
        ((mean G_BK = 0, std G_BK = 0),
         (mean G_BK = 0.5, std G_BK = 0.5),
         (mean G_BK = 1, std G_BK = 1)).
    """
    scaled_parameters = scale_tabak(A)

    burstiness_factors = []

    pool = mp.Pool(processes=mp.cpu_count() - 1)

    parameters = [[g_BK,
                   scaled_parameters["g_K"],
                   scaled_parameters["g_Ca"],
                   scaled_parameters["g_SK"],
                   scaled_parameters["g_l"],
                   A,
                   dt]]*burstiness_factor_reruns

    for burstiness_factor in tqdm(pool.imap(tabak_parallel, parameters),
                                  desc="Running model",
                                  total=burstiness_factor_reruns):

        burstiness_factors.append(burstiness_factor)

    pool.close()

    return np.mean(burstiness_factors), np.std(burstiness_factors)


def burstiness_factor_g_BK_statistics(A, dt):
    """
    Calculate the mean and standard deviation of the burstiness factor
    for a given area and timestep for three values of g_BK.

    Parameters
    ----------
    A : {int, float}
        Area to use in the model, in cm^2.
    dt : {int, float}
        Timestep used in the model, in ms.

    Returns
    ----------
    results : tuple
        A tuple of the results, on the form:
        ((mean G_BK = 0, std G_BK = 0),
         (mean G_BK = 0.5, std G_BK = 0.5),
         (mean G_BK = 1, std G_BK = 1)).
    """
    # G_BK = 0
    g_BK = scale_conductance(0, A)
    mean_0, std_0 = burstiness_factors(g_BK, A, dt)

    # G_BK = 0.5
    g_BK = scale_conductance(0.5, A)
    mean_05, std_05 = burstiness_factors(g_BK, A, dt)

    # G_BK = 1
    g_BK = scale_conductance(1, A)
    mean_1, std_1 = burstiness_factors(g_BK, A, dt)

    results = ((mean_0, std_0), (mean_05, std_05), (mean_1, std_1))

    return results


def burstiness_factor_statistics(areas, dts):
    """
    Calculate the mean and standard deviation of the burstiness factor
    for different areas and timesteps.

    Parameters
    ----------
    areas : list
        Areas to use in the model, in cm^2.
    dts : list
        Timesteps used in the model, in ms.

    Returns
    ----------
    area_results : array
        A collection of the mean and standard deviation of the burstiness factor
        for the different areas.
    dt_results : array
        A collection of the mean and standard deviation of the burstiness factor
        for the different timesteps.
    """
    dt_results = []
    area_results = []

    for i, A in enumerate(areas):
        print("Running for A = {}".format(A))

        results = burstiness_factor_g_BK_statistics(A, dt_default)
        area_results.append(results)

    np.save(os.path.join(data_folder, "area_statistics.npy"), area_results)

    for i, dt in enumerate(dts):
        print("Running for dt = {}".format(dt))

        results = burstiness_factor_g_BK_statistics(A_default, dt)
        dt_results.append(results)

    np.save(os.path.join(data_folder, "timestep_statistics.npy"), dt_results)

    return np.array(area_results), np.array(dt_results)




def generate_data_area(areas):
    """
    Calculate the burstiness factor for the model for various areas.

    Parameters
    ----------
    areas : list
        Areas to use in the model, in cm^2.

    Returns
    -------
    bins : list
        Bins for the binned burstiness factors for each area.
    binned_burstiness_factors : list
        Binned burstiness factors for each area.
    """
    bins = []
    binned_burstiness_factors = []

    with open(os.path.join(data_folder, output_file_area), "w") as output:
        for i, A in enumerate(areas):
            print("Running for A = {}".format(A))

            g_BK = scale_conductance(1, A)

            tmp_bins, tmp_binned_burstiness_factors, bursters, spikers \
                = robustness(g_BK=g_BK,
                             A=A,
                             dt=dt_default)

            bins.append(tmp_bins)
            binned_burstiness_factors.append(tmp_binned_burstiness_factors)

            output.write("A = {}\n".format(A))
            output.write("Spikers = {}\n".format(spikers))
            output.write("Bursters = {}\n\n".format(bursters))

            np.save(os.path.join(data_folder, "bins_{}_area".format(i)), tmp_bins)
            np.save(os.path.join(data_folder, "binned_burstiness_factors_{}_area".format(i)), tmp_binned_burstiness_factors)


    return bins, binned_burstiness_factors



def generate_data_dt(dts):
    """
    Calculate the burstiness factor for the model for various timesteps.

    Parameters
    ----------
    dts : list
        Timesteps to use in the model, in ms.

    Returns
    -------
    bins : list
        Bins for the binned burstiness factors for each timestep.
    binned_burstiness_factors : list
        Binned burstiness factors for each timestep.
    """
    bins = []
    binned_burstiness_factors = []

    with open(os.path.join(data_folder, output_file_dt), "w") as output:
        for i, dt in enumerate(dts):
            print("Running for dt = {}".format(dt))
            g_BK = scale_conductance(1, A_default)
            tmp_bins, tmp_binned_burstiness_factors, bursters, spikers \
                = robustness(g_BK=g_BK, A=A_default, dt=dt)

            bins.append(tmp_bins)
            binned_burstiness_factors.append(tmp_binned_burstiness_factors)

            output.write("dt = {}\n".format(dt))
            output.write("Spikers = {}\n".format(spikers))
            output.write("Bursters = {}\n\n".format(bursters))

            np.save(os.path.join(data_folder, "bins_{}_dt".format(i)), tmp_bins)
            np.save(os.path.join(data_folder, "binned_burstiness_factors_{}_dt".format(i)), tmp_binned_burstiness_factors)


    return bins, binned_burstiness_factors



def load_burstiness_factor_statistics():
    """
    Load results for the burstiness factor calculations.

    Returns
    ----------
    area_results : array
        A collection of the mean and standard deviation of the burstiness factor
        for the different areas.
    dt_results : array
        A collection of the mean and standard deviation of the burstiness factor
        for the different timesteps.
    """

    area_results = np.load(os.path.join(data_folder, "area_statistics.npy"))
    dt_results = np.load(os.path.join(data_folder, "timestep_statistics.npy"))

    return area_results, dt_results


def load_results_area(areas):
    """
    Load results for the robustness calculation using various areas.

    Parameters
    ----------
    areas : list
        Areas to use in the model, in cm^2.

    Returns
    -------
    all_areas : list
        All areas used in the model, the default area is added first.
    bins : list
        Bins for the binned burstiness factors for each area.
    """
    bins = []
    binned_burstiness_factors = []

    for i in range(len(areas)):
        tmp_bins = np.load(os.path.join(data_folder, "bins_{}_area.npy".format(i)))
        tmp_binned_burstiness_factors = np.load(os.path.join(data_folder, "binned_burstiness_factors_{}_area.npy".format(i)))

        bins.append(tmp_bins)
        binned_burstiness_factors.append(tmp_binned_burstiness_factors)

    return bins, binned_burstiness_factors



def load_results_dt(dts):
    """
    Load results for the robustness calculation using various timesteps.

    Parameters
    ----------
    dts : list
        Timersteps to use in the model, in ms.

    Returns
    -------
    all_dts : lists
        Timesteps used in the model, default timestep is added first.
    bins : list
        Bins for the binned burstiness factors for each timersteps.
    """
    bins = []
    binned_burstiness_factors = []

    for i in range(len(dts)):
        tmp_bins = np.load(os.path.join(data_folder, "bins_{}_dt.npy".format(i)))
        tmp_binned_burstiness_factors = np.load(os.path.join(data_folder, "binned_burstiness_factors_{}_dt.npy".format(i)))

        bins.append(tmp_bins)
        binned_burstiness_factors.append(tmp_binned_burstiness_factors)

    return bins, binned_burstiness_factors



def plot_burstiness_factor_statistics(areas, area_results, dts, dt_results):
    """
    Plot the burstines factor mean and standard deviation
    for different timesteps and areas.

    Parameters
    ----------
    areas : list
        Areas to use in the model, in cm^2.
    area_results : array
        A collection of the mean and standard deviation of the burstiness factor
        for the different areas.
    dts : list
        Timesteps used in the model, in ms.
    dt_results : array
        A collection of the mean and standard deviation of the burstiness factor
        for the different timesteps.
    """
    style = "seaborn-darkgrid"
    set_style(style)

    plt.rcParams.update({"axes.titlepad": 8,
                         "font.family": "serif"})

    fig, axes = plt.subplots(nrows=len(areas), ncols=2, figsize=(figure_width, 1.5*figure_width))

    width = 0.2
    xlabels = ["$G_{BK} = 0$", "$G_{BK} = 0.5$", "$G_{BK} = 1$"]

    j = 0
    for i, dt in enumerate(dts):
        ax = axes[i][0]

        increased_titlesize = titlesize + 2
        ax.text(label_x, label_y, string.ascii_uppercase[j], transform=ax.transAxes, fontsize=increased_titlesize, fontweight="bold")

        means = dt_results[i, :, 0]
        stds = dt_results[i, :, 1]

        index = np.arange(1, 4)*width

        prettyBar(means,
                  error=stds,
                  xlabels=xlabels,
                  nr_colors=3,
                  index=index,
                  palette="colorblind",
                  ax=ax,
                  style=style)

        for tick in ax.get_xticklabels():
            tick.set_rotation(-40)

        title = "dt = {} ms".format(dt)
        ax.set_title(title, fontsize=increased_titlesize, fontweight=fontweight)

        ax.set_ylim([0, 1.05])
        ax.set_ylabel("Burstiness factor", fontweight=fontweight, fontsize=titlesize)

        ax.tick_params(axis="both", which="major", labelsize=labelsize, labelcolor="black")

        j += 2

    j = 1
    for i, A in enumerate(areas):
        ax = axes[i][1]

        increased_titlesize = titlesize + 2
        ax.text(label_x, label_y, string.ascii_uppercase[j], transform=ax.transAxes, fontsize=increased_titlesize, fontweight="bold")

        means = area_results[i, :, 0]
        stds = area_results[i, :, 1]

        index = np.arange(1, 4)*width

        prettyBar(means,
                  error=stds,
                  xlabels=xlabels,
                  nr_colors=3,
                  index=index,
                  palette="colorblind",
                  ax=ax,
                  style=style)

        for tick in ax.get_xticklabels():
            tick.set_rotation(-40)

        title = "A = {:.2g} cm$^2$".format(A)
        ax.set_title(title, fontsize=increased_titlesize, fontweight=fontweight)

        ax.set_ylim([0, 1.05])
        ax.tick_params(axis="both", which="major", labelsize=labelsize, labelcolor="black")

        j += 2

    plt.tight_layout()

    plt.savefig(os.path.join(figure_folder, "burstiness_factor" + figure_format))



def plot_areas_dts(areas, bins_areas, binned_burstiness_factors_areas,
                   dts, bins_dts, binned_burstiness_factors_dts):
    """
    Plot the robustness for different timesteps and areas.

    Parameters
    ----------
    areas : list
        Areas used in the model, in cm^2.
    bins_areas : list
        Bins for the binned burstiness factors for each area.
    binned_burstiness_factors_areas : list
        Binned burstiness factors for each area.
    dts : list
        Timesteps used in the model, in ms.
    bins_dts : list
        Bins for the binned burstiness factors for each timestep.
    binned_burstiness_factors_dts : list
        Binned burstiness factors for each timestep.
    """
    plt.rcParams.update(params)

    fig, axes = plt.subplots(nrows=len(dts), ncols=2, figsize=(figure_width, 1.5*figure_width))

    xticks = np.arange(0, 1.1, 0.2)

    j = 0
    for i, dt in enumerate(dts):
        ax = axes[i][0]

        increased_titlesize = titlesize + 2
        ax.text(label_x, label_y, string.ascii_uppercase[j], transform=ax.transAxes, fontsize=increased_titlesize, fontweight="bold")

        ax.bar(bins_dts[i][:-1], binned_burstiness_factors_dts[i], width=(bins_dts[i][1] - bins_dts[i][0]), align="edge")
        title = "dt = {} ms".format(dt)
        ax.set_title(title, fontsize=increased_titlesize, fontweight=fontweight)

        ax.set_ylim([0, 450])
        ax.set_xticks(xticks)
        ax.set_ylabel("Number of models", fontweight=fontweight, fontsize=titlesize)

        ax.tick_params(axis="both", which="major", labelsize=labelsize, labelcolor="black")

        j += 2

    ax.set_xlabel("Burstiness", fontweight=fontweight, fontsize=titlesize)

    j = 1
    for i, A in enumerate(areas):
        ax = axes[i][1]

        increased_titlesize = titlesize + 2
        ax.text(label_x, label_y, string.ascii_uppercase[j], transform=ax.transAxes, fontsize=increased_titlesize, fontweight="bold")

        ax.bar(bins_areas[i][:-1], binned_burstiness_factors_areas[i], width=(bins_areas[i][1] - bins_areas[i][0]), align="edge")
        title = "A = {:.2g} cm$^2$".format(A)
        ax.set_title(title, fontsize=increased_titlesize, fontweight=fontweight)

        ax.set_ylim([0, 450])
        ax.set_xticks(xticks)
        ax.tick_params(axis="both", which="major", labelsize=labelsize, labelcolor="black")

        j += 2

    ax.set_xlabel("Burstiness", fontweight=fontweight, fontsize=titlesize)

    plt.tight_layout()

    plt.savefig(os.path.join(figure_folder, "dt_area" + figure_format))



def figure_burstiness_factor_statistics(areas, dts):
    """
    Perform the calculations and create the figure for the mean and standard
    deviation of the burstiness factor for various areas and timesteps.

    Parameters
    ----------
    areas : list
        Areas used in the model, in cm^2.
    dts : list
        Timesteps used in the model, in ms.
    """
    area_results, dt_results = burstiness_factor_statistics(areas, dts)
    # This can be used to load previously generated results
    # area_results, dt_results = load_burstiness_factor_statistics()

    plot_burstiness_factor_statistics(areas, area_results, dts, dt_results)


def figure_area_dt(areas, dts):
    """
    Perform the calculations and create the figure for the robustness for
    different timesteps and areas.

    Parameters
    ----------
    areas : list
        Areas used in the model, in cm^2.
    dts : list
        Timesteps used in the model, in ms.
    """
    bins_areas, binned_burstiness_factors_areas =  generate_data_area(areas)
    # This can be used to load previously generated results
    # bins_areas, binned_burstiness_factors_areas = load_results_area(areas)

    bins_dts, binned_burstiness_factors_dts = generate_data_dt(dts)

    # This can be used to load previously generated results
    # bins_dts, binned_burstiness_factors_dts = load_results_dt(dts)

    plot_areas_dts(areas, bins_areas, binned_burstiness_factors_areas,
                   dts, bins_dts, binned_burstiness_factors_dts)



if __name__ == "__main__":
    figure_area_dt(areas, dts)
    figure_burstiness_factor_statistics(areas_burstiness, dts)
