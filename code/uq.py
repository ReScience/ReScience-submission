# To work with Python 2.7
from __future__ import absolute_import, division, print_function, unicode_literals

import string
import os

import chaospy as cp
import matplotlib.pyplot as plt
import numpy as np
import uncertainpy as un

from mpl_toolkits.mplot3d import Axes3D, art3d
from matplotlib import colors
from tqdm import tqdm

from tabak import scale_conductance
from tabak import G_l, G_K, G_SK, G_Ca
from burstiness import burstiness_factor, duration, min_spike_amplitude

from uncertainpy.plotting.prettyplot.prettyplot import prettyBar, set_latex_font, set_style, spines_color


# Plotting parameters
figure_width = 7.08
titlesize = 12
labelsize = 10
fontsize = 8
fontweight = "medium"
plot_label_weight = "bold"
label_x = -0.16
label_y = 1.08
figure_format = ".eps"

# Simulation parameters
discard = 10000                # ms
simulation_time = 60000        # ms
noise_amplitude = 0            # mV

# Equivalent to ../article/figures
figure_folder = os.path.join(os.pardir, "article", "figures")
data_folder = "data"

# Various parameters
polynomial_order = 8
nr_points_exploration = 60
features_to_run = ["spike_rate", "average_AP_overshoot", "average_AHP_depth",
                   "burstiness_factor", "average_duration"]


def plot_sobol(data, filename):
    """
    Plot the total-order Sobol indices from the uncertainty quantification.

    Parameters
    ----------
    data : uncertainpy.Data object
        The data from the uncertainty quantification.
    filename : str
        Filename of the figure.
    """

    feature_titles = {"spike_rate": "Event rate",
                      "average_AP_overshoot": "Average event peak",
                      "average_AHP_depth": "Average AHP depth" ,
                      "burstiness_factor": "Burstiness factor",
                      "average_duration": "Average duration"}

    nr_plots = len(features_to_run)
    grid_size = np.ceil(np.sqrt(nr_plots))
    grid_x_size = int(grid_size)
    grid_y_size = int(np.ceil(nr_plots/float(grid_x_size)))


    style = "seaborn-darkgrid"

    set_style(style)

    plt.rcParams.update({"axes.titlepad": 8,
                         "font.family": "serif"})
    fig, axes = plt.subplots(nrows=grid_y_size,
                             ncols=grid_x_size,
                             squeeze=False,
                             figsize=(figure_width, figure_width*0.8))


    # Add a larger subplot to use to set a common xlabel and ylabel
    set_style("seaborn-white")
    ax = fig.add_subplot(111, zorder=-10)
    spines_color(ax, edges={"top": "None", "bottom": "None",
                            "right": "None", "left": "None"})
    ax.tick_params(top=False, bottom=False, left=False, right=False, labelsize=labelsize,
                   labelbottom=False, labelleft=False)
    ax.set_ylabel("Total-order Sobol indices", labelpad=35, fontsize=titlesize, fontweight=fontweight)


    width = 0.2
    index = np.arange(1, len(data.uncertain_parameters)+1)*width


    latex_labels = {"g_K": r"$g_\mathrm{K}$",
                    "g_Ca": r"$g_\mathrm{Ca}$",
                    "g_SK": r"$g_\mathrm{SK}$",
                    "g_Na": r"$g_\mathrm{Na}$",
                    "g_l": r"$g_\mathrm{l}$",
                    "g_BK": r"$g_\mathrm{BK}$"}

    xlabels = []
    for label in data.uncertain_parameters:
        xlabels.append(latex_labels[label])


    for i in range(0, grid_x_size*grid_y_size):
        nx = i % grid_x_size
        ny = int(np.floor(i/float(grid_x_size)))

        ax = axes[ny][nx]

        if i < nr_plots:
            title = feature_titles[features_to_run[i]]

            sensitivity = data[features_to_run[i]].sobol_total_average
            mean = data[features_to_run[i]].mean
            std = np.sqrt(data[features_to_run[i]].variance)
            unit = data[features_to_run[i]].labels

            if unit and "(" in unit[0]:
                unit = unit[0].split("(")[-1].strip(")")
            else:
                unit = ""


            # Convert spike rate from 1/ms to Hz (1/s)
            if features_to_run[i] == "spike_rate":
                mean *= 1000
                std *= 1000
                unit = "Hz"

            prettyBar(sensitivity,
                      xlabels=xlabels,
                      nr_colors=len(data.uncertain_parameters),
                      index=index,
                      palette="colorblind",
                      ax=ax,
                      style=style)

            for tick in ax.get_xticklabels():
                tick.set_rotation(-40)

            ax.set_ylim([0, 1.15])
            ax.set_title(title, fontsize=titlesize)
            ax.text(label_x,
                    label_y,
                    string.ascii_uppercase[i],
                    transform=ax.transAxes,
                    fontsize=titlesize,
                    fontweight=plot_label_weight)

            ax.text(0.1,
                    0.9,
                    "Mean = {mean:.2{c}} {unit}".format(mean=mean,
                                                        c="e" if abs(mean) < 1e-2 else "f",
                                                        unit=unit),
                    transform=ax.transAxes,
                    fontsize=labelsize,
                    fontweight=fontweight)

            ax.text(0.1,
                    0.8,
                    "Std. = {std:.2{c}} {unit}".format(std=std,
                                                       c="e" if abs(mean) < 1e-2 else "f",
                                                       unit=unit),
                    transform=ax.transAxes,
                    fontsize=labelsize,
                    fontweight=fontweight)

        else:
            ax.axis("off")

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.13, left=0.1)
    plt.savefig(filename + figure_format)

    # Reset matplotlib defaults
    plt.rcdefaults()




def plot_wireframe(ax, x, y, z, cmap, norm=None, **kwargs):
    """
    Plot a 3D wireframe with a colormap along the lines.

    Parameters
    ----------
    ax: Matplotlib Axis object
        Axis object where to plot
    x : array_like
        2D data values to plot.
    y : array_like
        2D data values to plot.
    z : array_like
        2D data values to plot.
    cmap: Matplotlib Colormap
        Colormap of the wires.
    norm : {Matplotlib Normalize, None}, optional
        Normalization for the colormap.
    **kwargs
        Other arguments are forwarded to ax.plot_wireframe().

    Returns
    -------
    wire: Matplotlib Line3DCollection
        The Line3DCollection object from ax.plot_wireframe.

    Notes
    -----
    Solution to create colored wireframe plot from:
    https://stackoverflow.com/questions/24909256/how-to-obtain-3d-colored-surface-via-python
    """

    if np.shape(x) != np.shape(y):
        raise ValueError("The shape of x and y must be the same")

    wire = ax.plot_wireframe(x, y, z, **kwargs)

    # # Retrive data from internal storage of plot_wireframe, then delete it
    nx, ny, _  = np.shape(wire._segments3d)
    wire_x = np.array(wire._segments3d)[:, :, 0].ravel()
    wire_y = np.array(wire._segments3d)[:, :, 1].ravel()
    wire_z = np.array(wire._segments3d)[:, :, 2].ravel()
    wire.remove()

    # create data for a LineCollection
    wire_x1 = np.vstack([wire_x, np.roll(wire_x, 1)])
    wire_y1 = np.vstack([wire_y, np.roll(wire_y, 1)])
    wire_z1 = np.vstack([wire_z, np.roll(wire_z, 1)])

    to_delete = np.arange(0, nx*ny, ny)
    wire_x1 = np.delete(wire_x1, to_delete, axis=1)
    wire_y1 = np.delete(wire_y1, to_delete, axis=1)
    wire_z1 = np.delete(wire_z1, to_delete, axis=1)
    scalars = np.delete(wire_z, to_delete)

    segs = [list(zip(xl, yl, zl)) for xl, yl, zl in zip(wire_x1.T, wire_y1.T, wire_z1.T)]

    # Plots the wireframe by a line3DCollection
    my_wire = art3d.Line3DCollection(segs, cmap=cmap, norm=norm)
    my_wire.set_array(scalars)
    ax.add_collection(my_wire)

    return wire



def perform_exploration():
    """
    Perform an limited parameter exploration of the average duration and
    burstiness of the model, varying g_SK and g_BK, as well as g_K and g_BK.

    Returns
    -------
    original_g_BKs : array
        The g_BK values used in the exploration, not scaled.
    original_g_SKs : array
        The g_SK values used in the exploration, not scaled.
    original_g_Ks : array
        The g_K values used in the exploration, not scaled.
    event_durations_SK : array
        The event durations when g_SK and g_BK is changed.
    event_durations_K : array
        The event durations when g_K and g_BK is changed.
    """
    model = un.NeuronModel(file="tabak.py", name="tabak")

    original_g_BKs = np.linspace(0, 1, nr_points_exploration)
    g_BKs = scale_conductance(original_g_BKs)

    # g_SK
    original_g_SKs = np.linspace(1, 3, nr_points_exploration)
    g_SKs = scale_conductance(original_g_SKs)

    event_durations_SK = np.zeros((len(original_g_BKs), len(original_g_SKs)))
    for i, g_BK in enumerate(tqdm(g_BKs, desc="Varying g_BK")):
        for j, g_SK in enumerate(tqdm(g_SKs, desc="Varying g_SK")):
            time, voltage, info = model.run(g_BK=g_BK,
                                            g_SK=g_SK,
                                            discard=discard,
                                            noise_amplitude=noise_amplitude,
                                            simulation_time=simulation_time)

            tmp_duration = np.mean(duration(time, voltage))

            if np.isnan(tmp_duration):
                tmp_duration = -1

            event_durations_SK[i, j] = tmp_duration


    original_g_Ks = np.linspace(1.5, 4.5, nr_points_exploration)
    g_Ks = scale_conductance(original_g_Ks)

    event_durations_K = np.zeros((len(original_g_BKs), len(original_g_Ks)))
    for i, g_BK in enumerate(tqdm(g_BKs, desc="Varying g_BK")):
        for j, g_K in enumerate(tqdm(g_Ks, desc="Varying g_K")):
            time, voltage, info = model.run(g_BK=g_BK,
                                            g_K=g_K,
                                            discard=discard,
                                            noise_amplitude=noise_amplitude,
                                            simulation_time=simulation_time)


            tmp_duration = np.mean(duration(time, voltage))

            if np.isnan(tmp_duration):
                tmp_duration = -1

            event_durations_K[i, j] = tmp_duration

    np.save(os.path.join(data_folder, "original_g_BKs"), original_g_BKs)
    np.save(os.path.join(data_folder, "original_g_SKs"), original_g_SKs)
    np.save(os.path.join(data_folder, "original_g_Ks"), original_g_Ks)
    np.save(os.path.join(data_folder, "event_durations_SK"), event_durations_SK)
    np.save(os.path.join(data_folder, "event_durations_K"), event_durations_K)

    return original_g_BKs, original_g_SKs, original_g_Ks, event_durations_SK, event_durations_K




def plot_exploration(original_g_BKs,
                     original_g_SKs,
                     original_g_Ks,
                     event_durations_SK,
                     event_durations_K):
    """
    Plot the results from the parameter exploration. Figure saved as
    durations.eps.

    Parameters
    ----------
    original_g_BKs : array
        The g_BK values used in the exploration, not scaled.
    original_g_SKs : array
        The g_SK values used in the exploration, not scaled.
    original_g_Ks : array
        The g_K values used in the exploration, not scaled.
    event_durations_SK : array
        The event durations when g_SK and g_BK is changed.
    event_durations_K : array
        The event durations when g_K and g_BK is changed.
    """
    set_style("seaborn-white")

    plt.rcParams.update({"font.family": "serif", "axes.titlepad": 10})

    label_x = 0.05
    label_y = 0.95

    increased_titlesize = titlesize + 4
    increased_labelsize = labelsize + 2
    increased_fontsize = fontsize + 1

    fig = plt.figure(figsize=(figure_width, 1.3*figure_width))
    ax1 = fig.add_subplot(2, 1, 1, projection='3d')
    ax2 = fig.add_subplot(2, 1, 2, projection='3d')

    colorblind = [(0.00392156862745098, 0.45098039215686275, 0.6980392156862745),
                  (0.8705882352941177, 0.5607843137254902, 0.0196078431372549),
                  (0.00784313725490196, 0.6196078431372549, 0.45098039215686275)]

    cmap = colors.ListedColormap(colorblind)

   # Plotting K
    x, y = np.meshgrid(original_g_Ks, original_g_BKs)

    bounds = [-1, 0, 60, np.max(event_durations_K) + 1]
    norm = colors.BoundaryNorm(bounds, 3)

    wire = plot_wireframe(ax1, x, y, event_durations_K,
                          cmap=cmap,
                          norm=norm,
                          rcount=nr_points_exploration,
                          ccount=nr_points_exploration)


    ax1.set_xlabel(r"$G_\mathrm{K}$ (nS)",
                   fontsize=increased_labelsize,
                   fontweight=fontweight,
                   labelpad=10)
    ax1.set_ylabel(r"$G_\mathrm{BK}$ (nS)",
                   fontsize=increased_labelsize,
                   fontweight=fontweight,
                   labelpad=10)
    ax1.set_zlabel("Event duration (ms)",
                   fontsize=increased_labelsize,
                   fontweight=fontweight,
                  labelpad=10)
    ax1.set_title("$G_\mathrm{K}$",
                  fontsize=increased_titlesize,
                  fontweight=fontweight)
    ax1.text2D(label_x,
               label_y,
               "A",
               transform=ax1.transAxes,
               fontsize=increased_titlesize,
               fontweight=plot_label_weight)

    ax1.set_xlim3d(min(original_g_Ks), max(original_g_Ks))
    ax1.set_ylim3d(min(original_g_BKs), max(original_g_BKs))

    ax1.tick_params(axis="both",
                    labelcolor="black",
                    labelsize=increased_fontsize)

    # Plotting SK
    x, y = np.meshgrid(original_g_SKs, original_g_BKs)

    bounds = [-1, 0, 60, np.max(event_durations_SK) + 1]
    norm = colors.BoundaryNorm(bounds, 3)

    wire = plot_wireframe(ax2, x, y, event_durations_SK,
                          cmap=cmap,
                          norm=norm,
                          rcount=nr_points_exploration,
                          ccount=nr_points_exploration)

    ax2.set_xlabel(r"$G_\mathrm{SK}$ (nS)",
                   fontsize=increased_labelsize,
                   fontweight=fontweight,
                   labelpad=10)
    ax2.set_ylabel(r"$G_\mathrm{BK}$ (nS)",
                   fontsize=increased_labelsize,
                   fontweight=fontweight,
                   labelpad=10)
    ax2.set_zlabel("Event duration (ms)",
                   fontsize=increased_labelsize,
                   fontweight=fontweight,
                   labelpad=10)
    ax2.set_title(r"$G_\mathrm{SK}$", fontsize=increased_titlesize)
    ax2.text2D(label_x,
               label_y,
               "B",
               transform=ax2.transAxes,
               fontsize=increased_titlesize,
               fontweight=plot_label_weight)

    ax2.set_xlim3d(min(original_g_SKs), max(original_g_SKs))
    ax2.set_ylim3d(min(original_g_BKs), max(original_g_BKs))

    ax2.tick_params(axis="both", labelcolor="black", labelsize=increased_fontsize)

    plt.tight_layout(pad=2.5)
    plt.savefig(os.path.join(figure_folder, "durations" + figure_format))


def load_exploration():
    """
    Load the results from a parameter exploration.

    Returns
    -------
    original_g_BKs : array
        The g_BK values used in the exploration, not scaled.
    original_g_SKs : array
        The g_SK values used in the exploration, not scaled.
    original_g_Ks : array
        The g_K values used in the exploration, not scaled.
    event_durations_SK : array
        The event durations when g_SK and g_BK is changed.
    event_durations_K : array
        The event durations when g_K and g_BK is changed.
    """
    original_g_BKs = np.load(os.path.join(data_folder, "original_g_BKs.npy"))
    original_g_SKs = np.load(os.path.join(data_folder, "original_g_SKs.npy"))
    original_g_Ks = np.load(os.path.join(data_folder, "original_g_Ks.npy"))
    event_durations_SK = np.load(os.path.join(data_folder, "event_durations_SK.npy"))
    event_durations_K = np.load(os.path.join(data_folder, "event_durations_K.npy"))

    return original_g_BKs, original_g_SKs, original_g_Ks, event_durations_SK, event_durations_K


def perform_uncertainty_analysis():
    """
    Perform an uncertainty quantification and sensitivity analysis of the model.

    Returns
    -------
    data : uncertainpy.Data object
        The results from the uncertainty quantification.
    """
    parameters = {"g_K": scale_conductance(G_K),
                  "g_Ca": scale_conductance(G_Ca),
                  "g_SK": scale_conductance(G_SK),
                  "g_l": scale_conductance(G_l),
                  "g_BK": scale_conductance(0.67)} # Temporary value

    parameters = un.Parameters(parameters)

    # Set all conductances to have a uniform distribution
    # within a +/- 50% interval around their original value
    parameters.set_all_distributions(un.uniform(1))

    parameters["g_BK"].distribution = cp.Uniform(scale_conductance(0),
                                                 scale_conductance(1))


    # Initialize features with the correct algorithm
    features = un.SpikingFeatures(new_features=burstiness_factor,
                                  strict=False,
                                  logger_level="error",
                                  features_to_run=features_to_run,
                                  threshold=0.55,
                                  end_threshold=-0.1,
                                  normalize=True,
                                  trim=False,
                                  min_amplitude=min_spike_amplitude)

    # Initialize the model and define default options
    model = un.NeuronModel(file="tabak.py",
                           name="tabak",
                           discard=discard,
                           noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time,
                           ignore=True)

    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model,
                                      parameters=parameters,
                                      features=features)

    # We set the seed to be able to reproduce the result
    data = UQ.quantify(seed=10,
                       plot=None,
                       data_folder=data_folder,
                       polynomial_order=polynomial_order)

    return data



def parameter_exploration():
    """
    Perform the parameter exploration and plot the results.
    """
    results = perform_exploration()
    # This can be used to load previously generated results
    # results = load_exploration()

    plot_exploration(*results)



def uncertainty_analysis():
    """
    Perform the uncertainty quantification and plot the results.
    """
    data = perform_uncertainty_analysis()
    # This can be used to load previously generated results
    # data = un.Data(os.path.join(data_folder, "tabak.h5"))

    plot_sobol(data, os.path.join(figure_folder, "sensitivity"))


if __name__ == "__main__":
    uncertainty_analysis()
    parameter_exploration()