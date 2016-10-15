"""
   Copyright 2016 Aaron R. Shifman

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

from scalebars import *
from ode_functions import *

from scipy.signal import argrelextrema
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

"""
Static model parameters (constant across all simulations)
"""
g_na_bar = 0.7
g_nap_bar = 0.05
g_k_bar = 1.3
g_p_bar = 0.05
g_leak = 0.005
g_a_bar = 1.0

e_na = 60
e_k = -80
e_leak = -50
e_ca = 40
e_h = -38.8

k1 = -0.0005
k2 = 0.04

cm = 0.040


def figure1():
    """
    Generate figure 1

    This is a 100ms simulation with "default model parameters"
    :return: None
    """

    plot_settings = {'y_limits': [-80, -50],
                     'x_limits': None,
                     'y_ticks': [-80, -70, -60, -50],
                     'locator_size': 5,
                     'y_label': 'Voltage (mV)',
                     'x_ticks': [],
                     'scale_size': 20,
                     'x_label': "",
                     'scale_loc': 4,
                     'figure_name': 'figure_1',
                     'legend_size': 8,
                     'legend': None,
                     'y_on': True}

    """
    Model parameters
    """
    t1 = 100
    dt = 1e-2
    g_t_bar = 0.1
    g_sk_bar = 0.3
    g_n_bar = 0.05
    t_start = 10
    duration = 1
    i_bias_on = 1
    g_h_bar = 0.005

    t, y = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar])  # Integrate solution

    plt.figure(figsize=(5, 2))  # Create figure
    plt.plot(t, y[:, 0], 'k-')  # Plot solution

    """
    Annotate plot with figures
    """
    plt.gca().annotate('fAHP', xy=(13.5, -65), xytext=(17, -60),
                       arrowprops=dict(facecolor='black', shrink=0, headlength=10, headwidth=5, width=1), )
    plt.gca().annotate('ADP', xy=(15.5, -66), xytext=(25, -65),
                       arrowprops=dict(facecolor='black', shrink=0, headlength=10, headwidth=5, width=1), )
    plt.gca().annotate('mAHP', xy=(38, -77), xytext=(43, -72),
                       arrowprops=dict(facecolor='black', shrink=0, headlength=10, headwidth=5, width=1), )
    alter_figure(plot_settings)  # Alter figure for publication


def figure2():
    """
    Generate figure 2

    This is a 100ms simulation with "default model parameters"
    :return: None
    """

    plot_settings = {'y_limits': [-25, 0],
                     'x_limits': None,
                     'y_ticks': [-25, -20, -15, -10, -5, 0],
                     'locator_size': 2.5,
                     'y_label': 'Current (nA)',
                     'x_ticks': [],
                     'scale_size': 0,
                     'x_label': "",
                     'scale_loc': 4,
                     'figure_name': 'figure_2',
                     'legend': ['I-Na', 'I-NaP'],
                     'legend_size': 8,
                     'y_on': True}

    t1 = 100
    dt = 1e-2
    g_t_bar = 0.1
    g_sk_bar = 0.3
    g_n_bar = 0.05
    t_start = 10
    duration = 1
    i_bias_on = 1
    g_h_bar = 0.005

    t, y = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar])
    t_short = np.where((t >= 8) & (t <= 18))[0]  # shorter time bounds for plots A and C

    v, m, h, m_nap, h_na_p, n, m_t, h_t, m_p, m_n, h_n, z_sk, m_a, h_a, m_h, ca = y[:, ].T  # Extract all variables

    """
    Explicitly calculate all currents
    """
    i_na = g_na_bar * (m ** 3) * h * (v - e_na)
    i_na_p = g_nap_bar * m_nap * h_na_p * (v - e_na)
    i_k = g_k_bar * (n ** 4) * (v - e_k)
    i_leak = g_leak * (v - e_leak)
    i_t = g_t_bar * m_t * h_t * (v - e_ca)
    i_n = g_n_bar * m_n * h_n * (v - e_ca)
    i_p = g_p_bar * m_p * (v - e_ca)
    i_sk = g_sk_bar * (z_sk ** 2) * (v - e_k)
    i_a = g_a_bar * m_a * h_a * (v - e_k)
    # i_h = g_h_bar * m_h * (v - e_h)

    plt.figure(figsize=(5, 3), dpi=96)  # Create figure

    plt.subplot(2, 2, 1)  # Generate subplot 1 (top left)

    plt.plot(t[t_short], i_na[t_short], 'k-')
    plt.plot(t[t_short], i_na_p[t_short], c='k', linestyle='dotted')

    alter_figure(plot_settings)  # Alter figure for publication

    plt.subplot(2, 2, 2)  # Generate subplot 2 (top right)

    plt.plot(t, i_t + i_n + i_p, 'k-')
    plt.plot(t, i_t, c='k', linestyle='dotted')
    plt.plot(t, i_n, 'k-.')
    plt.plot(t, i_p, 'k--')

    plot_settings['y_limits'] = [-2.5, 0]
    plot_settings['y_ticks'] = [-2.5, -2, -1.5, -1, -0.5, 0]
    plot_settings['locator_size'] = 0.25
    plot_settings['y_label'] = ""
    plot_settings['legend'] = ['I-Ca', 'I-T', 'I-P', 'I-N']
    alter_figure(plot_settings)  # Alter figure for publication

    plt.subplot(2, 2, 3)  # Generate subplot 3 (bottom left)

    plt.plot(t[t_short], i_k[t_short], 'k-')
    plt.plot(t[t_short], i_a[t_short], c='k', linestyle='dotted')
    plt.plot(t[t_short], i_leak[t_short], 'k-.')

    plot_settings['y_limits'] = [0, 25]
    plot_settings['y_ticks'] = [0, 5, 10, 15, 20, 25]
    plot_settings['locator_size'] = 2.5
    plot_settings['y_label'] = "Current (nA)"
    plot_settings['legend'] = ['I-K', 'I-A', 'I-leak']
    plot_settings['scale_size'] = 2
    plot_settings['scale_loc'] = 2
    alter_figure(plot_settings)  # Alter figure for publication

    plt.subplot(2, 2, 4)  # Generate subplot 4 (bottom left)

    plt.plot(t, i_sk, 'k-')

    plot_settings['y_limits'] = [0, 1]
    plot_settings['y_ticks'] = [0, 0.2, 0.4, 0.6, 0.8, 1]
    plot_settings['locator_size'] = 0.2
    plot_settings['y_label'] = ""
    plot_settings['legend'] = ['I-SK']
    plot_settings['scale_size'] = 20
    plot_settings['scale_loc'] = 2
    alter_figure(plot_settings)  # Alter figure for publication


def figure3():
    """
    Generate figure 3

    This is a set of 2 simulations
    The first simulation is 7ms in which the role of T type calcium current is explored
    The second simulation is 1020ms in which the neuron is held in a hyperpolarized
    condition prior to the current pulse
    In order to hold and hyperpolarize the neuron the simulation is split into 2 steps.
    in the first step (1000 ms) the simulation is a neuron being hyperpolarized by a negative current
    In the second step, the initial conditions are the end conditions of the first step,
    and a current pulse is applied
    :return: None
    """

    plot_settings = {'y_limits': [-75, -50],
                     'x_limits': None,
                     'y_ticks': [-75, -70, -65, -60, -55, -50],
                     'locator_size': 2.5,
                     'y_label': 'Voltage (mV)',
                     'x_ticks': [],
                     'scale_size': 2,
                     'x_label': "",
                     'scale_loc': 3,
                     'figure_name': 'figure_3',
                     'legend': ['0.1 $\mu$S', '0.125 $\mu$S', '0.15 $\mu$S'],
                     'legend_size': 8,
                     'y_on': True}
    line_styles = ['-', 'dotted', '-.']

    t1 = 7
    dt = 1e-2
    g_sk_bar = 0.3
    g_n_bar = 0.05
    t_start = 0.5
    duration = 1
    i_bias_on = 1
    g_h_bar = 0.005

    plt.figure(figsize=(5, 3))
    plt.subplot(1, 2, 1)  # Generate subplot 1 (left)

    for ix, g_t_bar in enumerate([0.1, 0.125, 0.15]):
        t, y = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar])
        plt.plot(t, y[:, 0], c='k', linestyle=line_styles[ix])

    alter_figure(plot_settings)  # Alter figure for publication

    plt.subplot(1, 2, 2)  # Generate subplot 2 (right)

    dt = 1e-2
    g_t_bar = 0.1
    g_sk_bar = 0.3
    g_n_bar = 0.05
    t_start = 0.5
    g_h_bar = 0.005

    for ix, i_bias_on in enumerate([0, -0.1, -0.2]):
        """
        First step (hyperpolarizing)
        """
        t1 = 1000
        duration = t1

        t, y_hold = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar])
        y_0 = y_hold[-1, :]  # Create new initial conditions

        """
        Second step (current pulse)
        """
        t1 = 20
        i_bias_on = 1
        duration = 1
        t_start = 0

        t, y = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar], y_hold=y_0)
        v = y[:, 0]

        """
        Append the end of the first step onto the second step (for plotting purposes)
        """
        len_pre = 100
        t = np.concatenate((np.linspace(-len_pre * dt, -dt, len_pre), t))
        v_bar = np.concatenate((y_hold[-len_pre:, 0], v))
        v_bar = v_bar - v_bar[0] + -7.16300325e+01  # Align solution to initial condition of the "standard simulation"

        plt.plot(t, v_bar, c='k', linestyle=line_styles[ix])

    plot_settings['y_ticks'] = []
    plot_settings['y_label'] = ""
    plot_settings['x_limits'] = [-1, 20]
    plot_settings['legend'] = ['-72 mV', '-75 mV', '-78 mV']
    plot_settings['scale_size'] = 5
    alter_figure(plot_settings)
    plt.savefig('figures/figure_3.pdf', dpi=1200)


def figure4():
    """
    Generate figure 4

    This is a 100ms simulation in which the role of the SK current is shown with simulated apamin application
    (gSK = 0)
    Furthermore a 1200ms simulation is run with standard SK current
    :return: None
    """

    plot_settings = {'y_limits': [-80, -50],
                     'x_limits': None,
                     'y_ticks': [-80, -70, -60, -50],
                     'locator_size': 5,
                     'y_label': 'Voltage (mV)',
                     'x_ticks': [],
                     'scale_size': 20,
                     'x_label': "",
                     'scale_loc': 4,
                     'figure_name': 'figure_4',
                     'legend': ['control', 'apamin'],
                     'legend_size': 8,
                     'y_on': True}
    line_styles = ['-', 'dotted']

    t1 = 100
    dt = 1e-2
    g_t_bar = 0.1
    g_n_bar = 0.05
    t_start = 10
    duration = 1
    i_bias_on = 1
    g_h_bar = 0.005

    plt.figure(figsize=(5, 3), dpi=96)

    plt.subplot(2, 1, 1)  # Generate figure 1 (top)

    for ix, g_sk_bar in enumerate([0.3, 0]):
        t, y = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar])
        plt.plot(t, y[:, 0], c='k', linestyle=line_styles[ix])

    alter_figure(plot_settings)  # Alter figure for publication

    t1 = 1200
    dt = 1e-2
    g_t_bar = 0.1
    g_sk_bar = 0.03
    g_n_bar = 0.05
    t_start = 50
    duration = t1
    i_bias_on = 0.33
    g_h_bar = 0.005

    plt.subplot(2, 1, 2)
    t, y = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar])
    plt.plot(t, y[:, 0], 'k-')

    plot_settings['y_limits'] = [-100, 30]
    plot_settings['x_limits'] = [0, t1]
    plot_settings['y_ticks'] = [-80, -60, -40, -20, 0, 20]
    plot_settings['locator_size'] = 10
    plot_settings['scale_size'] = 100
    plot_settings['legend'] = None
    alter_figure(plot_settings)  # Alter plot for publication


def figure5():
    """
    Generate figure 5

    This is a 2500ms simulation in which firing cessation is demonstrated following
    a constant bias current of 0.22nA
    :return: None
    """

    plot_settings = {'y_limits': [-100, 30],
                     'x_limits': None,
                     'y_ticks': [-80, -60, -40, -20, 0, 20],
                     'locator_size': 10,
                     'y_label': 'Voltage (mV)',
                     'x_ticks': [],
                     'scale_size': 500,
                     'x_label': "",
                     'scale_loc': 4,
                     'figure_name': 'figure_5',
                     'legend': None,
                     'legend_size': 8,
                     'y_on': True}

    t1 = 2500
    dt = 1e-2
    g_t_bar = 0.1
    g_sk_bar = 0.3
    g_n_bar = 0.05
    t_start = 60
    duration = 2000
    i_bias_on = 0.22
    g_h_bar = 0.005

    t, y = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar])

    plt.figure(figsize=(5, 2))

    plt.plot(t, y[:, 0], 'k-')

    ix_start = np.where(t < 40)[0][-1]  # Find index for beginning of inset
    ix_end = np.where(t < 160)[0][-1]  # Find index for end of inset

    alter_figure(plot_settings)  # Alter figure for publication

    plt.gca().add_patch(patches.Rectangle((40, -75), 120, 16, fill=False))  # Draw rectangle to highlight inset

    """
    Create inset of highlighted region
    """
    plt.axes([.75, .5, .25, .4], axisbg='y')
    v_highlighted = y[ix_start:ix_end, 0]
    plt.plot(t[ix_start:ix_end], v_highlighted, 'k')
    plt.ylim([-75, -55])
    plt.box('off')
    plt.xticks([])
    plt.yticks([])
    plt.xlim([t[ix_start], t[ix_end]])
    add_scalebar(plt.gca(), matchx=False, matchy=False, sizex=25, sizey=5, hidey=False, labelx='25', labely='5', loc=1)

    plt.savefig('figures/figure_5.pdf', dpi=1200)
    plt.close()


def figure6():
    """
    Generate figure 6

    This is a 240ms simulation to determine spike frequency stabilization with attenuated T current
    :return: None
    """

    plot_settings = {'y_limits': [-100, 30],
                     'x_limits': None,
                     'y_ticks': [-80, -60, -40, -20, 0, 20],
                     'locator_size': 10,
                     'y_label': 'Voltage (mV)',
                     'x_ticks': [],
                     'scale_size': 50,
                     'x_label': "",
                     'scale_loc': 3,
                     'figure_name': 'figure_6',
                     'legend': None,
                     'legend_size': 8,
                     'y_on': True}

    marker = ['o', 's', '^']
    line_styles = ['-', 'dotted', '--']

    t1 = 240
    dt = 1e-2
    g_t_bar = 0.1 / 10
    g_sk_bar = 0.3
    g_n_bar = 0.05
    t_start = 10
    duration = 250
    i_bias_on = 2
    g_h_bar = 0.005

    plt.figure(figsize=(5, 3), dpi=96)

    plt.subplot(2, 1, 1)  # Generate subplot 1 (top)

    t, y = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar])
    plt.plot(t, y[:, 0], 'k-')
    alter_figure(plot_settings)

    t1 = 240
    dt = 1e-2
    g_t_bar = 0.1 / 10
    g_sk_bar = 0.3
    g_n_bar = 0.05
    t_start = 10
    duration = 250
    g_h_bar = 0.005

    plt.subplot(2, 1, 2)  # Generate subplot 2 (bottom)

    for ix, i_bias_on in enumerate([2, 1.5, 1]):
        t, y = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar])
        t_spike, f = spike_times(t, y[:, 0])

        plt.plot(t_spike[0:-1], f, c='k', linestyle=line_styles[ix], marker=marker[ix], fillstyle='none')

    plot_settings['y_limits'] = [0, 200]
    plot_settings['y_ticks'] = [0, 50, 100, 150, 200]
    plot_settings['locator_size'] = 25
    plot_settings['y_label'] = 'Frequency (Hz)'
    plot_settings['legend'] = ['2.0 nA', '1.5 nA', '1.0 nA']
    plot_settings['scale_size'] = 0
    alter_figure(plot_settings)  # Alter figure for publication


def figure7():
    """
    Generate figure 7

    This is a 250ms simulation to determine spike frequency stabilization with full T current
    :return: None
    """
    plot_settings = {'y_limits': [-100, 30],
                     'x_limits': None,
                     'y_ticks': [-80, -60, -40, -20, 0, 20],
                     'locator_size': 10,
                     'y_label': 'Voltage (mV)',
                     'x_ticks': [],
                     'scale_size': 50,
                     'x_label': "",
                     'scale_loc': 3,
                     'figure_name': 'figure_7',
                     'legend': None,
                     'legend_size': 8,
                     'y_on': True}

    marker = ['o', 's', '^']
    line_styles = ['-', 'dotted', '--']

    plt.figure(figsize=(5, 3), dpi=96)
    plt.subplot(2, 1, 1)  # Generate subplot 1 (top)

    t1 = 250
    dt = 1e-2
    g_t_bar = 0.1
    g_sk_bar = 0.3
    g_n_bar = 0.05
    t_start = 10
    duration = 260
    i_bias_on = 2
    g_h_bar = 0.005

    t, y = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar])
    plt.plot(t, y[:, 0], 'k-')
    alter_figure(plot_settings)

    plt.subplot(2, 1, 2)  # Generate subplot 2 (bottom)

    t1 = 250
    dt = 1e-2
    g_t_bar = 0.1
    g_sk_bar = 0.3
    g_n_bar = 0.05
    t_start = 10
    duration = 260
    g_h_bar = 0.005

    for ix, i_bias_on in enumerate([2, 1.5, 1]):
        t, y = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar])
        t_spike, f = spike_times(t, y[:, 0])
        plt.plot(t_spike[0:-1], f, c='k', linestyle=line_styles[ix], marker=marker[ix], fillstyle='none')

    plot_settings['y_limits'] = [20, 40]
    plot_settings['y_ticks'] = [20, 25, 30, 35, 40]
    plot_settings['locator_size'] = 2.5
    plot_settings['y_label'] = 'Frequency (Hz)'
    plot_settings['legend'] = ['2.0 nA', '1.5 nA', '1.0 nA']
    plot_settings['scale_size'] = 0
    plot_settings['legend_location'] = 4
    alter_figure(plot_settings)


def figure8():
    """
    Generate figure 8

    200ms simulation in order to show the effects of T and N calcium currents on spike frequency adaptation
    :return: None
    """

    plot_settings = {'y_limits': [15, 60],
                     'x_limits': None,
                     'y_ticks': [20, 30, 40, 50, 60],
                     'locator_size': 5,
                     'y_label': 'ISI (ms)',
                     'x_ticks': [],
                     'scale_size': 0,
                     'x_label': "",
                     'scale_loc': 4,
                     'figure_name': 'figure_8',
                     'legend': ['First ISI', 'Second ISI'],
                     'legend_size': 8,
                     'y_on': True,
                     'legend_location': 3}

    t1 = 200
    dt = 1e-2
    g_sk_bar = 0.3
    g_n_bar = 0.05
    t_start = 15
    duration = 260
    i_bias_on = 1
    g_h_bar = 0.005

    g_t_bars = np.linspace(0.02, 0.2, 10)
    isi = np.zeros((len(g_t_bars), 2))

    for ix, g_t_bar in enumerate(g_t_bars):
        t, y = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar])
        t_spike, f = spike_times(t, y[:, 0])

        isi[ix, 0] = t_spike[1] - t_spike[0]
        isi[ix, 1] = t_spike[2] - t_spike[1]

    plt.subplot(2, 2, 1)  # Generate subplot 1 (top left)

    plt.plot(g_t_bars, isi[:, 0], c='k', marker='o', fillstyle='none', linestyle='-')
    plt.plot(g_t_bars, isi[:, 1], c='k', marker='s', fillstyle='none', linestyle='dotted')

    """
    Annotate plot
    """
    plt.gca().arrow(g_t_bars[3], 29, 0, 22, head_width=0, head_length=0, fc='k', ec='k')
    plt.gca().arrow(g_t_bars[3], 51, -0.01, 0, head_width=2, head_length=0.01, fc='k', ec='k')
    plt.gca().arrow(g_t_bars[3], 29, 0.01, 0, head_width=2, head_length=0.01, fc='k', ec='k')
    plt.gca().annotate("Acceleration", (0.1, 29), fontsize=8)
    plt.gca().annotate("Adaptation", (0.01, 51), fontsize=8)
    alter_figure(plot_settings)

    plt.subplot(2, 2, 2)  # Generate subplot 2 (top right)

    t1 = 200
    dt = 1e-2
    g_sk_bar = 0.3
    g_t_bar = 0.02
    t_start = 15
    duration = 260
    i_bias_on = 1
    g_h_bar = 0.005

    g_n_bars = np.linspace(0.02, 0.2, 10)
    isi = np.zeros((len(g_t_bars), 2))
    for ix, g_n_bar in enumerate(g_n_bars):
        t, y = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar])
        t_spike, f = spike_times(t, y[:, 0])

        isi[ix, 0] = t_spike[1] - t_spike[0]
        isi[ix, 1] = t_spike[2] - t_spike[1]

    plt.plot(g_t_bars, isi[:, 0], c='k', marker='o', fillstyle='none', linestyle='-')
    plt.plot(g_t_bars, isi[:, 1], c='k', marker='s', fillstyle='none', linestyle='dotted')

    """
    Annotate plot
    """
    plt.gca().arrow(g_n_bars[3], 25, 0, 20, head_width=0, head_length=0, fc='k', ec='k')
    plt.gca().arrow(g_n_bars[3], 45, -0.01, 0, head_width=2, head_length=0.01, fc='k', ec='k')
    plt.gca().arrow(g_n_bars[3], 25, 0.01, 0, head_width=2, head_length=0.01, fc='k', ec='k')
    plt.gca().annotate("Acceleration", (0.1, 25), fontsize=8)
    plt.gca().annotate("Adaptation", (0.015, 45), fontsize=8)
    plot_settings['y_ticks'] = []
    plot_settings['y_label'] = ""
    plot_settings['y_on'] = False
    plot_settings['legend_location'] = 4
    alter_figure(plot_settings)

    plt.subplot(2, 2, 3)  # Generate subplot 3 (bottom left)

    t1 = 200
    dt = 1e-2
    g_sk_bar = 0.3
    g_n_bar = 0.05
    t_start = 15
    duration = 260
    i_bias_on = 1
    g_h_bar = 0.005

    g_t_bars = np.linspace(0.02, 0.16, 8)
    isi = np.zeros((len(g_t_bars), 2))
    for ix, g_t_bar in enumerate(g_t_bars):
        t, y = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar], ca_type=1)
        t_spike, f = spike_times(t, y[:, 0])

        isi[ix, 0] = t_spike[1] - t_spike[0]
        isi[ix, 1] = t_spike[2] - t_spike[1]

    plt.plot(g_t_bars, isi[:, 0], c='k', marker='o', fillstyle='none', linestyle='-')
    plt.plot(g_t_bars, isi[:, 1], c='k', marker='s', fillstyle='none', linestyle='dotted')

    plt.gca().arrow(g_t_bars[2], 25, -0.02, 0, head_width=2, head_length=0.01, fc='k', ec='k')
    plt.gca().arrow(g_t_bars[4], 25, 0.02, 0, head_width=2, head_length=0.01, fc='k', ec='k')
    plt.gca().annotate("Adaptation", (0.06, 25), fontsize=8)

    plot_settings['y_limits'] = [0, 45]
    plot_settings['y_ticks'] = [0, 10, 20, 30, 40]
    plot_settings['locator_size'] = 5
    plot_settings['y_label'] = 'ISI (ms)'
    plot_settings['y_on'] = True
    plot_settings['legend_location'] = 3
    alter_figure(plot_settings)

    plt.subplot(2, 2, 4)

    t1 = 200
    dt = 1e-2
    g_sk_bar = 0.3
    g_t_bar = 0.02
    t_start = 15
    duration = 260
    i_bias_on = 1
    g_h_bar = 0.005

    g_n_bars = np.linspace(0.02, 0.16, 8)
    isi = np.zeros((len(g_t_bars), 2))
    for ix, g_n_bar in enumerate(g_n_bars):
        t, y = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar], ca_type=2)
        t_spike, f = spike_times(t, y[:, 0])

        isi[ix, 0] = t_spike[1] - t_spike[0]
        isi[ix, 1] = t_spike[2] - t_spike[1]

    plt.plot(g_t_bars, isi[:, 0], c='k', marker='o', fillstyle='none', linestyle='-')
    plt.plot(g_t_bars, isi[:, 1], c='k', marker='s', fillstyle='none', linestyle='dotted')

    plt.gca().arrow(g_n_bars[2], 20, -0.02, 0, head_width=2, head_length=0.01, fc='k', ec='k')
    plt.gca().arrow(g_n_bars[4], 20, 0.02, 0, head_width=2, head_length=0.01, fc='k', ec='k')
    plt.gca().annotate("Adaptation", (0.06, 20), fontsize=8)

    plot_settings['y_ticks'] = []
    plot_settings['y_label'] = ''
    plot_settings['y_on'] = False
    alter_figure(plot_settings)


def figure9():
    """
    Generate figure 9

    2000ms simulation in order to show the effects of H current on mAHP length
    :return: None
    """
    plot_settings = {'y_limits': [70, 115],
                     'x_limits': [0.005, 0.05],
                     'y_ticks': [70, 80, 90, 100, 110],
                     'locator_size': 5,
                     'y_label': 'mAHP Duration (ms)',
                     'x_ticks': [],
                     'scale_size': 0,
                     'x_label': "",
                     'scale_loc': 4,
                     'figure_name': 'figure_9',
                     'legend': None,
                     'legend_size': 8,
                     'y_on': True}

    t1 = 5000  # long to allow V to stabilize to new rest location
    dt = 1e-3
    g_sk_bar = 0.3
    g_t_bar = 0.1
    g_n_bar = 0.05
    t_start = 3000  # long to allow V to stabilize to new rest location
    duration = 1
    i_bias_on = 1

    g_h_bars = np.linspace(0.005, 0.05, 10)
    length_after = np.zeros((len(g_h_bars),))

    for ix, gHBar in enumerate(g_h_bars):
        t, y = solver(t1, dt, [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, gHBar])
        v = y[:, 0]

        pk = np.where(v == np.max(v))[0][0]
        v_clipped = v[pk:]
        v_rest = v[np.where(t <= t_start)[0][-1]]  # v_rest is v immediately before the stimulus turns on
        v_clipped -= v_rest
        crossings = np.where((v_clipped[:-1] * v_clipped[1:]) < 0)[0]
        ix_after = np.where(t[crossings] < 200)[0][-1]

        length_after[ix] = t[crossings[ix_after]-crossings[ix_after-1]]

    plt.figure(figsize=(5, 3), dpi=96)
    """
    x is digitized from figure 9 in the original manuscript
    """
    x = [108.44444444444443, 97.03703703703704, 89.7037037037037, 85.55555555555556, 82.22222222222223,
         80.2962962962963, 78.22222222222223, 77.18518518518519, 76.81481481481481, 74.07407407407408]
    plt.plot(g_h_bars, length_after, c='k', marker='o', fillstyle='none')
    #plt.plot(g_h_bars, x)

    ellipse = patches.Ellipse(xy=(0.017, 105), width=0.01, height=4, angle=0)
    plt.gca().add_artist(ellipse)
    ellipse.set_facecolor((1, 1, 1))
    plt.gca().annotate("Neonatal", (0.017, 105), fontsize=8, ha="center", va="center")

    ellipse = patches.Ellipse(xy=(0.04, 72), width=0.005, height=4, angle=0)
    plt.gca().add_artist(ellipse)
    ellipse.set_facecolor((1, 1, 1))
    plt.gca().annotate("Adult", (0.04, 72), fontsize=8, ha="center", va="center")

    alter_figure(plot_settings)
    plt.xticks(g_h_bars[0::2])
    plt.xlabel('$\\bar{g_H}$ ($\mu S$)')
    plt.tight_layout()
    plt.gca().spines['bottom'].set_visible(True)
    plt.gca().xaxis.set_ticks_position('bottom')

    plt.savefig('figures/figure_9.pdf', dpi=1200)


def alter_figure(settings):
    plt.ylim(settings['y_limits'])
    plt.yticks(settings['y_ticks'])
    plt.xticks(settings['x_ticks'])
    plt.gca().tick_params(axis='y', direction='out')
    plt.gca().tick_params(axis='x', direction='out')
    minor_locator = MultipleLocator(settings['locator_size'])
    plt.gca().yaxis.set_minor_locator(minor_locator)
    plt.ylabel(settings['y_label'])
    plt.xlabel(settings['x_label'])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().xaxis.set_ticks_position('none')

    if 'y_on' in settings:
        plt.gca().yaxis.set_ticks_position('left')
    else:
        plt.gca().yaxis.set_ticks_position('none')
        plt.gca().spines['left'].set_visible(False)

    plt.xlim(settings['x_limits'])
    if settings['scale_size'] != 0:
        add_scalebar(plt.gca(), matchx=False, matchy=False, sizex=settings['scale_size'], hidey=False,
                     labelx=str(settings['scale_size']) + ' ms', loc=settings['scale_loc'])
    if settings['legend'] is not None:
        if 'legend_location' not in settings:
            plt.legend(settings['legend'], frameon=False, fontsize=settings['legend_size'])
        else:
            plt.legend(settings['legend'], frameon=False, fontsize=settings['legend_size'],
                       loc=settings['legend_location'])
    plt.tight_layout()
    plt.savefig('figures/' + settings['figure_name'] + '.pdf', dpi=1200)


def spike_times(t, v):
    spike = argrelextrema(v, np.greater)[0]  # Get spike indices
    spike = spike[v[spike] > 0]  # Threshold all spikes >0 (filter for true action potential peak)
    t_spike = t[spike]  # Get spike times

    return t_spike, 1000 / np.diff(t_spike)
