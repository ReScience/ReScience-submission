#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This script implements the Mihalas-Niebur Neuron model based on the paper
# S. Mihalas and E. Niebur, "A Generalized Linear Integrate-and-Fire Neural
# Model Produces Diverse Spiking Behaviors", Neural Computation 21, 2009.
# Copyright (C) 2017  Georgios Is. Detorakis
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
import numpy as np
from copy import copy
from matplotlib import gridspec
import matplotlib.pylab as plt
from neuron_model import mnn_model


# Define simulation and model parameters dictionary
params = {'k1': 0.2,             # I₁ decay term
          'k2': 0.02,            # I₂ decay term
          'b': 0.01,             # Θ decay term
          'R1': 0.0,             # I₁ update constant
          'R2': 1.0,             # I₂ update constant
          'El': -70.0,           # Reverse potential
          'Vr': -70.0,           # Resting potential
          'Thetar': -60.0,       # Resting threshold
          'a': 0.000,            # Specifier of spike response
          'A1': 0.0,             # I₁ additive update constant
          'A2': 0.0,             # I₂ additive update constant
          'G': 0.05,             # Membrane potential decay term
          'C': 1.0,              # Membrane capacitance
          'ThetaInf': -50.0}       # Reverse threshold


def plot_preparation(V, Theta, theta_, thr=-3.5):
    """
        This function preprocesses the data for pretty plots. It removes not
        useful data from the input arrays.
        It removes all the values lower than a threshold for the membrane
        potential variable (V), and all the zeros from threshold potential
        at the spike events.

        Args:
            V  (ndarray)        : Membrane potential
            Theta (ndarray)     : Threshold potential
            theta_ (ndarray)    : Threshold potential before spike event
            thr (float)         : Derivative threshold

        Return:
            v_plot  (ndarray)   : Prettify membrane potential
            theta   (ndarray)   : Prettify threshold potential
            theta_  (ndarray)   : Prettify threshold potential
    """
    # Copy the data
    v_plot, theta = V.copy(), Theta.copy()

    # Find where the derivative of V is less than a threshold
    pos = np.where((np.diff(v_plot) < thr))[0]

    # Set all the elements at pos to NaN (matplotlib ignores them)
    v_plot[pos] = np.nan

    # Remove all the reset values of the threshold variable
    theta[(theta_ != 0)] = np.nan
    theta_[theta_ == 0] = np.nan
    return v_plot, theta, theta_


def plot_figure2(V_, V, Theta, Spk, Theta_, pms=params, savefig=False):
    """
        This function plots Figure 2 in the paper (phasic spiking
        phase space).

        Args:
            V_  (ndarray)       : Filtered membrane potential
            V   (ndarray)       : Membrane potential
            Theta (ndarray)     : Threshold potential
            Spk (ndarray)       : Spike events
            Theta_  (ndarray)   : Threshold potentials before spike events
            savefig (boolean)   : Enables figure saving in pdf format

        Return:

    """
    fig2 = plt.figure(figsize=(10, 10))
    gs1 = gridspec.GridSpec(1, 1)
    # ax = fig2.add_subplot(111)
    ax = plt.subplot(gs1[0])

    ax.plot(V_, Theta)
    ax.plot(V[Spk-1], Theta[Spk-1], 'k.', ms=10)

    p = np.polyfit(V[Spk-1], Theta[Spk-1], deg=1)
    x_the = np.linspace(-75, -25, 100)
    y_the = np.poly1d(p)(x_the)
    ax.plot(x_the, y_the, 'k--', alpha=0.5)

    X, Y = np.meshgrid(np.linspace(-75, -28, 8),
                       np.linspace(-55, -28, 8))
    V = 1.5 - pms['G'] * (X - pms['El'])
    U = pms['a'] * (X - pms['El']) - pms['b'] * (Y - pms['ThetaInf'])
    ax.quiver(X, Y, V, U, width=.003, alpha=0.5, headwidth=7)

    ax.set_xlim([-75, -28])
    ax.set_ylim([-55, -28])
    ax.set_title('Phasic spiking', fontsize=16, weight='bold')
    ax.set_xlabel('Membrane potential (mV)', fontsize=16, weight='bold')
    ax.set_ylabel('Threshold potential (mV)', fontsize=16, weight='bold')
    ax.set_xticklabels(ax.get_xticks().astype('i'), fontsize=14, weight='bold')
    ax.set_yticklabels(ax.get_yticks().astype('i'), fontsize=14, weight='bold')
    if savefig is True:
        fig2.show()
        plt.savefig('Figure02.pdf', axis='tight')
        plt.close()


def plot_figure3(I1, I2, V, V_, Spk, savefig=False):
    """
        This function plots Figure 3 in the paper (phasic bursting phase
        space).

        Args:
            I1  (ndarray)       : First internal current
            I2   (ndarray)      : Second internal current
            V (ndarray)         : Membrane potential
            V_ (ndarray)        : Filtered membrane potential
            Spk (ndarray)       : Spike events
            savefig (boolean)   : Enables figure saving in pdf format

        Return:

    """
    fig3 = plt.figure(figsize=(11, 11))
    from mpl_toolkits.mplot3d import Axes3D
    ax = fig3.add_subplot(111)
    ax = fig3.gca(projection='3d')

    ax.plot(I1, I2, V_)
    ax.plot(I1[Spk-1], I2[Spk-1], V[Spk-1], 'k.', ms=15)

    ax.set_xlim([0, 10])
    ax.set_ylim([-3, 0])
    ax.set_title('Phasic bursting', fontsize=16, weight='bold')
    ax.set_ylabel(r'$I_2 (V/s)$', fontsize=16, weight='bold')
    ax.set_xlabel(r'$I_1 (V/s)$', fontsize=16, weight='bold')
    ax.set_zlabel('Membrane potential (mV)', fontsize=16, weight='bold')
    ax.set_xticklabels(ax.get_xticks().astype('i'), fontsize=14, weight='bold')
    ax.xaxis.labelpad = 20
    ax.set_yticklabels(ax.get_yticks(), fontsize=14, weight='bold')
    ax.yaxis.labelpad = 30
    ax.set_zticklabels(ax.get_zticks().astype('i'), fontsize=14, weight='bold')
    ax.zaxis.labelpad = 30
    if savefig is True:
        fig3.show()
        plt.savefig('Figure03.pdf', axis='tight')
        plt.close()


def run_and_plot(Iext, params, figure, subplots_id,
                 IC=(0.01, 0.001, -70.0, -50.0)):
    """
        This function runs a simulation, prettyfies the results and
        plots them.

        Args:
            Iext    (ndarray)       : External input to the model
            params  (dict)          : Dictionary containing simulation
                                      parameters
            figure  ()              : Matplotlib plot container
            subplots_id (tuple)     : Subplots configuration
            IC      (tuple)         : Initial conditions

        Return:
            axes_inst (matplotlib object)   : Axes instance
    """
    axes_inst = plt.subplot(subplots_id)

    # Run the simulation
    sol, theta_, spk = mnn_model(params, Iext, dt=0.1, IC=IC)

    # Prettify the results
    v, theta, _ = plot_preparation(sol[3], sol[4], theta_)

    # Plot the results
    axes_inst.plot(sol[0], v, 'k', lw=1.5)
    axes_inst.plot(sol[0], sol[4], 'r--', lw=1.5)
    axes_inst.plot(sol[0], theta_, 'k.', lw=1.5, ms=10)
    axes_inst.plot(sol[0], params['Vr'] + Iext, 'k', lw=1.5, alpha=0.5)
    axes_inst.set_xlim([0, int(Iext.shape[0]*0.1)])

    return axes_inst, sol, v, theta_, spk


# --------------------------------------------------------------------------
if __name__ == '__main__':
    """
        This script plots all the figures in the paper (Figures 1, 2, 3).
    """
    axx = []

    # Subplots titles
    titles = ['Tonic spiking', 'Class 1', 'Spike frequency \n adaptation',
              'Phasic spiking', 'Accomodation', 'Threshold variability',
              'Rebound spike', 'Class 2', 'Integrator', 'Input bistability',
              'Hyperpolarization induced \n spiking',
              'Hyperpolarization induced \n bursting', 'Tonic bursting',
              'Phasic bursting', 'Rebound burst',
              'Mixed mode', 'Afterpotentials', 'Basal \n bistability',
              'Preferred frequency', 'Spike latency']

    # Subplots numbering
    case = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
            'N', 'O', 'P', 'Q', 'R', 'S', 'T']

    fig = plt.figure(figsize=(12, 20))
    gs = gridspec.GridSpec(10, 2)
    gs.update(wspace=0.2, hspace=0.8)
    # fig.subplots_adjust(wspace=1.5, hspace=.5)

    # Tonic spiking
    Iext = 1.5 * np.ones((2000, ))
    axx.append(run_and_plot(Iext, params, figure=fig, subplots_id=gs[0])[0])

    # Class 1
    Iext = np.ones((5000, )) + 10**-6
    axx.append(run_and_plot(Iext, params, figure=fig, subplots_id=gs[1])[0])

    # Spike frequency adaptation
    pms = copy(params)
    pms['a'] = 0.005
    Iext = np.ones((2000, )) * 2
    axx.append(run_and_plot(Iext, pms, figure=fig, subplots_id=gs[2])[0])

    # Phasic spiking
    Iext = 1.5 * np.ones((5000, ))
    ax, sol1, v1, theta_1, spk1 = run_and_plot(Iext, pms, figure=fig,
                                               subplots_id=gs[3])
    pms1 = pms.copy()
    axx.append(ax)

    # Accomodation
    Iext = np.zeros((10000,))
    Iext[:1000] = 1.5
    Iext[6000:7000] = 0.5
    Iext[7000:8000] = 1
    Iext[8000:9000] = 1.5
    axx.append(run_and_plot(Iext, pms, figure=fig, subplots_id=gs[4])[0])

    # Threshold variability
    Iext = np.zeros((4000,))
    Iext[:200] = 1.5
    Iext[2000:2250] = -1.5
    Iext[2500:2750] = 1.5
    axx.append(run_and_plot(Iext, pms, figure=fig, subplots_id=gs[5])[0])

    # Rebound spike
    Iext = np.zeros((10000,))
    Iext[500:8050] = -3.5
    axx.append(run_and_plot(Iext, pms, figure=fig, subplots_id=gs[6])[0])

    # Class 2
    Iext = 2 * (1 + 10**-6) * np.ones((3000,))
    axx.append(run_and_plot(Iext, pms, figure=fig, subplots_id=gs[7],
                            IC=(0.01, 0.001, -30.0, -30.0))[0])

    # Integrator
    Iext = np.zeros((4000,))
    Iext[:200] = 1.5
    Iext[300:500] = 1.5
    Iext[3000:3200] = 1.5
    Iext[3400:3600] = 1.5
    axx.append(run_and_plot(Iext, pms, figure=fig, subplots_id=gs[8])[0])

    # Input bistability
    Iext = np.zeros((10000,))
    Iext[:1000] = 1.5
    Iext[1000:5000] = 1.7
    Iext[5000:6000] = 1.5
    Iext[6000:] = 1.7
    axx.append(run_and_plot(Iext, pms, figure=fig, subplots_id=gs[9])[0])

    # Hyperpolarization induced spiking
    pms = copy(params)
    pms['a'] = 0.03
    Iext = -1.0 * np.ones((4000, ))
    axx.append(run_and_plot(Iext, pms, figure=fig, subplots_id=gs[10])[0])

    # Hyperpolizration induced bursting
    pms = copy(params)
    pms['a'] = 0.03
    pms['A1'] = 10
    pms['A2'] = -.6
    Iext = -1.0 * np.ones((4000,))
    axx.append(run_and_plot(Iext, pms, figure=fig, subplots_id=gs[11])[0])

    # Tonic bursting
    pms = copy(params)
    pms['a'] = 0.005
    pms['A1'] = 10
    pms['A2'] = -.6
    Iext = 2.0 * np.ones((5000,))
    axx.append(run_and_plot(Iext, pms, figure=fig, subplots_id=gs[12])[0])

    # Phasic Bursting
    pms = copy(params)
    pms['a'] = 0.005
    pms['A1'] = 10
    pms['A2'] = -.6
    Iext = 1.5 * np.ones((5000,))
    ax, sol2, v2, theta_2, spk2 = run_and_plot(Iext, pms, figure=fig,
                                               subplots_id=gs[13])
    axx.append(ax)

    # Rebound burst
    pms = copy(params)
    pms['a'] = 0.005
    pms['A1'] = 10
    pms['A2'] = -.6
    Iext = np.zeros((10000,))
    Iext[1000:6000] = -3.5
    axx.append(run_and_plot(Iext, pms, figure=fig, subplots_id=gs[14])[0])

    # Mixed mode
    pms = copy(params)
    pms['a'] = 0.005
    pms['A1'] = 5
    pms['A2'] = -0.3
    Iext = 2 * np.ones((5000,))
    axx.append(run_and_plot(Iext, pms, figure=fig, subplots_id=gs[15])[0])

    # Afterpotentials
    pms = copy(params)
    pms['a'] = 0.005
    pms['A1'] = 5
    pms['A2'] = -0.3
    Iext = np.zeros((2000,))
    Iext[:150] = 2
    Iext[150:] = 0
    axx.append(run_and_plot(Iext, pms, figure=fig, subplots_id=gs[16])[0])

    # Basal bistability
    pms = copy(params)
    pms['a'] = 0.0
    pms['A1'] = 8
    pms['A2'] = -0.1
    Iext = np.zeros((2000,))
    Iext[:100] = 5
    Iext[1000:1100] = 5
    axx.append(run_and_plot(Iext, pms, figure=fig, subplots_id=gs[17])[0])

    # Preferred frequency
    pms = copy(params)
    pms['a'] = 0.005
    pms['A1'] = -3
    pms['A2'] = 0.5
    Iext = np.zeros((8000, ))
    Iext[:50] = 5
    Iext[100:150] = 4
    Iext[4000:4050] = 5
    Iext[4500:4550] = 4
    axx.append(run_and_plot(Iext, pms, figure=fig, subplots_id=gs[18])[0])

    # Spike latency
    pms = copy(params)
    pms['a'] = -0.08
    pms['A1'] = 0
    pms['A2'] = 0
    Iext = np.zeros((500,))
    Iext[:20] = 8
    axx.append(run_and_plot(Iext, pms, figure=fig, subplots_id=gs[19])[0])

    # Figure prettification
    for i, ax in enumerate(axx):
        ax.set_ylim([-95.0, -25.0])
        text_pos_y = -12.0
        if i == 6 or i == 14:
            ax.set_ylim([-145, -25])
            text_pos_y = -2.0

        ax.text(0, text_pos_y, case[i],
                va='top',
                ha='left',
                fontsize=12,
                weight='bold',
                color='k')

        ax.set_xticklabels((ax.get_xticks() / 1000),
                           fontsize=12,
                           weight='bold')

        ax.set_yticklabels(ax.get_yticks().astype('i'),
                           fontsize=12,
                           weight='bold')

        ax.set_title(titles[i], fontsize=10, weight='bold')
        if i > 17:
            ax.set_xlabel('Time (s)', fontsize=12, weight='bold')
        if i in [0, 2, 4, 6, 8, 10,  12, 14, 16, 18]:
            ax.set_ylabel('Potential (mV)', fontsize=12,
                          weight='bold')

    # plt.savefig('Figure01.pdf', axis='tight')
    plot_figure2(v1, sol1[3], sol1[4], spk1, theta_1, pms1, savefig=False)
    plot_figure3(sol2[1], sol2[2], sol2[3], v2, spk2, savefig=False)
    plt.show()
