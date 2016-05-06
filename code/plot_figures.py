# -*- coding: utf-8 -*-
# plot_figures: Plots figures for the reference implementation of
# [1]
# [1]: "Multiple dynamical modes of thalamic relay neurons: rhythmic
# bursting and intermittent phase-locking", Wang, X-J, Neuroscience,
# 59(1), pg.21-31, 1994.
#
# Copyright (C) 2016 Georgios Is. Detorakis (gdetor@protonmail.com)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
import numpy as np
import matplotlib.cm as cmx
import matplotlib.pylab as plt
import matplotlib.colors as colors
from matplotlib import collections  as mc


def plot_figure1():
    """ Plots Figure 1. See text for more details about Figure 1.
    """
    base = "../data/"

    N, M = 10000, 20000         # Time window
    v, t = [], []               # Voltage and time lists
    # Load the data
    for i in range(5):
        v.append(np.load(base+"Fig1_V"+str(i)+".npy")[N:M])
        t.append(np.load(base+"Fig1_T"+str(i)+".npy")[N:M])

    # Plot the data
    color = ['b', 'k', 'r', 'c', 'm']
    fig = plt.figure(figsize=(10, 10))
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    for i in range(1, 6):
        ax = fig.add_subplot(5, 1, i)
        ax.plot(t[i-1], v[i-1][:, 0], color=color[i-1], lw=2.0)
        ax.set_ylim([-100, 25])
        if i == 3:
            ax.set_ylabel("Voltage (mV)", fontsize=20, weight='bold')
        if i == 5:
            ticks = ax.get_xticks()
            ticks = [(j/1000.0) for j in ticks]
            ax.set_xticklabels(ticks, fontsize=13, weight='bold')
            ax.set_xlabel("Time (s)", fontsize=16, weight='bold')
        else:
            ax.set_xticks([])
        # Configure the figure
        ax.get_xaxis().set_tick_params(which='both', direction='out')
        ax.get_yaxis().set_tick_params(which='both', direction='out')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_yticks([-100, 0, 20])
        ticks = ax.get_yticks().astype('i')
        ax.set_yticklabels(ticks, fontsize=13, weight='bold')


def plot_figure2():
    base = "../data/"
    # idx is a list with all the external currents 
    idx = [3.0, 0.0, -0.5, -0.55, -0.6, -0.8, -1.3, -2.1]

    # Load data
    v, t = [], []
    for i in idx:
        v.append(np.load(base+"Fig2_V"+str(i)+".npy"))
        t.append(np.load(base+"Fig2_T"+str(i)+".npy"))

    jet = cm = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=min(idx), vmax=max(idx))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    # Plot data
    fig = plt.figure(figsize=(12, 15))
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    for i in range(1, len(v)+1):
        ax = fig.add_subplot(8, 1, i)
        colorVal = scalarMap.to_rgba(idx[i-1])
        colorText = ('$I_{app}$: %1.2f'%(idx[i-1]))
        ax.plot(v[i-1][5000:, 0], c=colorVal, lw=1, label=colorText)
        ax.get_xaxis().set_tick_params(which='both', direction='out')
        ax.get_yaxis().set_tick_params(which='both', direction='out')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_yticks([-100, 0, 20])
        ticks = ax.get_yticks().astype('i')
        ax.set_yticklabels(ticks, fontsize=12, weight='bold')
        if i == 8:
            ticks = ax.get_xticks()
            ticks = [(j/1000.0) for j in ticks]
            ax.set_xticklabels(ticks, fontsize=13, weight='bold')
            ax.set_xlabel("Time (s)", fontsize=16, weight='bold')
        else:
            ax.set_xticks([])
        handles,labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='upper right')

def plot_figure3():
    """ Plot Figure 3, see text for more details regarding Figure 3.
    """
    N = 5000                            # Time window
    # Load data
    base = "../data/"
    v = np.load(base+"Fig3_V.npy")[N:]
    t = np.load(base+"Fig3_T.npy")[N:]

    # Plot phase plane
    fig = plt.figure(figsize=(10, 10))
    fig.subplots_adjust(wspace=0.5, hspace=.5)
    ax = plt.subplot2grid((4,2), (0, 0), colspan=2, rowspan=2)
    ax.plot(v[:, 1], v[:, 0], 'm', lw=1)
    ax.set_ylabel("V [Voltage (mV)]", fontsize=15, weight='bold')
    ax.set_xlabel("h", fontsize=15, weight='bold')
    ticks = ax.get_xticks()
    ax.set_xticklabels(ticks, fontsize=15, weight='bold')
    ticks = ax.get_yticks()
    ax.set_yticklabels(ticks, fontsize=15, weight='bold')

    # Plot membrane potential
    ax = plt.subplot2grid((4,2), (2, 0), colspan=2)
    ax.plot(t, v[:, 0], 'k', lw=1)
    ax.set_xlabel("Time (s)", fontsize=15, weight='bold')
    ax.set_ylabel("V (mV)", fontsize=15, weight='bold')
    ticks = ax.get_xticks()
    ticks = [(i/1000.0) for i in ticks]
    ax.set_xticklabels(ticks, fontsize=15, weight='bold')
    ax.set_xlabel("Time (s)", fontsize=20, weight='bold')
    ax.set_yticks([-100, 0, 20])
    ticks = ax.get_yticks().astype('i')
    ax.set_yticklabels(ticks, fontsize=15, weight='bold')
    ax.get_xaxis().set_tick_params(which='both', direction='out')
    ax.get_yaxis().set_tick_params(which='both', direction='out')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Plot h 
    ax = plt.subplot2grid((4,2), (3, 0), colspan=2)
    ax.plot(t, v[:, 1], 'r', lw=1)
    ax.set_xlabel("Time (s)", fontsize=15, weight='bold')
    ax.set_ylabel("h", fontsize=15, weight='bold')
    ticks = ax.get_xticks()
    ticks = [(i/1000.0) for i in ticks]
    ax.set_xticklabels(ticks, fontsize=15, weight='bold')
    ax.set_xlabel("Time (s)", fontsize=20, weight='bold')
    ax.set_yticks([0, 0.07])
    ticks = ax.get_yticks()
    ax.set_yticklabels(ticks, fontsize=15, weight='bold')
    ax.get_xaxis().set_tick_params(which='both', direction='out')
    ax.get_yaxis().set_tick_params(which='both', direction='out')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

if __name__ == '__main__':
    plot_figure1()
    plt.savefig('../article/figs/Figure1.pdf', axis='tight')

    plot_figure2()
    plt.savefig('../article/figs/Figure2.pdf', axis='tight')

    plot_figure3()
    plt.savefig('../article/figs/Figure3.pdf', axis='tight')

    plt.show()
