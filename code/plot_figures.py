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
import matplotlib.patches as patches
from matplotlib import collections  as mc
from skimage.feature import peak_local_max


def plot_figure1():
    base = '../data/'

    data = []
    data.append(np.load(base+'test_dopri_V.npy')[:, 1])   # dopri5
    data.append(np.load(base+'test_adams_V.npy')[:, 1])   # adams
    data.append(np.load(base+'test_bdf_V.npy')[:, 1])     # BDF

    labels = ['Dopri5', 'Adams', 'BDF']
    case = ['A', 'B']

    fig = plt.figure(figsize=(15, 7))
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    for i in range(2):
        ax = fig.add_subplot(1, 2, i+1, aspect=1)
        ax.plot(data[0], data[i+1], 'k.')
        plt.fill_between([-80, 20], [-82, 18], [-78, 22],
                          edgecolor='none',
                          facecolor='k', alpha=0.2)
        ax.set_xlim([-80, 20])
        ax.set_ylim([-80, 20])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ticks = ax.get_yticks()
        ax.set_yticklabels(ticks, fontsize=14, weight='bold')
        ticks = ax.get_xticks()
        ax.set_xticklabels(ticks, fontsize=14, weight='bold')
        ax.set_xlabel(labels[0]+' V (mV)', fontsize=16, weight='bold')
        ax.set_ylabel(labels[i+1]+' V (mV)', fontsize=17, weight='bold')
        ax.get_xaxis().set_tick_params(which='both', direction='out')
        ax.get_yaxis().set_tick_params(which='both', direction='out')
        if i == 0:
            ax.yaxis.set_ticks_position('left')
            ax.xaxis.set_ticks_position('bottom')
        ax.text(-78, 19, case[i],
                va='top',
                ha='left',
                fontsize=23,
                weight='bold')

    plt.savefig('../article/figs/Figure1.pdf', axis='tight')


def find_spikes(data, f, dt=0.05, th1=0, th2=-65, dist=30):
    win = int(((1.0 / f) * 1000) / dt)
    bins = int(data.shape[0] / win)
    msum = 0
    for i in range(1, bins):
        tmp = data[i*win:(i+1)*win]
        spikes = peak_local_max(tmp, min_distance=dist)
        if any(spikes > th1):
            spikes = spikes[tmp[spikes] > th1]
        else:
            spikes = spikes[tmp[spikes] > th2]
        msum += spikes.shape[0]
    return msum / bins   


def plot_figure2(diagram=False):
    """ Plots Figure 2. See text for more details about Figure 1.
    """
    base = "../data/"
    N = 10000
    v, t = [], []               # Voltage and time lists
    # Load the data
    for i in range(3):
        v.append(np.load(base+"Fig2B_V"+str(i)+".npy")[:, 1])
        t.append(np.load(base+"Fig2B_V"+str(i)+".npy")[:, 0])


    if diagram is True:
        fig = plt.figure(figsize=(8, 15))
        fig.subplots_adjust(wspace=0.5, hspace=.8)
    else:
        fig = plt.figure(figsize=(10, 10))

    if diagram is True:
        # mat_ = np.load('mat.npy')
        # Plot parameters diagram
        ax = plt.subplot2grid((8,4), (3, 0), colspan=4, rowspan=4)
        samples = 100
        freq = np.linspace(0.1, 15, samples)
        mat_ = np.empty((samples, samples))
        idx = 0
        for i in range(samples):
            period = 1./freq[i] * 1000
            dur = np.linspace(0, period, samples)
            for j in range(samples):
                data = np.load(base+'Fig2A_V'+str(idx)+'.npy')[:, 1]
                mat_[i, j] = find_spikes(data, freq[i])
                idx += 1
        # np.save('mat', mat_)
        im = ax.imshow(mat_.T, interpolation='nearest', cmap=plt.cm.hot,
                       origin='lower', aspect='auto')
        ax.contour(mat_.T, levels=[0, 0.5, 1, 2, 3, 4], linewidths=2, colors='w')
        ax.set_xlim([0, samples-1])
        ax.set_ylim([0, samples-1])
        ax.set_xlabel(r"Stimulation Frequecy $\frac{1}{P_0}$ (Hz)",
                      fontsize=16,
                      weight='bold')
        ax.set_ylabel(r"Current Pulse Duration / Period ($\frac{p}{P_0}$)",
                      fontsize=16,
                      weight='bold')
        ax.get_xaxis().set_tick_params(which='both', direction='out')
        ax.get_yaxis().set_tick_params(which='both', direction='out')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        m = mat_.shape[0]
        ax.set_xticks([0, 33, 66, 99])
        ax.set_xticklabels(['0', '5', '10', '15'],
                           fontsize=15,
                           weight='bold')
        ax.set_yticks([0, (samples-1)//2, (samples-1)])
        ax.set_yticklabels(['0', '0.5', '1.0'] ,fontsize=15, weight='bold')
        ax.text(20, 50, '4',
                va='top',
                ha='left',
                fontsize=14,
                weight='bold')
        ax.text(35, 35, '3',
                va='top',
                ha='left',
                fontsize=14,
                weight='bold')
        ax.text(42, 30, '2',
                va='top',
                ha='left',
                fontsize=14,
                weight='bold')
        ax.text(48, 23, '1',
                va='top',
                ha='left',
                fontsize=14,
                weight='bold')
        ax.text(55, 15, '0',
                va='top',
                ha='left',
                fontsize=14,
                color='w',
                weight='bold')
        ax.text(80, 65, '(0,1)',
                va='top',
                ha='left',
                fontsize=14,
                color='w',
                weight='bold')

    # Plot the data
    color = ['b', 'k', 'r']
    P = [200, 2000]
    f = [5, 0.5]         
    ms = 0.001
    p = [120, 1200]
    M, N = 0, 0
    for i in range(2):
        if diagram is True:
            ax = plt.subplot2grid((8,4), (0, i*2), colspan=2, rowspan=2)
        else:
            ax = plt.subplot2grid((4,1), (i*2, 0), colspan=1, rowspan=2)
        if i == 0:
            M, N = int(1600 / 0.05), int(2000 / 0.05)
            X, Y = 1700, 70
        if i == 1:
            M, N = int(300 / 0.05), int(3500 / 0.05)
            X, Y = 1000, 55

        ax.plot(t[i+1][M:N], v[i+1][M:N], color=color[i], lw=1.5)
        if i == 1:
            ax.set_xlim([300, 3500])
        ax.text(X, Y,
                "$1/P_0$ = "+str(f[i])+" (Hz) \n $p/P_0$ = "+str(p[i]/P[i])
                +" \n $p$ = "+str(p[i])+" (ms)",
                va='top',
                ha='left',
                fontsize=13,
                weight='bold')
        if i == 0:
            ax.set_ylim([-100, 25])
            ax.set_ylabel("Voltage (mV)", fontsize=16, weight='bold')
            ax.get_yaxis().set_tick_params(which='both', direction='out')
            ax.set_yticks([-100, 0, 20])
            ticks = ax.get_yticks().astype('i')
            ax.set_yticklabels(ticks, fontsize=15, weight='bold')
        if i == 1:
            ax.set_yticks([])
        
        if i == 0:
            ax.set_xticks([1600, 1800, 2000])
            ticks = ax.get_xticks()
            ticks = [i / 1000 for i in ticks]
            ax.set_xticklabels(ticks, fontsize=13, weight='bold')
        else:
            ticks = ax.get_xticks()
            ticks = [i / 1000 for i in ticks]
            ax.set_xticklabels(ticks, fontsize=13, weight='bold')
            ax.spines['left'].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
        ax.set_xlabel("Time (s)", fontsize=16, weight='bold')
        ax.get_xaxis().set_tick_params(which='both', direction='out')

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

    plt.savefig('../article/figs/Figure2.pdf', axis='tight')


def plot_figure3():
    """ Plots Figure 3. See text for more details about Figure 1.
    """
    dt = 0.05
    base = "../data/"
    # idx is a list with all the external currents 
    idx = [3.0, 0.0, -0.45, -0.455, -0.47, -0.55, -0.6, -0.8, -1.3, -1.4, -2.0]

    # Load data
    v, t = [], []
    for i in range(len(idx)):
        v.append(np.load(base+"Fig3_V"+str(i)+".npy"))

    # Plot data
    N = 0
    fig = plt.figure(figsize=(17, 15))
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    for i in range(1, len(v)+1):
        ax = fig.add_subplot(11, 1, i)
        ax.plot(v[i-1][N:, 0], v[i-1][N:, 1], 'k', lw=1)
        ax.text(3000, 5, '$I_{ext} = %1.3f$'%idx[i-1],
                va='top',
                ha='left',
                fontsize=18,
                weight='bold')
        ax.set_xlim([2000, 3000])
        ax.get_xaxis().set_tick_params(which='both', direction='out')
        ax.get_yaxis().set_tick_params(which='both', direction='out')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_yticks([-100, 0, 20])
        ticks = ax.get_yticks().astype('i')
        ax.set_yticklabels(ticks, fontsize=12, weight='bold')
        if i == 11:
            ticks = ax.get_xticks().astype('i')
            ticks = [j/1000 for j in ticks]
            ax.set_xticklabels(ticks, fontsize=13, weight='bold')
            ax.set_xlabel("Time (s)", fontsize=16, weight='bold')
        else:
            ax.set_xticks([])

    plt.savefig('../article/figs/Figure3.pdf', axis='tight')


def plot_figure4():
    """ Plots Figure 4. See text for more details about Figure 1.
    """
    base = '../data/'
    
    a = np.load(base+'Fig4B1_V.npy')
    s = np.load(base+'Fig4B1_Vstim.npy')
    b = np.load(base+'Fig4B2_V.npy')
    z = np.load(base+'Fig4B2_Vstim.npy')

    fig = plt.figure(figsize=(8.5, 12.5))
    fig.subplots_adjust(wspace=.5, hspace=.8)
    ax = plt.subplot2grid((5,2), (2, 0), colspan=2, rowspan=1)
    ax.plot(a[:, 0], a[:, 1], 'k', lw=2, alpha=0.5, zorder=10)
    ax.plot(b[:, 0], b[:, 1], 'k', lw=2, zorder=0)
    ax.set_ylim([-80, 20])
    ax.set_xticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ticks = ax.get_yticks().astype('i')
    ax.set_yticklabels(ticks, fontsize=13, weight='bold')
    ax.set_ylabel('V (mV)', fontsize=16, weight='bold')
    ax.get_xaxis().set_tick_params(which='both', direction='out')
    ax.get_yaxis().set_tick_params(which='both', direction='out')
    ax.text(950, 12, 'B',
            va='top',
            ha='left',
            fontsize=26,
            weight='bold')
    
    ax = plt.subplot2grid((5,2), (4, 0), colspan=2, rowspan=1)
    ax.plot(a[:, 0], a[:, 3], 'k', lw=2, alpha=0.5)
    ax.plot(b[:, 0], b[:, 3], 'k', lw=2)
    ax.set_ylim([-0.2, 0.5])
    ticks = ax.get_xticks()
    ticks = [(i/1000) for i in ticks]
    ax.set_xticklabels(ticks, fontsize=13, weight='bold')
    ax.set_xlabel('Time (s)', fontsize=16, weight='bold')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ticks = ax.get_yticks()
    ax.set_yticklabels(ticks, fontsize=13, weight='bold')
    ax.set_ylabel('H', fontsize=16, weight='bold')
    ax.get_xaxis().set_tick_params(which='both', direction='out')
    ax.get_yaxis().set_tick_params(which='both', direction='out')
    ax.text(950, 0.5, 'D',
            va='top',
            ha='left',
            fontsize=26,
            weight='bold')

    ax = plt.subplot2grid((5,2), (3, 0), colspan=2, rowspan=1)
    ax.plot(a[:, 0], s[:a.shape[0]], 'k', lw=2, alpha=0.5)
    ax.plot(a[:, 0], z[:b.shape[0]], 'k', lw=2)
    ax.set_ylim([-1.5, 0.5])
    ax.set_xticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_yticks([1.2, 0.25, 0, -0.47, -1.2])
    ticks = ax.get_yticks()
    ax.set_yticklabels(ticks, fontsize=13, weight='bold')
    ax.set_ylabel('$I_{app}\, (\mathrm{\mu A / cm^2})$', fontsize=20,
                  weight='bold')
    ax.get_xaxis().set_tick_params(which='both', direction='out')
    ax.get_yaxis().set_tick_params(which='both', direction='out')
    ax.text(950, 0.5, 'C',
            ha='left',
            fontsize=26,
            weight='bold')


    def extrema(x):
        spikes = peak_local_max(x, min_distance=150)
        return spikes

    N = 10000
    maxima, minima = [], []
    Iext = np.round(np.linspace(-0.433, -0.55, 20), 3)
    for i in range(Iext.shape[0]):
        data = np.load(base+'Fig4A_V'+str(i)+'.npy')[N:, 1]
        time = np.load(base+'Fig4A_V'+str(i)+'.npy')[N:, 0]
        ext = extrema(data)
        maxima.append(data[ext].max())
        minima.append(data[ext].min())

    maxima = np.array(maxima)
    maxima[maxima > 0] = 10
    maxima[maxima < 5] = np.nan
    minima = np.array(minima)
    minima[minima < 0] = -50
    minima[minima > 0] = np.nan

    ax = plt.subplot2grid((5,2), (0, 0), colspan=2, rowspan=2)
    ax.plot(maxima, 'k', lw=2)
    ax.plot(minima, 'k', lw=2)
    ax.plot((3, 3), (-50, 10), 'k--')
    ax.plot((11, 11), (-50, 10), 'k--')
    ax.set_xlabel(r'$I_{app}\, (\mathrm{\mu A /cm^2})$', fontsize=28)
    ax.set_ylabel(r'$V_{max}\, (\mathrm{mV})$', fontsize=28)
    ax.set_ylim([-80, 20])
    ax.set_yticks([-50, +10])
    ticks = ax.get_yticks().astype('i')
    ax.set_yticklabels(ticks, fontsize=14, weight='bold')
    ax.set_xticks(np.arange(0, 24, 4))
    ticks = [np.round(Iext[i], 3) for i in range(0, 20, 4)]
    ax.set_xticklabels(ticks, fontsize=14, weight='bold')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.get_xaxis().set_tick_params(which='both', direction='out')
    ax.get_yaxis().set_tick_params(which='both', direction='out')
    ax.arrow(3, -20, 0, -1, head_width=0.4, head_length=3, fc='k', ec='k')
    ax.arrow(11, -25, 0, +1, head_width=0.4, head_length=3, fc='k', ec='k')
    ax.arrow(6.5, 10, -0.01, 0, head_width=3.1, head_length=.45, fc='k', ec='k')
    ax.arrow(6.5, -50, 0.01, 0, head_width=3.1, head_length=.45, fc='k', ec='k')
    ax.text(3, -52, 'Subthreshold',
            va='top',
            ha='left',
            fontsize=15,
            weight='bold')
    ax.text(5, 17, 'Bursting',
            va='top',
            ha='left',
            fontsize=15,
            weight='bold')
    ax.text(19, 0, 'A',
            va='top',
            ha='left',
            fontsize=26,
            weight='bold')
    plt.savefig('../article/figs/Figure4.pdf', axis='tight')


def plot_figure5():
    """ Plot Figure 5, see text for more details regarding Figure 3.
    """
    N, M = 10000, 70000  # Temporal window - get rid of transient

    # Load data
    base = "../data/"
    t = np.load(base+"Fig5_V.npy")[N:, 0]
    v = np.load(base+"Fig5_V.npy")[N:, 1]
    h = np.load(base+"Fig5_V.npy")[N:, 2]

    # Plot phase plane
    fig = plt.figure(figsize=(8, 10))
    fig.subplots_adjust(wspace=1., hspace=1.)
    ax = plt.subplot2grid((7,3), (0, 0), colspan=3, rowspan=4)
    ax.plot(h, v, 'k', lw=1)
    ax.set_ylim([-80, 20])
    ax.set_xlim([0, 0.1])
    ax.set_ylabel("V (mV)", fontsize=16, weight='bold')
    ax.set_xlabel("h", fontsize=15, weight='bold')
    ticks = ax.get_xticks()
    ax.set_xticklabels(ticks, fontsize=15, weight='bold')
    ticks = ax.get_yticks()
    ax.set_yticklabels(ticks, fontsize=15, weight='bold')

    # Plot membrane potential
    ax = plt.subplot2grid((7,3), (4, 0), colspan=3, rowspan=2)
    ax.plot(t, v, 'k', lw=1)
    ax.set_ylim([-80, 20])
    ax.set_xticks([])
    ax.set_ylabel("V (mV)", fontsize=16, weight='bold')
    ticks = ax.get_yticks().astype('i')
    ax.set_yticklabels(ticks, fontsize=15, weight='bold')
    ax.get_yaxis().set_tick_params(which='both', direction='out')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Plot h 
    ax = plt.subplot2grid((7,3), (6, 0), colspan=3, rowspan=1)
    ax.plot(t, h, 'b', lw=1)
    ax.set_ylim([0.017, 0.077])
    ax.set_xlabel("Time (s)", fontsize=16, weight='bold')
    ax.set_ylabel("h", fontsize=16, weight='bold')
    ticks = ax.get_xticks().astype('i')
    ticks = [i/1000 for i in ticks]
    ax.set_xticklabels(ticks, fontsize=15, weight='bold')
    ax.set_yticks([0.02, 0.07])
    ticks = ax.get_yticks()
    ax.set_yticklabels(ticks, fontsize=15, weight='bold')
    ax.get_xaxis().set_tick_params(which='both', direction='out')
    ax.get_yaxis().set_tick_params(which='both', direction='out')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    plt.savefig('../article/figs/Figure5.pdf', axis='tight')


def freq_analysis(data, dt=0.05):
    from numpy.fft import rfft 
    from scipy.signal import detrend, blackmanharris

    signal = data.copy()
    N = len(signal)
    signal = detrend(signal)
    windowed = signal * blackmanharris(N)
    Sf = rfft(windowed)
    Xf = np.linspace(0, 1/(2*dt), N/2)
    # Xf = np.linspace(0, 10, N/2)
    index = np.argmax(np.abs(Sf[:N//2]))
    return Xf[index] * 1000


def plot_figure6():
    """ Plot Figure 6, see text for more details regarding Figure 4.
    """
    N, M = 10000, 30000  # Temporal window - get rid of transient
    base = "../data/"

    # Process the data - compute frequencies
    # IextReal_A = np.genfromtxt(base+'data1.dat')
    IextReal_A = np.genfromtxt(base+'data1or.dat')
    IextOr_A = np.genfromtxt(base+'data1or.dat')
    FreqReal_A, IRealA = [], []
    for i in range(IextReal_A.shape[0]):
        data = np.load(base+'Fig6AReal_V'+str(i)+".npy")[N:, 1]
        FreqReal_A.append(freq_analysis(data))

    # IextReal_B = np.genfromtxt(base+'data2.dat')
    IextReal_B = np.genfromtxt(base+'data2or.dat')
    IextOr_B = np.genfromtxt(base+'data2or.dat')
    FreqReal_B, IRealB = [], []
    for i in range(IextReal_B.shape[0]):
        data = np.load(base+'Fig6BReal_V'+str(i)+".npy")[N:, 1]
        FreqReal_B.append(freq_analysis(data))

    I = [np.array(IextReal_A).astype('f'), np.array(IextReal_B).astype('f')]
    Freq = [np.array(FreqReal_A).astype('f'), np.array(FreqReal_B).astype('f')]

    # Plot the data
    fig = plt.figure(figsize=(18, 8))
    fig.subplots_adjust(wspace=0.5, hspace=.5)
    for i in range(2):
        ax1 = plt.subplot2grid((2, 4), (0, i*2), colspan=2, rowspan=2)
        ax2 = ax1.twinx()
        x = I[i]
        xx = x.copy()
        x[x == 0] = np.nan
        y = Freq[i]
        yy = y.copy()
        y[y == 0] = np.nan
        ax1.plot(x, y, 'k.-', alpha=0.7, lw=2.5, zorder=0, ms=15)
        ax2.plot(x, 1.0 / y,'bo-', alpha=0.7, lw=2.5, zorder=5, ms = 10,
                 mfc='None', mew=2)
        if i == 0:
            ax1.scatter(IextOr_A[:, 0], IextOr_A[:, 1], s=65, c='c',
                        marker='x', linewidths=2)
            ax2.scatter(IextOr_A[:, 0], 1.0/IextOr_A[:, 1], marker='p',
                        s=65, c='m')
        if i == 1:
            ax1.scatter(IextOr_B[:, 0], IextOr_B[:, 1], marker='x', s=65,
                        c='c')
            ax2.scatter(IextOr_B[:, 0], 1.0/IextOr_B[:, 1], marker='p', s=65,
                        c='m')


        ax1.grid(ls='--', color='k')
        ax1.set_ylim([0, 20])
        ax2.set_xlim([0, -2])
        ax1.set_xlim([0, -2])
        ax2.set_ylim([0, 2])
        ax1.invert_xaxis()
        ax2.invert_xaxis()
        ticks = ax1.get_xticks()
        ax1.set_xticklabels(ticks, fontsize=15, weight='bold')
        ticks = ax2.get_xticks()
        ax2.set_xticklabels(ticks, fontsize=15, weight='bold')
        ticks = ax1.get_yticks()
        ax1.set_yticklabels(ticks, fontsize=15, weight='bold')
        ticks = ax2.get_yticks()
        ax2.set_yticklabels(ticks, fontsize=15, weight='bold')
        ax1.get_xaxis().set_tick_params(which='both',
                                        direction='in',
                                        width=1,
                                        length=10)
        ax1.get_yaxis().set_tick_params(which='both',
                                        direction='in',
                                        width=1,
                                        length=10)
        ax2.get_xaxis().set_tick_params(which='both',
                                        direction='in',
                                        width=1,
                                        length=10)
        ax2.get_yaxis().set_tick_params(which='both',
                                        direction='in',
                                        width=1,
                                        length=10)
        ax1.set_xlabel(r'$I_{app} (\mu A / cm^2)$',
                       fontsize=22,
                       weight='bold',
                       color='k')

        if i == 0:
            ax2.set_yticks([])
            ax1.set_ylabel('Frequency (Hz)',
                           fontsize=18,
                           weight='bold',
                           color='k')

        if i == 1:
            ax2.grid(ls='--', color='k')
            ax1.set_yticks([])
            ax2.set_ylabel('Period (sec)',
                            fontsize=16, 
                            weight='bold',
                            color='b')
        if i == 0:
            ax1.text(-0.05, 19, 'A',
                     va='top',
                     ha='left',
                     fontsize=22,
                     weight='bold')
        if i == 1:
            ax1.text(-0.05, 19, 'B',
                     va='top',
                     ha='left',
                     fontsize=22,
                     weight='bold')
    plt.savefig('../article/figs/Figure6.pdf', axis='tight')


def symbolic_analysis(data, freq=10, dist=10, th=0, dt=0.05):
    """ Compute the symbolic patterns for Figure 2 of the
        original article. 

        |  :param data: Array contains all the data
        |  :param freq: Frequency of input pulse current
        |  :param dist: Peak detector neighbor distance
        |  :param th  : Spike detection threshold
        |  :param dt  : Temporal integration step
        |  :return: A symbolic word contains 0, 1
    """ 
    # Definition of temporal window
    win = int(((1.0 / freq) * 1000) / dt)
    tmp = data                # Use tmp from now on

    # Peaks detection
    peaks = peak_local_max(tmp, min_distance=dist)
    a_idx = peaks[tmp[peaks] > th]

    # Remove local maxima or minima that does not correspond to 
    # neither supra- not to subthreshold spikes
    peaks_r = []
    for i in range(1, peaks.shape[0]):
        if (tmp[peaks[i-1]] > th) and (tmp[peaks[i]] < th):
            d = np.abs(peaks[i-1] - peaks[i])
            if d < 300:
                peaks_r.append(peaks[i])

    def fun(a, b):
        """ Check for matching indices and return them. """
        idx = []
        for j, aa in enumerate(a):
            if aa not in b:
                idx.append(j)
        return np.array(idx)

    peaks_r = np.array(peaks_r)
    # If no ghost spikes appear in peaks_r 
    if peaks_r.shape[0] == 0:
        z_idx = peaks[tmp[peaks] < th]
    # If ghost spikes exist
    else:
        idxs = fun(peaks, peaks_r)      # Check for matching indices
        peaks = peaks[idxs]
        z_idx = peaks[tmp[peaks] < th]

    word = -1 * np.ones((int(dt*tmp.shape[0]), ))     # Initalization of symbolic word
    if len(a_idx) == 0:
        word[(dt*z_idx).astype('i')] = 0
        return word[word != -1]
    elif len(z_idx)==0:
        word[(dt*a_idx).astype('i')] = 1
        return word[word != -1]
    else:
        word[(dt*z_idx).astype('i')] = 0
        word[(dt*a_idx).astype('i')] = 1
        return word[word != -1]


def spike_train(data):
    tmp = np.zeros(data.shape)
    spikes = peak_local_max(data, min_distance=20)
    spikes = spikes[data[spikes] > 0]
    tmp[spikes] = 1
    tmp[tmp==0] = np.nan
    return tmp


def add_subplot_axes(ax, rect, axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    subax.set_xticks([])
    subax.set_yticks([])
    subax.spines['right'].set_visible(False)
    subax.spines['left'].set_visible(False)
    subax.spines['top'].set_visible(False)
    subax.spines['bottom'].set_visible(False)
    return subax


def plot_figure7():
    """ Plot Figure 7, see text for more details regarding Figure 5.
    """
    N = 0  # Temporal window - get rid of transient
    base = '../data/'
    stim = np.load(base+'Fig7A_V0stim.npy')

    fig = plt.figure(figsize=(22, 16))
    fig.subplots_adjust(wspace=1.0, hspace=0.2)

    Iext = [-1.4, -1.5, -1.6, -1.2, -1.5, -1.8]
    time1, data1, stim1, sym1 = [], [], [], []
    for i in  range(3):
        data1.append(np.load(base+'Fig7A_V'+str(i)+'.npy')[N:, 1])
        time1.append(np.load(base+'Fig7A_V'+str(i)+'.npy')[N:, 0])
        stim1.append(np.load(base+'Fig7A_V'+str(i)+'stim.npy')[N:])
        sym1.append(symbolic_analysis(data1[i]))

    time2, data2, stim2, sym2 = [], [], [], []
    for i in  range(3):
        data2.append(np.load(base+'Fig7B_V'+str(i)+'.npy')[N:, 1])
        time2.append(np.load(base+'Fig7B_V'+str(i)+'.npy')[N:, 0])
        stim2.append(np.load(base+'Fig7B_V'+str(i)+'stim.npy')[N:])
        sym2.append(symbolic_analysis(data2[i]))

    d = 150
    k = 0
    K, L = 80000, 120000
    M, N = 120000, 160000
    rect = [0.0, 0.8, 1.0, 0.2]
    for i in range(3):
        ax = plt.subplot2grid((6,6), (i, 0), colspan=6, rowspan=1)

        if i == 0:
            K, L = 196000, 236000
        if i == 1:
            K, L = 428000, 468000
        if i == 2:
            K, L = 220000, 260000
        ax.plot(data1[k][K:L], 'k', lw=.8)
        ax.plot(20*stim1[k][K:L]-100, 'r', lw=.8)
        tmp = symbolic_analysis(data1[k][K:L])

        subax = add_subplot_axes(ax, rect)
        x = np.linspace(0, tmp.shape[0], tmp.shape[0])
        subax.plot(x[tmp==0], 0.5+tmp[tmp==0], 'ok', ms=11, alpha=0.6)
        subax.plot(x[tmp==1], tmp[tmp==1]-0.5, '|g', ms=15, mew=1.5, alpha=1.8)

        # ax.yaxis.tick_right()
        ax.set_ylim([-160, 80])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.set_xticks([])
        ax.text(40000, d-150, str(Iext[i])+'$\mathrm{\mu A/cm^2}$', 
                 va='top',
                 ha='left',
                 fontsize=17,
                 weight='bold')

        if i == 0:
            M, N = 60000, 100000
        if i == 1:
            M, N = 50000, 90000
        ax1 = plt.subplot2grid((6,6), (i+3, 0), colspan=6, rowspan=1)
        ax1.plot(data2[k][M:N], 'b', lw=.8)
        ax1.plot(20*stim2[k][M:N]-100, 'r', lw=.8)
        tmp = symbolic_analysis(data2[k][M:N])

        subax1 = add_subplot_axes(ax1, rect)
        x = np.linspace(0, tmp.shape[0], tmp.shape[0])
        subax1.plot(x[tmp==0], 0.5+tmp[tmp==0], 'ok', ms=11, alpha=0.6)
        subax1.plot(x[tmp==1], tmp[tmp==1]-0.5, '|g', ms=15, mew=1.5, alpha=1.8)

        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        # ax1.yaxis.tick_right()
        ax1.set_ylim([-160, 80])
        ax1.set_xticks([])
        k += 1
        if i == 2:
            ax1.spines['bottom'].set_visible(True)
            ax1.set_xticks([0, 20000, 40000])
            ax1.set_xticklabels(['0', '1', '2'], fontsize=15, weight='bold')
            ax1.set_xlabel('Time (s)', fontsize=16, weight='bold')
        ax.set_yticks([-100, 0, 40])
        ax.set_yticklabels(['-100', '0', '40'], fontsize=15, weight='bold')
        ax1.set_yticks([-100, 0, 40])
        ax1.set_yticklabels(['-100', '0', '40'], fontsize=15, weight='bold')

        ax1.text(40000, d-150, str(Iext[i+3])+'$\mathrm{\mu A/cm^2}$', 
                 va='top',
                 ha='left',
                 fontsize=17,
                 weight='bold')

    plt.savefig('../article/figs/Figure7.pdf', axis='tight')


if __name__ == '__main__':
    plot_figure1()
    plot_figure2(diagram=True)
    plot_figure3()
    plot_figure4()
    plot_figure5()
    plot_figure6()
    plot_figure7()
    plt.show()
