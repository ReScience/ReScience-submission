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


def detect_spikes(data, freq=0, dt=0.05, dist=40, th=0):
    win = int(((1.0 / freq) * 1000) / dt)
    tmp = data[4*win:5*win, 0]
    spikes = peak_local_max(tmp, min_distance=dist)
    spikes = spikes[tmp[spikes] > th]
    # Check if the model behaves properly
    if spikes.shape[0] > 4:
        raise ValueError('Bad number of spikes!')
    else:
        return spikes.shape[0]


def plot_figure1(diagram=False):
    """ Plots Figure 1. See text for more details about Figure 1.
    """
    base = "../data/"
    v, t = [], []               # Voltage and time lists
    # Load the data
    for i in range(3):
        v.append(np.load(base+"Fig1B_V"+str(i)+".npy")[:, 0])
        t.append(np.load(base+"Fig1B_T"+str(i)+".npy"))


    if diagram is True:
        fig = plt.figure(figsize=(19.5, 7.5))
        fig.subplots_adjust(wspace=0.5, hspace=.5)
    else:
        fig = plt.figure(figsize=(10, 10))

    if diagram is True:
        # Plot parameters diagram
        ax = plt.subplot2grid((3,7), (0, 0), colspan=4, rowspan=3)
        samples = 150
        freq = np.linspace(0.1, 15, samples)
        mat_ = np.empty((samples, samples))
        idx = 0
        for i in range(samples):
            for j in range(samples):
                data = np.load(base+'Fig1A_V'+str(idx)+'.npy')[0:, :]
                mat_[i, j] = detect_spikes(data, freq[i])
                idx += 1
        im = ax.imshow(mat_.T, interpolation='nearest', cmap=plt.cm.gray,
                       origin='lower', aspect='auto')
        cb = plt.colorbar(im, ticks=[0, 1, 2, 3, 4])
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
        ax.set_xticks([0, (samples-1)//4, (samples-1)//2, 3*(samples-1)//4, (samples-1)])
        ax.set_xticklabels(['0', '3.75', '7.5', '11.25', '15'],
                           fontsize=15,
                           weight='bold')
        ax.set_yticks([0, (samples-1)//2, (samples-1)])
        ax.set_yticklabels(['0', '0.5', '1.0'] ,fontsize=15, weight='bold')

    # Plot the data
    color = ['b', 'k', 'r']
    P = [100, 200, 2000]
    f = [10, 5, 0.5]         
    ms = 0.001
    p = [10, 120, 1200]
    for i in range(3):
        if diagram is True:
            ax = plt.subplot2grid((3,7), (i, 4), colspan=3)
        else:
            ax = plt.subplot2grid((3,1), (i, 0), colspan=1)
        ax.plot(t[i], v[i], color=color[i-1], lw=1.5)
        ax.add_patch(patches.Rectangle((0.0, -100.0), 200, 130, fill=False,
                     hatch='x'))
        ax.text(2000, 0.15,
                "$P_0$ = "+str(f[i])+" (Hz) \n $p/P_0$ = "+str(p[i]/P[i]),
                va='top',
                ha='left',
                fontsize=13,
                weight='bold')
        ax.set_ylim([-100, 25])
        if i == 1:
            ax.set_ylabel("Voltage (mV)", fontsize=16, weight='bold')
        if i == 2:
            ticks = ax.get_xticks()
            ticks = [(j/2000.0) for j in ticks]
            ax.set_xticklabels(ticks, fontsize=13, weight='bold')
            ax.set_xlabel("Time (s)", fontsize=16, weight='bold')
        else:
            ax.set_xticks([])
            ax.spines['bottom'].set_visible(False)
        ax.get_xaxis().set_tick_params(which='both', direction='out')
        ax.get_yaxis().set_tick_params(which='both', direction='out')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_yticks([-100, 0, 20])
        ticks = ax.get_yticks().astype('i')
        ax.set_yticklabels(ticks, fontsize=15, weight='bold')


def plot_figure2():
    base = "../data/"
    # idx is a list with all the external currents 
    idx = [3.0, 0.0, -0.47, -0.6, -0.8, -1.3, -1.4, -2.0]

    # Load data
    v, t = [], []
    for i in range(len(idx)):
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
        ax.plot(v[i-1][20000:, 0], c=colorVal, lw=1, label=colorText)
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
            ticks = [(j/10000.0) for j in ticks]
            ax.set_xticklabels(ticks, fontsize=13, weight='bold')
            ax.set_xlabel("Time (s)", fontsize=16, weight='bold')
        else:
            ax.set_xticks([])
        handles,labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='upper right')

def plot_figure3():
    """ Plot Figure 3, see text for more details regarding Figure 3.
    """
    N = 20000  # Temporal window - get rid of transient

    # Load data
    base = "../data/"
    v = np.load(base+"Fig3_V.npy")[N:]
    t = np.load(base+"Fig3_T.npy")[N:]

    # Plot phase plane
    fig = plt.figure(figsize=(10, 10))
    fig.subplots_adjust(wspace=0.5, hspace=.5)
    ax = plt.subplot2grid((4,2), (0, 0), colspan=2, rowspan=2)
    ax.plot(v[:, 1], v[:, 0], 'k', lw=1)
    ax.set_ylabel("V [Voltage (mV)]", fontsize=16, weight='bold')
    ax.set_xlabel("h", fontsize=15, weight='bold')
    ticks = ax.get_xticks()
    ax.set_xticklabels(ticks, fontsize=15, weight='bold')
    ticks = ax.get_yticks()
    ax.set_yticklabels(ticks, fontsize=15, weight='bold')

    # Plot membrane potential
    ax = plt.subplot2grid((4,2), (2, 0), colspan=2)
    ax.plot(t, v[:, 0], 'k', lw=1)
    ax.set_xticks([])
    ax.set_ylabel("V (mV)", fontsize=16, weight='bold')
    ticks = [(i/1000.0) for i in ticks]
    ax.set_yticks([-100, 0, 20])
    ticks = ax.get_yticks().astype('i')
    ax.set_yticklabels(ticks, fontsize=15, weight='bold')
    ax.get_yaxis().set_tick_params(which='both', direction='out')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Plot h 
    ax = plt.subplot2grid((4,2), (3, 0), colspan=2)
    ax.plot(t, v[:, 1], 'b', lw=1)
    ax.set_xlabel("Time (s)", fontsize=16, weight='bold')
    ax.set_ylabel("h", fontsize=16, weight='bold')
    ticks = ax.get_xticks()
    ticks = [(i/1000.0) for i in ticks]
    ax.set_xticklabels(ticks, fontsize=15, weight='bold')
    ax.set_yticks([0, 0.07])
    ticks = ax.get_yticks()
    ax.set_yticklabels(ticks, fontsize=15, weight='bold')
    ax.get_xaxis().set_tick_params(which='both', direction='out')
    ax.get_yaxis().set_tick_params(which='both', direction='out')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')


def freq_analysis(data, dt=0.05):
    from numpy.fft import rfft
    from scipy.signal import detrend, blackmanharris

    signal = data.copy()
    N = len(signal)
    signal = detrend(signal)
    windowed = signal * blackmanharris(N)
    Sf = rfft(windowed)
    Xf = np.linspace(0, 1/(2*dt), N/2)
    index = np.argmax(np.abs(Sf[:N/2]))
    return Xf[index]*1000


def plot_figure4():
    """ Plot Figure 4, see text for more details regarding Figure 4.
    """
    N = 20000  # Temporal window - get rid of transient
    base = "../data/"

    # Process the data - compute frequencies
    IextReal_A = np.genfromtxt(base+'data.dat')
    FreqReal_A = []
    for i in range(IextReal_A.shape[0]):
        data = np.load(base+'Fig4AReal_V'+str(i)+".npy")[:, 0]
        FreqReal_A.append(freq_analysis(data))

    IextReal_B = np.genfromtxt(base+'data1.dat')
    FreqReal_B = []
    for i in range(IextReal_B.shape[0]):
        data = np.load(base+'Fig4BReal_V'+str(i)+".npy")[:, 0]
        FreqReal_B.append(freq_analysis(data))

    samples = 100
    Iext = np.linspace(0, -2.0, samples)
    Freq_A, Freq_B = [], []
    for i in range(samples):
        data = np.load(base+'Fig4A_V'+str(i)+'.npy')[:, 0]
        Freq_A.append(freq_analysis(data))
        data = np.load(base+'Fig4B_V'+str(i)+'.npy')[:, 0]
        Freq_B.append(freq_analysis(data))

    I = [np.array(IextReal_A),
         np.array(IextReal_B),
         np.array(Iext),
         np.array(Iext)]
    Freq = [np.array(FreqReal_A),
            np.array(FreqReal_B),
            np.array(Freq_A),
            np.array(Freq_B)]

    # Plot the data
    fig = plt.figure(figsize=(13, 13))
    fig.subplots_adjust(wspace=0.5, hspace=.5)
    ii = 0
    for i in range(0, 4, 2):
        for j in range(0, 4, 2):
            ax1 = plt.subplot2grid((4, 4), (i, j), colspan=2, rowspan=2)
            ax2 = ax1.twinx()
            x = I[ii]
            x[x == 0] = np.nan
            y = Freq[ii]
            y[y == 0] = np.nan
            ax1.plot(x, y, 'k', alpha=0.7, lw=2.5, zorder=0)
            ax2.plot(x, 1.0 / y,'b', alpha=0.7, lw=2.5, zorder=5)
            ax1.grid(ls='--', color='k')
            ii += 1
            if i == 2:
                ax1.set_xlabel(r'$I_{app} (\mu A / cm^2)$',
                               fontsize=22,
                               weight='bold',
                               color='k')
            if j == 0:
                ax1.set_ylabel('Frequency (Hz)',
                               fontsize=16,
                               weight='bold',
                               color='k')
            if j == 2:
                ax2.set_ylabel('Period (sec)',
                               fontsize=16, 
                               weight='bold',
                               color='b')
            ax1.set_xlim([-2, 0])
            ax2.set_xlim([-2, 0])
            ax1.set_yticks(np.arange(0, 20, 2))
            ax1.set_xticks(np.round(np.linspace(0, -2.0, 5), 2))
            ax1.invert_xaxis()
            ax1.set_xticklabels(('0', '-0.5', '-1.0', '-1.5', '-2.0'),
                                fontsize=15,
                                weight='bold')

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


def plot_figure5():
    """ Plot Figure 5, see text for more details regarding Figure 5.
    """
    N = 0  # Temporal window - get rid of transient
    base = '../data/'
    stim = np.load('stimulus.npy')

    fig = plt.figure(figsize=(22, 16))
    fig.subplots_adjust(wspace=1.0, hspace=0.2)

    Iext = [-1.4, -1.5, -1.6]

    data1, stim1, sym1 = [], [], []
    for i in  range(len(Iext)):
        data1.append(np.load(base+'Fig5A_V'+str(i)+'.npy')[N:, 0])
        stim1.append(np.load(base+'Fig5A_V'+str(i)+'stimulus.npy')[N:])
        sym1.append(symbolic_analysis(data1[i]))

    d = 150
    ax = plt.subplot2grid((6,9), (0, 0), colspan=3, rowspan=3)
    for i in range(len(data1)):
        ax.plot(d+data1[i][1000:], 'k')
        ax.text(data1[0].shape[0], d-10, str(Iext[i])+'$\mu A/cm^2$',
                 va='top',
                 ha='left',
                 fontsize=13,
                 weight='bold')
        d += 150
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['top'].set_color(None)
    ax.spines['bottom'].set_visible(False)

    Iext = [-1.2, -1.5, -1.8]

    data2, stim2, sym2 = [], [], []
    for i in  range(len(Iext)):
        data2.append(np.load(base+'Fig5B_V'+str(i)+'.npy')[N:, 0])
        stim2.append(np.load(base+'Fig5B_V'+str(i)+'stimulus.npy')[N:])
        sym2.append(symbolic_analysis(data2[i]))

    d = 150
    ax = plt.subplot2grid((6,9), (3, 0), colspan=3, rowspan=3)
    for i in range(len(data2)):
        ax.plot(d+data2[i][1000:], 'b')
        d += 150
        plt.text(data2[0].shape[0], d-190, str(Iext[i])+'$\mu A/cm^2$', 
                 va='top',
                 ha='left',
                 fontsize=13,
                 weight='bold')
    ax.tick_params(axis='both',
                   which='both',
                   bottom='on',
                   top='off',
                   left='off',
                   right='off')
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_color('none')
    ax.spines['top'].set_visible(False)
    ax.set_yticks([])
    ax.set_xticks([0, (data2[0][:].shape[0] - 1)//2, (data2[0][:].shape[0]-1)])
    ax.set_xticklabels(['0', '4', '8'], fontsize=15, weight='bold')
    ax.set_xlabel('Time (s)', fontsize=16, weight='bold')

    k = 2
    M, N = 120000, 160000
    rect = [0.0, 1.0, 1.0, 0.2]
    for i in range(6):
        ax = plt.subplot2grid((6,9), (i, 3), colspan=6, rowspan=1)
        if i < 3:
            ax.plot(data1[k][M:N], 'k', lw=.8)
            ax.plot(20*stim1[k][M:N]-100, 'r', lw=.8)
            tmp = symbolic_analysis(data1[k][M:N])

            subax = add_subplot_axes(ax, rect)
            x = np.linspace(0, tmp.shape[0], tmp.shape[0])
            subax.plot(x[tmp==0], 0.5+tmp[tmp==0], 'ok', ms=11, alpha=0.6)
            subax.plot(x[tmp==1], tmp[tmp==1]-0.5, '|g', ms=15, mew=1.5, alpha=1.8)

            ax.yaxis.tick_right()
            ax.set_ylim([-160, 60])
            ax.spines['left'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
        else:
            ax.plot(data2[k][120000:160000], 'b', lw=.8)
            ax.plot(20*stim2[k][120000:160000]-100, 'r', lw=.8)
            tmp = symbolic_analysis(data2[k][M:N])

            subax = add_subplot_axes(ax, rect)
            x = np.linspace(0, tmp.shape[0], tmp.shape[0])
            subax.plot(x[tmp==0], 0.5+tmp[tmp==0], 'ok', ms=11, alpha=0.6)
            subax.plot(x[tmp==1], tmp[tmp==1]-0.5, '|g', ms=15, mew=1.5, alpha=1.8)

            ax.spines['left'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.yaxis.tick_right()
            ax.set_ylim([-160, 60])
        k -= 1
        if i == 5:
            ax.spines['bottom'].set_visible(True)
            ax.set_xticks([0, 20000, 40000])
            ax.set_xticklabels(['6', '7', '8'], fontsize=15, weight='bold')
            ax.set_xlabel('Time (s)', fontsize=16, weight='bold')
        else:
            ax.set_xticks([])
        ax.set_yticks([-100, 0, 40])
        ax.set_yticklabels(['-100', '0', '40'], fontsize=15, weight='bold')


if __name__ == '__main__':
    # plot_figure1(diagram=True)
    # plt.savefig('../article/figs/Figure1.pdf', axis='tight')

    # plot_figure2()
    # plt.savefig('../article/figs/Figure2.pdf', axis='tight')

    # plot_figure3()
    # plt.savefig('../article/figs/Figure3.pdf', axis='tight')

    # plot_figure4()
    # plt.savefig('../article/figs/Figure4.pdf', axis='tight')

    plot_figure5()
    plt.savefig('../article/figs/Figure5.pdf', axis='tight')

    plt.show()
