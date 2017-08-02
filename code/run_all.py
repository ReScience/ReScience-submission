# This script implements the Mihalas-Niebur Neuron model based on the paper
# S. Mihalas and E. Niebur, "A Generalized Linear Integrate-and-Fire Neural
# Model Produces Diverse Spiking Behaviors", Neural Computation 21, 2009.
#
# Copyright (C) 2017  Georgios Is. Detorakis (gdetor@protonmail.com)
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details. You should have received a copy of the GNU General Public
# License along with this program; if not, write to
# the Free Software Foundation, Inc.,
# 51 Franklin Street,
# Fifth Floor, Boston, MA  02110-1301, USA.
import numpy as np
from copy import copy
from neuron_model import mnn_model
import matplotlib.pylab as plt


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
          'Iext': np.zeros((200,)),    # External current
          'ThetaInf': -50.0}       # Reverse threshold


if __name__ == '__main__':
    """
        This script plots all the figures in the paper (Figures 1, 2, 3).
    """
    axx = []

    # Subplots titles
    titles = ['Tonic spiking', 'Class 1', 'Spike frequency adaptation',
              'Phasic spiking', 'Accomodation', 'Threshold variability',
              'Rebound spike', 'Class 2', 'Integrator', 'Input bistability',
              'Hyperpolarization \n induced spiking',
              'Hyperpolarization \n induced bursting', 'Tonic bursting',
              'Phasic bursting', 'Rebound burst',
              'Mixed mode', 'Afterpotentials', 'Basal bistability',
              'Preferred frequency', 'Spike latency']

    case = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
            'N', 'O', 'P', 'Q', 'R', 'S', 'T']

    # Figure 1
    fig = plt.figure(figsize=(17, 11))
    fig.subplots_adjust(wspace=.5, hspace=.5)

    # Tonic spiking
    ax = fig.add_subplot(5, 4, 1)
    axx.append(ax)
    pms = copy(params)
    pms['Iext'] = 1.5 * np.ones((2000, ))
    sol, v, s, _ = mnn_model(pms, time=200, dt=0.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    # Class 1
    ax = fig.add_subplot(5, 4, 2)
    axx.append(ax)
    pms = copy(params)
    pms['Iext'] = np.ones((5000, )) + 10**-6
    sol, v, s, _ = mnn_model(pms, time=500, dt=0.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    # Spike frequency adaptation
    ax = fig.add_subplot(5, 4, 3)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.005
    pms['Iext'] = np.ones((2000, )) * 2
    sol, v, s, _ = mnn_model(pms, time=200, dt=0.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    # Phasic spiking
    ax = fig.add_subplot(5, 4, 4)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.005
    pms['Iext'] = 1.5 * np.ones((5000, ))
    sol, v, p, s = mnn_model(pms, time=500, dt=0.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], p, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)
    v_ = v.copy()
    x1 = sol[3].copy()
    y1 = sol[4].copy()
    pp = p.copy()
    ss = s.copy()

    # Accomodation
    ax = fig.add_subplot(5, 4, 5)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.005
    pms['Iext'] = np.zeros((10000,))
    pms['Iext'][:1000] = 1.5
    pms['Iext'][6000:7000] = 0.5
    pms['Iext'][7000:8000] = 1
    pms['Iext'][8000:9000] = 1.5
    sol, v, s, _ = mnn_model(pms, time=1000, dt=0.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    # Threshold variability
    ax = fig.add_subplot(5, 4, 6)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.005
    pms['Iext'] = np.zeros((4000,))
    pms['Iext'][:200] = 1.5
    pms['Iext'][2000:2250] = -1.5
    pms['Iext'][2500:2750] = 1.5
    sol, v, s, _ = mnn_model(pms, time=400, dt=0.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    # Rebound spike
    ax = fig.add_subplot(5, 4, 7)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.005
    pms['Iext'] = np.zeros((10000,))
    pms['Iext'][500:8050] = -3.5
    sol, v, s, _ = mnn_model(pms, time=1000, dt=0.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    # Class 2
    ax = fig.add_subplot(5, 4, 8)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.005
    pms['Iext'] = 2 * (1 + 10**-6) * np.ones((3000,))
    sol, v, s, _ = mnn_model(pms, time=300, dt=0.1, IC=(0.01, 0.001, -30, -30))
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    # Integrator
    ax = fig.add_subplot(5, 4, 9)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.005
    pms['Iext'] = np.zeros((4000,))
    pms['Iext'][:200] = 1.5
    pms['Iext'][300:500] = 1.5
    pms['Iext'][3000:3200] = 1.5
    pms['Iext'][3400:3600] = 1.5
    sol, v, s, _ = mnn_model(pms, time=400, dt=0.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    # Input bistability
    ax = fig.add_subplot(5, 4, 10)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.005
    pms['Iext'] = np.zeros((10000,))
    pms['Iext'][:1000] = 1.5
    pms['Iext'][1000:5000] = 1.7
    pms['Iext'][5000:6000] = 1.5
    pms['Iext'][6000:] = 1.7
    sol, v, s, _ = mnn_model(pms, time=1000, dt=0.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    # Hyperpolarization induced spiking
    ax = fig.add_subplot(5, 4, 11)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.03
    pms['Iext'] = -1.0 * np.ones((4000, ))
    sol, v, s, _ = mnn_model(pms, time=400, dt=0.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    # Hyperpolizration induced bursting
    ax = fig.add_subplot(5, 4, 12)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.03
    pms['A1'] = 10
    pms['A2'] = -.6
    pms['Iext'] = -1.0 * np.ones((4000,))
    sol, v, s, _ = mnn_model(pms, time=400, dt=0.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    # Tonic bursting
    ax = fig.add_subplot(5, 4, 13)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.005
    pms['A1'] = 10
    pms['A2'] = -.6
    pms['Iext'] = 2.0 * np.ones((5000,))
    sol, v, s, _ = mnn_model(pms, time=500, dt=0.1)       # Run a simulation
    ax.plot(sol[0], sol[3], 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    # Phasic Bursting
    ax = fig.add_subplot(5, 4, 14)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.005
    pms['A1'] = 10
    pms['A2'] = -.6
    pms['Iext'] = 1.5 * np.ones((5000,))
    sol, v, p, s = mnn_model(pms, time=500, dt=0.1)       # Run a simulation
    ax.plot(sol[0], sol[3], 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], p, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)
    x2 = sol[1].copy()
    y2 = sol[2].copy()
    z2 = sol[3].copy()
    z = v.copy()
    ss2 = s.copy()

    # Rebound burst
    ax = fig.add_subplot(5, 4, 15)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.005
    pms['A1'] = 10
    pms['A2'] = -.6
    pms['Iext'] = np.zeros((10000,))
    pms['Iext'][1000:6000] = -3.5
    sol, v, s, _ = mnn_model(pms, time=1000, dt=0.1)       # Run a simulation
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)
    ax.plot(sol[0], sol[3], 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)

    # Mixed mode
    ax = fig.add_subplot(5, 4, 16)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.005
    pms['A1'] = 5
    pms['A2'] = -0.3
    pms['Iext'] = 2 * np.ones((5000,))
    sol, v, s, _ = mnn_model(pms, time=500, dt=0.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    # Afterpotentials
    ax = fig.add_subplot(5, 4, 17)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.005
    pms['A1'] = 5
    pms['A2'] = -0.3
    pms['Iext'] = np.zeros((2000,))
    pms['Iext'][:150] = 2
    pms['Iext'][150:] = 0
    sol, v, s, _ = mnn_model(pms, time=200, dt=0.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    # Basal bistability
    ax = fig.add_subplot(5, 4, 18)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.0
    pms['A1'] = 8
    pms['A2'] = -0.1
    pms['Iext'] = np.zeros((2000,))
    pms['Iext'][:100] = 5
    pms['Iext'][1000:1100] = 5
    sol, v, s, _ = mnn_model(pms, time=200, dt=.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    # Preferred frequency
    ax = fig.add_subplot(5, 4, 19)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = 0.005
    pms['A1'] = -3
    pms['A2'] = 0.5
    pms['Iext'] = np.zeros((8000, ))
    pms['Iext'][:50] = 5
    pms['Iext'][100:150] = 4
    pms['Iext'][4000:4050] = 5
    pms['Iext'][4500:4550] = 4
    sol, v, s, _ = mnn_model(pms, time=800, dt=.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    # Spike latency
    ax = fig.add_subplot(5, 4, 20)
    axx.append(ax)
    pms = copy(params)
    pms['a'] = -0.08
    pms['A1'] = 0
    pms['A2'] = 0
    pms['Iext'] = np.zeros((500,))
    pms['Iext'][:20] = 8
    sol, v, s, _ = mnn_model(pms, time=50, dt=0.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)

    for i, ax in enumerate(axx):
        ax.set_ylim([-75, -35])
        if i == 5:
            ax.set_ylim([-95, -35])
        if i == 6:
            ax.set_ylim([-145, -35])
        if i == 7:
            ax.set_ylim([-75, -20])
        if i == 10:
            ax.set_ylim([-95, -35])
        if i == 11:
            ax.set_ylim([-115, -35])
        if i == 14:
            ax.set_ylim([-145, -35])
        if i == 17:
            ax.set_ylim([-95, -35])
        ax.set_xticklabels((ax.get_xticks() / 1000),
                           fontsize=12,
                           weight='bold')
        ax.set_yticklabels(ax.get_yticks().astype('i'),
                           fontsize=12,
                           weight='bold')
        ax.set_title(titles[i], fontsize=12, weight='bold')
        if i > 15:
            ax.set_xlabel('Time (s)', fontsize=12, weight='bold')
        if i in [0, 4, 8, 12, 16]:
            ax.set_ylabel('Potential (mV)', fontsize=12,
                          weight='bold')
        ax.text(0, -37, case[i],
                va='top',
                ha='left',
                fontsize=12,
                weight='bold',
                color='orange')
    # plt.savefig('Figure01.pdf', axis='tight')

    # ----------------------------------------------------------------------
    # Figure 2 -- Phase space of phasic spiking simulation
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(v_, y1)
    ax.plot(x1[ss-1], y1[ss-1], 'k.', ms=10)
    p = np.polyfit(x1[ss-1], y1[ss-1], deg=1)
    x_the = np.linspace(-75, -25, 100)
    y_the = np.poly1d(p)(x_the)
    ax.plot(x_the, y_the, 'k--', alpha=0.5)
    ax.set_xlim([-75, -28])
    ax.set_ylim([-55, -28])
    ax.set_title('Phasic spiking', fontsize=16, weight='bold')
    ax.set_xlabel('Membrane potential (mV)', fontsize=16, weight='bold')
    ax.set_ylabel('Threshold potential (mV)', fontsize=16, weight='bold')
    ax.set_xticklabels(ax.get_xticks().astype('i'), fontsize=14, weight='bold')
    ax.set_yticklabels(ax.get_yticks().astype('i'), fontsize=14, weight='bold')
    # plt.savefig('Figure02.pdf', axis='tight')

    # ----------------------------------------------------------------------
    # Figure 2 -- phase space of phasic bursting simulation
    fig = plt.figure(figsize=(11, 11))
    from mpl_toolkits.mplot3d import Axes3D
    ax = fig.add_subplot(111)
    ax = fig.gca(projection='3d')
    ax.plot(x2, y2, z)
    ax.plot(x2[ss2-1], y2[ss2-1], z2[ss2-1], 'k.', ms=15)
    ax.set_xlim([0, 10])
    ax.set_ylim([-3, 0])
    ax.set_title('Phasic bursting', fontsize=16, weight='bold')
    ax.set_ylabel('I2 (V/s)', fontsize=16, weight='bold')
    ax.set_xlabel('I1 (V/s)', fontsize=16, weight='bold')
    ax.set_zlabel('Membrane potential (mV)', fontsize=16, weight='bold')
    ax.set_xticklabels(ax.get_xticks().astype('i'), fontsize=14, weight='bold')
    ax.xaxis.labelpad = 20
    ax.set_yticklabels(ax.get_yticks(), fontsize=14, weight='bold')
    ax.yaxis.labelpad = 30
    ax.set_zticklabels(ax.get_zticks().astype('i'), fontsize=14, weight='bold')
    ax.zaxis.labelpad = 30
    # plt.savefig('Figure03.pdf', axis='tight')

    plt.show()
