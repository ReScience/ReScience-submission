# -*- coding: utf-8 -*-
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
# from memory_profiler import profile


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
          'ThetaInf': -50.0}     # Reverse threshold


# @profile
def mnn_model(pms=params, time=10, dt=1.0, IC=(0.01, 0.001, -70.0, -50.0)):
    """ Main simulation function. Actual implementation of MNN model.

        Args:
            time (int) : Total simulation duration
            dt (float) : Euler's step
            IC (float) : Initial conditions tuple (I1, I2, V, Theta)

        Return:
            sol (ndarray)    : A numpy array containing the time (T), the two
                               intrinsic currents (I1 and I2), the membrane
                               potential (V) and the instantaneous threshold
                               potential (Theta).
            v_plot (ndarray) : The membrane potential processed for pretty
                               plots.
            P (ndarray)      : The values of the instantaneous threshold
                               potential at the spike-event times.
            spikes (ndarray) : The spike trains.
    """
    sim_time = int(time / dt)   # Total simulation time

    # Initial conditions
    I1 = IC[0] * np.ones((sim_time, ))          # Internal current 1
    I2 = IC[1] * np.ones((sim_time, ))          # Internal current 2
    V = IC[2] * np.ones((sim_time, ))           # Membrane potential
    Theta = IC[3] * np.ones((sim_time, ))       # Threshold potential
    P = np.zeros((sim_time, ))                  # Threshold track

    T, spikes = [], []          # Time and spike-event times
    # Forward Euler integration of MNN
    T.append(0)     # Time zero
    for t in range(1, sim_time):
        T.append(dt * t)        # Keep track of time
        I1[t] = I1[t-1] + dt * (-pms['k1'] * I1[t-1])
        I2[t] = I2[t-1] + dt * (-pms['k2'] * I2[t-1])
        V[t] = V[t-1] + dt * 1.0/pms['C'] * ((pms['Iext'][t-1] + I1[t-1]
                                              + I2[t-1] - pms['G'] *
                                              (V[t-1] - pms['El'])))
        Theta[t] = Theta[t-1] + dt * (pms['a'] * (V[t-1] -
                                      pms['El']) -
                                      pms['b'] *
                                      (Theta[t-1] - pms['ThetaInf']))
        # Spike event and update rules
        if V[t] >= Theta[t]:
            I1[t] = pms['R1'] * I1[t] + pms['A1']
            I2[t] = pms['R2'] * I2[t] + pms['A2']
            V[t] = pms['Vr']
            Theta[t] = max(pms['Thetar'], Theta[t])
            P[t] = Theta[t-1]
            spikes.append(t)

    # Prepare data for plotting
    V_plot = V.copy()
    pos = np.where((np.diff(V_plot) < -3.5))[0]
    V_plot[pos] = np.nan
    Theta[(P != 0)] = np.nan
    P[P == 0] = np.nan

    # Convert lists into numpy arrays
    T = np.array(T)
    spikes = np.array(spikes)

    # Pack the solution to a numpy array
    sol = np.array([T, I1, I2, V, Theta])
    return sol, V_plot, P, spikes


if __name__ == '__main__':
    import matplotlib.pylab as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pms = params
    pms['a'] = 0.005
    pms['Iext'] = 2 * (1 + 10**-8) * np.ones((3000,))
    sol, v, s, _ = mnn_model(pms, time=300, dt=0.1)       # Run a simulation
    ax.plot(sol[0], v, 'k', lw=1.5)
    ax.plot(sol[0], sol[4], 'r--', lw=1.5)
    ax.plot(sol[0], s, 'k.', lw=1.5, ms=10)
    ax.plot(sol[0], pms['Vr'] + pms['Iext'], 'k', lw=1.5, alpha=0.5)
    plt.show()
