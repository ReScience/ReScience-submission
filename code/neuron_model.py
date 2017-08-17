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
# from memory_profiler import profile


def default_parameters(**kwargs):
    """ Sets and stores the default values for simulations in a dictionary.

        Args:
            **kwargs            : Keyword arguments for defining specific
                                  values for the parameters

        Return:
            params (dict)       : Simulation parameters dictionary
    """
    params = dict()
    k1, k2 = 0.2, 0.02
    a, b = 0.0, 0.01
    A1, A2 = 0.0, 0.0
    R1, R2 = 0.0, 1.0
    G, C = 0.05, 1.0
    El, Vr = -70.0, -70.0
    ThetaR, ThetaInf = -60.0, -50.0

    params['k1'] = k1
    params['k2'] = k2
    params['a'] = a
    params['b'] = b
    params['A1'] = A1
    params['A2'] = A2
    params['R1'] = R1
    params['R2'] = R2
    params['El'] = El
    params['Vr'] = Vr
    params['Thetar'] = ThetaR
    params['ThetaInf'] = ThetaInf
    params['G'] = G
    params['C'] = C

    if kwargs is not None:
        for key, value in kwargs.items():
            if key in params:
                params[key] = value
            else:
                print("Key {} not found! Default values are used!".format(key))
    return params


# @profile
def mnn_model(pms=None, Iext=None, dt=1.0,
              IC=(0.01, 0.001, -70.0, -50.0)):
    """ Main simulation function. Actual implementation of MNN model.

        Args:
            params (dict)   : Parameters dictionary
            Iext (ndarray)  : External current array
            dt (float)      : Euler's step
            IC (float)      : Initial conditions tuple (I1, I2, V, Theta)

        Return:
            sol (ndarray)    : A numpy array containing the time (T), the two
                               intrinsic currents (I1 and I2), the membrane
                               potential (V) and the instantaneous threshold
                               potential (Theta).
            theta_ (ndarray)  : The values of the instantaneous threshold
                               potential at the spike-event times.
            spikes (ndarray)  : The spike trains.
    """
    if Iext is None:
        Iext = np.zeros((200, ))

    sim_steps = Iext.shape[0]                    # Total simulation steps

    if pms is None:
        pms = default_parameters()

    # ODEs Initial conditions
    I1 = IC[0] * np.ones((sim_steps, ))          # Internal current 1
    I2 = IC[1] * np.ones((sim_steps, ))          # Internal current 2
    V = IC[2] * np.ones((sim_steps, ))           # Membrane potential
    Theta = IC[3] * np.ones((sim_steps, ))       # Threshold potential

    # Auxiliary vectors
    theta_ = np.zeros((sim_steps, ))                 # Threshold track
    time_ = np.zeros((sim_steps, ))                  # Time vector

    spikes = []                                 # Spike-event times
    # Forward Euler integration of MNN
    for t in range(1, sim_steps):
        # ODEs
        I1[t] = I1[t-1] + dt * (-pms['k1'] * I1[t-1])
        I2[t] = I2[t-1] + dt * (-pms['k2'] * I2[t-1])
        V[t] = V[t-1] + dt * 1.0/pms['C'] * ((Iext[t-1] + I1[t-1]
                                              + I2[t-1] - pms['G'] *
                                              (V[t-1] - pms['El'])))
        Theta[t] = Theta[t-1] + dt * (pms['a'] * (V[t-1] -
                                      pms['El']) -
                                      pms['b'] *
                                      (Theta[t-1] - pms['ThetaInf']))

        # Spike event and update rules
        if V[t] >= Theta[t]:
            # Updates
            I1[t] = pms['R1'] * I1[t] + pms['A1']
            I2[t] = pms['R2'] * I2[t] + pms['A2']
            V[t] = pms['Vr']
            Theta[t] = max(pms['Thetar'], Theta[t])

            # Auxiliary variables - Threshold in previous time-step
            theta_[t] = Theta[t-1]

            # Spike events
            spikes.append(t)

        # Keep track of time
        time_[t] = dt * t

    # Convert lists into numpy arrays
    spikes = np.array(spikes)

    # Pack the solution to a numpy array
    sol = np.array([time_, I1, I2, V, Theta])
    return sol, theta_, spikes
