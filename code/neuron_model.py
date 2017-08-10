#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This script implements the Mihalas-Niebur Neuron model based on the paper
# S. Mihalas and E. Niebur, "A Generalized Linear Integrate-and-Fire Neural
# Model Produces Diverse Spiking Behaviors", Neural Computation 21, 2009.
#
# Copyright (c) 2017 Georgios Is. Detorakis
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. The name of the author may not be used to endorse or promote products
#    derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
          'ThetaInf': -50.0}     # Reverse threshold


# @profile
def mnn_model(pms=params, Iext=np.zeros((200,)), dt=1.0,
              IC=(0.01, 0.001, -70.0, -50.0)):
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
            theta_ (ndarray)  : The values of the instantaneous threshold
                               potential at the spike-event times.
            spikes (ndarray)  : The spike trains.
    """
    sim_time = Iext.shape[0]                    # Simulation time

    # ODEs Initial conditions
    I1 = IC[0] * np.ones((sim_time, ))          # Internal current 1
    I2 = IC[1] * np.ones((sim_time, ))          # Internal current 2
    V = IC[2] * np.ones((sim_time, ))           # Membrane potential
    Theta = IC[3] * np.ones((sim_time, ))       # Threshold potential

    # Auxiliary vectors
    theta_ = np.zeros((sim_time, ))                 # Threshold track
    time_ = np.zeros((sim_time, ))                  # Time vector

    spikes = []                                 # Spike-event times
    # Forward Euler integration of MNN
    for t in range(1, sim_time):
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
