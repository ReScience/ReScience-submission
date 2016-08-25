# -*- coding: utf-8 -*-
# 
# neuron_model:
# This is a reference implementation of: "Multiple dynamical modes
# of thalamic relay neurons: rhythmic bursting and intermittent
# phase-locking", Wang, X-J, Neuroscience, 59(1), pg.21-31, 1994.
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
from scipy.integrate import ode
from skimage.feature import peak_local_max


def loadParameters(fname):
    """ Load all the necessary paremeters from a file.

        |  :param fname: File name containing all parameters
        |  :return: A dictionary with all the parameters
    """
    import os
    if os.path.isfile(fname) is False:
        print('File {} does not exist!'.format(fname))
    else:
        import configparser

        config = configparser.RawConfigParser()
        config.read(fname)

        p = {
            'V0': config.getfloat('Membrane', 'V0'),
            'phi_h': config.getfloat('Ttype', 'phi_h'),
            'gT': config.getfloat('Ttype', 'gT'),
            'VCa': config.getfloat('Ttype', 'VCa'),
            'theta_h': config.getfloat('Ttype', 'theta_h'),
            'k_h': config.getfloat('Ttype', 'k_h'),
            'phi_H': config.getfloat('Sag', 'phi_H'),
            'gH': config.getfloat('Sag', 'gH'),
            'VH': config.getfloat('Sag', 'VH'),
            'phi_K': config.getfloat('Potassium', 'phi_K'),
            'sigma_K': config.getfloat('Potassium', 'sigma_K'),
            'gK': config.getfloat('Potassium', 'gK'),
            'VK': config.getfloat('Potassium', 'VK'),
            'gNa': config.getfloat('Sodium', 'gNa'),
            'VNa': config.getfloat('Sodium', 'VNa'),
            'sigma_Na': config.getfloat('Sodium', 'sigma_Na'),
            'gNaP': config.getfloat('PersistentSodium', 'gNaP'),
            'VNaP': config.getfloat('PersistentSodium', 'VNaP'),
            'sigma_NaP': config.getfloat('PersistentSodium', 'sigma_NaP'),
            'gL': config.getfloat('Leak', 'gL'),
            'VL': config.getfloat('Leak', 'VL'),
            'Cm': config.getfloat('Neuron', 'Cm')
        }
        return p


class WangModel():
    def __init__(self, fname,
                 dt=0.05,   
                 frequency=1,
                 duration=1000,
                 inp_type='const',
                 amplitude=1,
                 pulse_ampl=0,
                 chunks=1,
                 store=True,
                 store_stim=False):
        """ Constructor method

        |  :param fname: File name containing all parameters
        |  :param dt: Time discretization step
        |  :param frequency: Stimulation frequency (in Hz)
        |  :param duration: Duration (ON) of stimulation (in ms)
        |  :param inp_type: Input signal type ['const' - constant current
                                               'periodic' - periodic pulses]
                                               'pulse' - one pulse at the 
                                               beginning of simulation
                                               following by a constant current.
        |  :param amplitude: Amplitude of injected current (in μA/cm²)
        |  :param pulse_ampl: Amplitude of initial pulse (input must be set to
                                                          'pulse')
        |  :param chunks: Number of repeated periods of pulses
        |  :param store: Store results of simulation on npy files
        |  :param store_stim: Store the stimulus signal
        |  :return:
        """
        self.p = loadParameters(fname)      # Load the parameters

        self.dt = dt                        # Integration time step
        self.freq = frequency               # Pulse current frequency
        self.P0 = (1.0/self.freq)*1000      # Period in ms
        self.P = duration                   # Duration in ms
        self.amplitude = amplitude          # Amplitude in μA/cm²
        self.pulse_amplitude = pulse_ampl   # Amplitude of a pulse
        self.chunks = chunks                # Number of perios (time bins)
        self.inp_type = inp_type            # Input type (constant or pulse)
        self.IsStore = store                # Store neuron activity and time
        self.IsStoreStim = store_stim       # Store stimulus

    def s_inf(self, v):
        """ T-type calcium channel gating kinetics """
        return 1.0/(1.0 + np.exp(-(v + 65.) / 7.8))

    def h_inf(self, v, theta, k):
        """ T-type calcium channel gating kinetics """
        return 1.0/(1.0 + np.exp((v - theta) / k))

    def tau_h(self, v, theta, k):
        """ T-type calcium channel gating kinetics """
        return self.h_inf(v, theta, k)*np.exp((v + 162.3) / 17.8) + 20.0

    def H_inf(self, v):
        """ Sag channel gating kinetics """
        return 1.0/(1.0 + np.exp((v + 69.0) / 7.1))

    def tau_H(self, v):
        """ Sag channel gating kinetics """
        return 1000.0/(np.exp((v + 66.4)/9.3) + np.exp(-(v + 81.6)/13.))

    def alpha_n(self, v, sigma):
        """ Potassium (K) channel gating kinetics """
        return (-0.01*(v + 45.7 - sigma)/(np.exp(-0.1 *
                (v + 45.7 - sigma)) - 1.))

    def beta_n(self, v, sigma):
        """ Potassium (K) channel gating kinetics """
        return 0.125*np.exp(-(v + 55.7 - sigma)/80.0)

    def n_inf(self, v, sigma):
        """ Potassium (K) channel gating kinetics """
        return (self.alpha_n(v, sigma) /
                (self.alpha_n(v, sigma) + self.beta_n(v, sigma)))

    def tau_n(self, v, sigma):
        """ Potassium (K) channel gating kinetics """
        return 1.0/(self.alpha_n(v, sigma) + self.beta_n(v, sigma))

    def alpha_m(self, v, sigma):
        """ Sodium (Na) channel gating kinetics """
        return (-0.1*(v + 29.7 - sigma)/(np.exp(-0.1 *
                (v + 29.7 - sigma)) - 1.))

    def beta_m(self, v, sigma):
        """ Sodium (Na) channel gating kinetics """
        return 4.*np.exp(-(v + 54.7 - sigma)/18.0)

    def m_inf(self, v, sigma):
        """ Sodium (Na) channel gating kinetics """
        return (self.alpha_m(v, sigma)/(self.alpha_m(v, sigma) +
                self.beta_m(v, sigma)))

    def I_T(self, v, h):
        """
        Membrane current (in μA/cm²)
        T-Type Calcium

        |  :param V: Membrane voltage (mV)
        |  :param h: Gating variable
        |  :return: T-Type current
        """
        return self.p['gT'] * self.s_inf(v)**3 * h * (v - self.p['VCa'])

    def I_H(self, v, H):
        """
        Membrane current (in μA/cm²)
        Sag

        |  :param V: Membrane voltage (mV)
        |  :param H: Gating variable
        |  :return: Sag current
        """
        return self.p['gH'] * H**2 * (v - self.p['VH'])

    def I_K(self, v, n):
        """
        Membrane current (in μA/cm²)
        Potassium (K)

        |  :param V: Membrane voltage (mV)
        |  :param n: Gating variable
        |  :return: Potassium current
        """
        return self.p['gK'] * n**4 * (v - self.p['VK'])

    def I_Na(self, v, n):
        """
        Membrane current (in μA/cm²)
        Sodium (Na)

        |  :param v: Membrane voltage (mV)
        |  :param n: Gating variable
        |  :return: Sodium current
        """
        return (self.p['gNa'] * self.m_inf(v, self.p['sigma_Na'])**3 *
                (0.85 - n) * (v - self.p['VNa']))

    def I_NaP(self, v):
        """
        Membrane current (in μA/cm²)
        Persistent Sodium (Na)

        |  :param v: Membrane voltage (mV)
        |  :param h: Gating variable
        |  :return: Persistent sodium current
        """
        return (self.p['gNaP'] * self.m_inf(v, self.p['sigma_NaP'])**3 *
                (v - self.p['VNa']))

    def I_L(self, v):
        """
        Membrane current (in μA/cm²)
        Leak current

        |  :param v: Membrane voltage (mV)
        |  :return: Leak current
        """
        return self.p['gL'] * (v - self.p['VL'])

    def I_ext(self, period, duration, time):
        """
        External current (in μA/cm²)

        |  :param period: Period of current on pulses
        |  :param duration: Duration of current on pulses
        |  :param amplitude: Pulse amplitude
        |  :param time: Total simulation time (ms)
        |  :param numChunks: Number of chunks according to period 
        """
        tmp = np.zeros((period, ))
        tmp[:duration] = self.amplitude
        # repeats = int(self.chunks + (time - period) / period)
        repeats = int(self.chunks + (time / period))
        return np.tile(tmp, repeats)

    def rhs(self, t, X, Iapp):
        """
        Right-hand side of Wang [1] model.

        |  :param X: States of the system
        |  :param t: Time
        |  :param Iapp: Injected current (stimulation)
        |  :return: Calculate membrane potential & activation variables
        """
        V, h, H, n = X[0], X[1], X[2], X[3]

        dV = ((-self.I_T(V, h) - self.I_H(V, H) - self.I_Na(V, n) -
              self.I_K(V, n) - self.I_NaP(V) - self.I_L(V) +
              Iapp)/self.p['Cm'])
        dh = (self.p['phi_h']*(self.h_inf(V, self.p['theta_h'], self.p['k_h'])
              - h)/self.tau_h(V, self.p['theta_h'], self.p['k_h']))
        dH = self.p['phi_H'] * (self.H_inf(V) - H) / self.tau_H(V)
        dn = (self.p['phi_K'] * (self.n_inf(V, self.p['sigma_K']) - n) /
              self.tau_n(V, self.p['sigma_K']))

        return np.array([dV, dh, dH, dn])

    def run(self, tf, vname, int_name='vode', *args):
        """ 
        Run the model by integrating equations.

        |  :param tf: Total time (ms)
        |  :param vname: Filename for storing time and membrane potential
                         Time - first column, Neuron states - following four
                         columns. [Second column - membrane potential (V)
                                   Third column - Inactivation variable h
                                   Fourth column - Activation variable H
                                   Fifth column - Gate variable n]
        |  :param int_name: Integration method (please refer to 
            http://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html)
        |  :param *args: If int_method is one of: vode, zvode, or lsoda, the user
                         can choose one of Adams or BDF methods
        """
        t0 = 0.0         # Initial time (ms)
        dt = self.dt     # Integration time step (ms)

        # Some temporal conversions
        t_sim = int(tf / dt)        # Convert tf (ms) to simulation time
        P0_sim = int(self.P0 / dt)  # Convert period (ms) to simulation period
        p_sim = int(self.P / dt)    # Convert duration (ms) to sim. duration

        # Define the input current
        if self.inp_type == "const":
            Iapp = np.ones((t_sim, )) * self.amplitude
        elif self.inp_type == "pulse":
            Iapp = np.zeros((t_sim, ))
            dur = int(self.P / dt)
            Iapp[:dur] = self.pulse_amplitude
            Iapp[dur:] = self.amplitude
        elif self.inp_type == "periodic":
            Iapp = self.I_ext(P0_sim, p_sim, t_sim)
            if t_sim > Iapp.shape[0]:
                raise ValueError('Simulation time exceeds input signal dimensions!')
        else:
            raise ValueError('Bad input signal argument!')

        V0 = self.p['V0']
        h0 = self.h_inf(V0, self.p['theta_h'], self.p['k_h'])
        H0 = self.H_inf(V0)
        n0 = self.n_inf(V0, self.p['sigma_K'])
        x0 = [V0, h0, H0, n0]

        # Set parameters and initialize the solver
        if len(args) == 0:
            r = ode(self.rhs).set_integrator(int_name)
        else:
            r = ode(self.rhs).set_integrator(int_name, method=args[0])
        r.set_initial_value(x0, t0).set_f_params(Iapp[0])

        v, t = [], []
        i = 0
        while r.successful() and r.t < (tf-1):
            r.set_f_params(Iapp[i])
            v.append(r.integrate(r.t+dt))
            t.append(r.t)
            i += 1

        # Store the neuron states and time in the same file if the 
        # flag store is true.
        v = np.array(v)
        t = np.array(t)
        v = np.hstack([t.reshape(t.shape[0], 1), v])

        if self.IsStore is True:
            np.save("../data/"+vname, v)
        if self.IsStoreStim is True:
            np.save("../data/"+vname+"stim", Iapp)
        return v
