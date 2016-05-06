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


def loadParameters(fname):
    """ Load all the necessary paremeters from a file.

        |  :param fname: File name containing all parameters
        |  :return: A dictionaty with all the parameters
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
            'phi_K': config.getfloat('Potasium', 'phi_K'),
            'sigma_K': config.getfloat('Potasium', 'sigma_K'),
            'gK': config.getfloat('Potasium', 'gK'),
            'VK': config.getfloat('Potasium', 'VK'),
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
    def __init__(self, fname, frequency, duration, amplitude,
                 chunks):
        """ Constructor method

        |  :param fname: File name containing all parameters
        |  :param frequency: Stimulation frequency (in Hz)
        |  :param duration: Duration (ON) of stimulation (in ms)
        |  :param amplitude: Amplitude of injected current (in μA/cm²)
        |  :param chunks: Number of repeated periods of pulses
        |  :return:
        """
        self.p = loadParameters(fname)      # Load the parameters

        self.freq = frequency
        self.P0 = (1.0/self.freq)*1000      # Period in ms
        self.P = duration                   # Duration in ms
        self.Amp = amplitude                # Amplitude in μA/cm²

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
                (v+29.7-sigma)) - 1.))

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

    def I_ext(self, period, duration, amplitude, time, numChunks):
        tmp = np.zeros((period, ))
        tmp[:duration] = amplitude
        repeats = int(numChunks + (time - period) / period)
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

    def run(self, tf, vname, tname):
        t0 = 0.0         # Initial time (ms)
        dt = 0.05        # Integration time step (ms)

        # Some temporal conversions
        t_sim = int(tf / dt)        # Convert tf (ms) to simulation time
        P0_sim = int(self.P0 / dt)  # Convert period (ms) to simulation period
        p_sim = int(self.P / dt)    # Convert duration (ms) to sim. duration

        # Define the input current
        Iapp = self.I_ext(P0_sim, p_sim, self.Amp, t_sim, 5)

        V0 = self.p['V0']
        h0 = self.h_inf(V0, self.p['theta_h'], self.p['k_h'])
        H0 = self.H_inf(V0)
        n0 = self.n_inf(V0, self.p['sigma_K'])
        x0 = [V0, h0, H0, n0]

        # Set parameters and initialize the solver
        r = ode(self.rhs).set_integrator('vode', method='BDF')
        r.set_initial_value(x0, t0).set_f_params(Iapp[0])

        v, t = [], []
        i = 0
        while r.successful() and r.t < tf:
            r.set_f_params(Iapp[i])
            v.append(r.integrate(r.t+dt))
            t.append(r.t)
            i += 1

        v = np.array(v)
        t = np.array(t)

        np.save("../data/"+vname, v)
        np.save("../data/"+tname, t)


if __name__ == '__main__':
    base = "params/"
    
    # Build data for Figure 1
    freq = [10, 10, 5, 8, 5]        # Pulse frequency in Hz
    dur = [10, 40, 40, 75, 120]     # Pulse duration in ms
    Iext = -1.0                     # Pulse current in mA/cm^2

    name = "Fig1"
    for i in range(len(freq)):
        simulation = WangModel(base+"params_figure1.cfg",
                               freq[i],
                               dur[i],
                               Iext,
                               100)
        simulation.run(6000, name+"_V"+str(i), name+"_T"+str(i))

    # Build data for Figure 2
    Iext = [3.0, 0.0, -0.5, -0.55, -0.6, -0.8, -1.3, -2.1]
    name = "Fig2"
    for i in Iext:
        simulation = WangModel(base+"params_figure2.cfg", 1, 1000, i, 1)
        simulation.run(2000, name+"_V"+str(i), name+"_T"+str(i))
    
    # Build data for Figure 3
    simulation = WangModel(base+"params_figure3.cfg", 1, 1000, -0.91, 2)
    simulation.run(6000, "Fig3_V", "Fig3_T")
