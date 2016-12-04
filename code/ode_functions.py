"""
   Copyright 2016 Aaron R. Shifman

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

import numpy as np
from scipy.integrate import odeint


def m_inf(v):
    """
    Sodium gating dynamics
    :param v: Membrane potential (mV)
    :return: m_inf(v)
    """
    return 1. / (1 + np.exp(-(v + 36.) / 8.5))


def tau_m(v):
    """
    Sodium gating dynamics
    :param v: Membrane potential (mV) (unused)
    :return: tau_m(v)
    """
    return 0.1


def h_inf(v):
    """
    Sodium gating dynamics
    :param v: Membrane potential (mV)
    :return: h_inf(v)
    """
    return 1. / (1 + np.exp((v + 44.1) / 7.))


def tau_h(v):
    """
    Sodium gating dynamics
    :param v: Membrane potential (mV)
    :return: tau_h(v)
    """
    return 3.5 / (np.exp((v + 35) / 4.) + np.exp(-(v + 35) / 25.)) + 1


def m_nap_inf(v):
    """
    Persistent sodium dynamics
    :param v: Membrane potential (mV)
    :return: m_nap_inf(v)
    """
    return 1. / (1 + np.exp(-(v + 47.1) / 4.1))


def tau_m_nap(v):
    """
    Persistent sodium dynamics
    :param v: Membrane potential (mV) (unused)
    :return: tau_m_nap(v)
    """
    return 0.1


def h_nap_inf(v):
    """
    Persistent sodium dynamics
    :param v: Membrane potential (mV)
    :return: h_nap_inf(v)
    """
    return 1. / (1 + np.exp((v + 65.) / 5.))


def tau_h_nap(v):
    """
    Persistent sodium dynamics
    :param v: Membrane potential (mV) (unused)
    :return: tau_h_nap(v)
    """
    return 150.


def n_inf(v):
    """
    Delayed rectifier potassium dynamics
    :param v: Membrane potential (mV)
    :return: n_inf(v)
    """
    return 1. / (1 + np.exp(-(v + 30) / 25.))


def tau_n(v):
    """
    Delayed rectifier potassium dynamics
    :param v: Membrane potential (mV)
    :return: tau_n(v)
    """
    return 2.5 / (np.exp((v + 30.) / 40.) + np.exp(-(v + 30.) / 50.)) + 0.01


def m_t_inf(v):
    """
    LVA calcium dynamics
    :param v: Membrane potential (mV)
    :return: m_t_inf(v)
    """
    return 1. / (1. + np.exp(-(v + 38.) / 5.))


def tau_m_t(v):
    """
    LVA calcium dynamics
    :param v: Membrane potential (mV)
    :return: tau_m_inf(v)
    """
    return 5. / (np.exp((v + 28) / 25.) + np.exp(-(v + 28) / 70.)) + 2


def h_t_inf(v):
    """
    LVA calcium dynamics
    :param v: Membrane potential (mV)
    :return: h_t_inf(v)
    """
    return 1. / (1. + np.exp((v + 70.1) / 7.))


def tau_h_t(v):
    """
    LVA calcium dynamics
    :param v: Membrane potential (mV)
    :return: tau_h_t(v)
    """
    return 20. / (np.exp((v + 70.) / 65.) + np.exp(-(v + 70.) / 65.)) + 1


def m_n_inf(v):
    """
    HVA calcium dynamics
    :param v: Membrane potential (mV)
    :return: m_n_inf(v)
    """
    return 1. / (1. + np.exp(-(v + 30.) / 6.))


def tau_m_n(v):
    """
    HVA calcium dynamics
    :param v: Membrane potential (mV) (unused)
    :return: tau_m_n(v)
    """
    return 5.


def h_n_inf(v):
    """
    HVA calcium dynamics
    :param v: Membrane potential (mV)
    :return: h_n_inf(v)
    """
    return 1. / (1. + np.exp((v + 70.) / 3.))


def tau_h_n(v):
    """
    HVA calcium dynamics
    :param v: Membrane potential (mV)
    :return: tau_h_n(v)
    """
    return 25.


def m_p_inf(v):
    """
    HVA calcium dynamics
    :param v: Membrane potential (mV)
    :return: m_p_inf(v)
    """
    return 1. / (1 + np.exp(-(v + 17.) / 3.))


def tau_m_p(v):
    """
    HVA calcium dynamics
    :param v: Membrane potential (mV) (unused)
    :return: tau_m_p(v)
    """
    return 10.


def z_sk_inf(ca):
    """
    Calcium dependent potassium dynamics
    :param ca: [Ca] \mu M
    :return: z_sk_inf(ca)
    """
    return 1. / (1. + ((0.003 / ca) ** 2))


def tau_z_sk(ca):
    """
    Calcium dependent potassium dynamics
    :param ca: [Ca] \mu M
    :return: tau_z_sk(ca)
    """
    return 1.


def m_a_inf(v):
    """
    Fast transient potassium dynamics
    :param v: Membrane potential (mV)
    :return: m_a_inf(v)
    """
    return 1. / (1 + np.exp(-(v + 27) / 16.))


def tau_m_a(v):
    """
    Fast transient potassium dynamics
    :param v: Membrane potential (mV)
    :return: tau_m_a(v)
    """
    return 1. / (np.exp((v + 40) / 5.) + np.exp(-(v + 74) / 7.5)) + 0.37


def h_a_inf(v):
    """
    Fast transient potassium dynamics
    :param v: Membrane potential (mV)
    :return: h_a_inf(v)
    """
    return 1. / (1 + np.exp((v + 80) / 11.))


def tau_h_a(v):
    """
    Fast transient potassium dynamics
    :param v: Membrane potential (mV) (unused)
    :return: tau_h_a(v)
    """
    return 20.


def m_h_inf(v):
    """
    Hyperpolarization activated dynamics
    :param v: Membrane potential (mV)
    :return: m_h_inf(v)
    """
    return 1. / (1. + np.exp((v + 79.8) / 5.3))


def tau_m_h(v):
    """
    Hyperpolarization activated dynamics
    :param v: Membrane potential (mV)
    :return: tau_m_h(v)
    """
    return 475. / (np.exp((v + 70.) / 11.) + np.exp(-(v + 70.) / 11.)) + 50.


def solver(t1, g_t_bar=0.1, g_sk_bar=0.3, g_n_bar=0.05, t_start=10, duration=1,
           i_bias_on=1, g_h_bar=0.005, y_hold=None, ca_type=None):
    """
    Centralized function to conglomerate all parameters and settings to solve the simulation
    Simulation is solved with odeint from scipy.integrate

    :param t1: Experiment end time (ms)
    :param t1: End time of simulation
    :param g_t_bar: Maximal conductance for T-type calcium current. Defaults to 0.1
    :param g_sk_bar: Maximal conductance for SK-type current. Defaults to 0.3
    :param g_n_bar: Maximal conductance for N current. Defaults to 0.05
    :param g_h_bar: Maximal conductance for H current. Defaults to 0.005
    :param t_start: Current step/pulse start time
    :param duration: Length of current step/pulse
    :param i_bias_on: Current amplitude
    :param y_hold: Optional, if there is a different set of initial conditions
    :param ca_type: Optional, if there is a condition where ICa different: ca_type=1, Ica=IP+IN. ca_type=2, ICa=IP + IT
    :return: t,Y
    """

    params = [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar]
    params = np.concatenate((params, [ca_type]))  # Add ca_type to params

    if y_hold is not None:
        y0 = y_hold
    else:
        y0 = np.array([-7.16300326e+01, 1.48943267e-02, 9.80788764e-01,
                       2.51507410e-03, 7.90179298e-01, 1.59065019e-01,
                       1.19787622e-03, 5.54427500e-01, 1.23450475e-08,
                       9.68939709e-04, 6.32590017e-01, 2.03398207e-03,
                       5.79009277e-02, 3.18449568e-01, 1.76316847e-01,
                       1.35436879e-04])  # steady state conditions for "typical" experiment
    dt = 0.05
    t = np.linspace(0, t1, int(t1 / dt))  # create vector of times

    y = odeint(dydt, y0, t, args=(params,), hmax=0.05, atol=1e-10, rtol=1e-10)
    return t, y


def dydt(y, t, args):
    """
    Function to calculate the RHS of the set of differential equations
    :param y: [v, m, h, m_nap, h_nap, n, m_t, h_t, m_p, m_n, h_n, z_sk, m_a, h_a, m_h, ca]
    :param t: Current time point
    :param args: Model parameters [t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar, ICaType]
    :return: Returns dydt in the same order as y
    """

    """
    Static model parameters (constant across all simulations)
    """
    g_na_bar = 0.7
    g_nap_bar = 0.05
    g_k_bar = 1.3
    g_p_bar = 0.05
    g_leak = 0.005
    g_a_bar = 1.0

    e_na = 60
    e_k = -80
    e_leak = -50
    e_ca = 40
    e_h = -38.8

    k1 = -0.0005
    k2 = 0.04

    cm = 0.040

    v, m, h, m_nap, h_nap, n, m_t, h_t, m_p, m_n, h_n, z_sk, m_a, h_a, m_h, ca = y  # extract from vector
    t_start, g_t_bar, duration, i_bias_on, g_sk_bar, g_n_bar, g_h_bar, ca_type = args  # extract parameters

    """
    Calculate all ionic (and other) currents
    """
    i_stim = (t >= t_start) * i_bias_on * (t < (t_start + duration))  # create current pulse (current is i_bias_on or 0)

    i_na = g_na_bar * (m ** 3) * h * (v - e_na)  # Na current
    i_nap = g_nap_bar * m_nap * h_nap * (v - e_na)  # NaP current

    i_k = g_k_bar * (n ** 4) * (v - e_k)  # K current

    i_leak = g_leak * (v - e_leak)  # leak current

    i_t = g_t_bar * m_t * h_t * (v - e_ca)  # T Calcium current
    i_n = g_n_bar * m_n * h_n * (v - e_ca)  # N Calcium current
    i_p = g_p_bar * m_p * (v - e_ca)  # P Calcium current

    i_sk = g_sk_bar * (z_sk ** 2) * (v - e_k)  # SK Current
    i_a = g_a_bar * m_a * h_a * (v - e_k)  # A current
    i_h = g_h_bar * m_h * (v - e_h)  # H current

    """
    Calculate derivatives f'(t) for all dynamic variables
    """
    dm = (m_inf(v) - m) / tau_m(v)
    dh = (h_inf(v) - h) / tau_h(v)
    dn = (n_inf(v) - n) / tau_n(v)
    dh_nap = (h_nap_inf(v) - h_nap) / tau_h_nap(v)
    dm_nap = (m_nap_inf(v) - m_nap) / tau_m_nap(v)
    dm_t = (m_t_inf(v) - m_t) / tau_m_t(v)
    dh_t = (h_t_inf(v) - h_t) / tau_h_t(v)
    dm_n = (m_n_inf(v) - m_n) / tau_m_n(v)
    dh_n = (h_n_inf(v) - h_n) / tau_h_n(v)
    dm_p = (m_p_inf(v) - m_p) / tau_m_p(v)
    dz_sk = (z_sk_inf(ca) - z_sk) / tau_z_sk(ca)
    dm_a = (m_a_inf(v) - m_a) / tau_m_a(v)
    dh_a = (h_a_inf(v) - h_a) / tau_h_a(v)
    dm_h = (m_h_inf(v) - m_h) / tau_m_h(v)

    i_ca = i_t + i_n + i_p  # Explicit definition of i_ca as sum of calcium currents (T, N, P)
    if ca_type == 1:
        i_ca = i_p + i_n
    elif ca_type == 2:
        i_ca = i_p + i_t

    d_ca = k1 * i_ca - k2 * ca

    iion = i_na + i_nap + i_k + i_leak + i_t + i_n + i_p + i_sk + i_a + i_h  # Explicit expression of ionic currents
    d_v = (1 / cm) * (-iion + i_stim)

    return np.array([d_v, dm, dh, dm_nap, dh_nap, dn, dm_t, dh_t, dm_p, dm_n, dh_n, dz_sk, dm_a, dh_a, dm_h, d_ca])
