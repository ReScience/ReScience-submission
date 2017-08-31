# -*- coding: utf-8 -*-
"""
Code for Diem AK (2017) Chemical Signalling in the Neurovascular Unit, The
ReScience Journal

https://github.com/akdiem/ReScience-submission/tree/Diem-2017

This is a reimplementation of Witthoft A, Karniadakis GE (2012) A bidirectional
model for communication in the neurovascular unit. Journal of Theoretical
Biology 311: 80-93.

This file contains the ODE describing chemical transports and blood vessel
mechanics in the neurovascular unit.
"""

import numpy as np
from scipy.integrate import ode
from scipy.interpolate import interp1d
import matplotlib.pylab as plt


def synapse(t, potassium_s, Jrho_IN, JSigKkNa, KKoa, **kwargs):
    """Implements the equations for chemical transport in the synaptic space
    
    Input variables are the K+ (uM/s) concentration and ratio of bound/unbound
    Glu receptors. They are realised as step functions with a steep linear increase
    to a maximum sustained value and steep linear decrease back to 0 once the
    stimulus is removed.

    Parameters
    --------------
    t : float
        Time.
    potassium_s : float
        Value of synaptic K+ at time t.
    Jrho_IN : array
        Array containing input values.

    Returns
    ---------
    float
        Value of ODE for synaptic K+ at time t.
    float
        Value of eq 2 at time t.
    """
    JSigK = JSigKkNa * potassium_s/(potassium_s + KKoa) # eq 2
    JKss = interp1d(Jrho_IN[:,0], Jrho_IN[:,1], bounds_error=False, fill_value=0)
    potassium_s_dt = JKss(t) - JSigK # eq 1
    return potassium_s_dt, JSigK


def astrocyte(t, ip3, calcium_a, h, ss, Vk, calcium_p, x, eet, nbk, Jrho_IN,
              x_rel, JSigK, delta, KG, rh, kdeg, Jmax, Ki, Kact, calcium_er,
              Vmax, Kp, Pl, gtrpv, vtrpv, Castr, gamma, beta, kon, Kinh, 
              tautrpv, uM, gammacai, gammacae, eps12, kappa, v1trpv, v2trpv,
              Veet, calcium_a_min, keet, v5bk, Ca3bk, Ca4bk, v6bk, psibk,
              v4bk, eetshift, gbk, vbk, gleak, vleak, trpv=True, **kwargs):
    """Implements the equations for chemical transport in the astrocytic space

    Parameters
    --------------
    t : float
        Time.
    ip3 : float
        Value of IP3 at time t.
    calcium_a : float
        Value of astrocyte Ca2+ at time t.
    h : float
        Value of gating variable h at time t.
    ss : float
        Value of gating variable ss at time t.
    Vk : float
        Value of astrocyte membrane potential at time t.
    calcium_p : float
        Value of perivascular Ca2+ at time t.
    x : float
        Value of vessel circumference at time t.
    eet : float
        Value of EET at time t.
    nbk : float
        Value of gating variable nbk at time t.
    Jrho_IN : array
        Array containing input values.
    x_rel : float
        Relaxed value of blood vessel.

    Returns
    ---------
    float
        Value of ODE for IP3 at time t.
    float
        Value of ODE for astrocyte Ca2+ at time t.
    float
        Value of ODE for gating variable h at time t.
    float
        Value of ODE for gating variable ss at time t.
    float
        Value of ODE for EET at time t.
    float
        Value of ODE for gating variable nbk at time t.
    float
        Value of ODE for astrocyte membrane potential at time t.
    float
        Value of ODE for IP3 at time t.
    float
        Value of Ibk (eq. 15) at time t.
    float
        Value of Jtrpv (eq 10) at time t.
    """
    rhos = interp1d(Jrho_IN[:,0], Jrho_IN[:,2], bounds_error=False,
                    fill_value=0)
    G = (rhos(t) + delta)/(KG + rhos(t) + delta) # eq 3
    ip3_dt = rh*G - kdeg*ip3 # eq 4
    Jip3 = Jmax * (ip3/(ip3+Ki) * calcium_a/(calcium_a+Kact) * h)**3 *\
        (1 - calcium_a/calcium_er) # eq 6
    Jpump = Vmax * calcium_a**2 / (calcium_a**2 + Kp**2) # eq 8
    Jleak = Pl * (1 - calcium_a/calcium_er) # eq 9
    Itrpv = gtrpv * ss * (Vk - vtrpv) # eq 10
    Jtrpv = Itrpv/(Castr*gamma) # eq 10
    if trpv:
        calcium_a_dt = beta * (Jip3 - Jpump + Jleak + Jtrpv) # eq 5
    else:
        calcium_a_dt = beta * (Jip3 - Jpump + Jleak)
    h_dt = kon * (Kinh - (calcium_a + Kinh) * h) # eq 7
    tauca = tautrpv / (calcium_p/uM) # eq 11
    eps = (x - x_rel)/x_rel # eq 11
    Hca = calcium_a/gammacai + calcium_p/gammacae # eq 13
    sinf = (1/(1 + np.exp(-(eps-eps12)/kappa))) * (1/(1+Hca) *\
        (Hca + np.tanh((Vk - v1trpv)/v2trpv))) # eq 12
    ss_dt = 1/tauca * (sinf - ss) # eq 11
    eet_dt = Veet * (calcium_a - calcium_a_min) - keet*eet # eq 14
    v3bk = -(v5bk/2) * np.tanh((calcium_a-Ca3bk)/Ca4bk) + v6bk # eq 19
    phibk = psibk * np.cosh((Vk-v3bk)/(2*v4bk)) # eq 17
    ninfbk = 0.5 * (1 + np.tanh((Vk + eetshift*eet - v3bk)/v4bk)) # eq 18
    nbk_dt = phibk * (ninfbk - nbk) # eq 16
    Ibk = gbk * nbk * (Vk - vbk) # eq 15
    Ileak = gleak * (Vk - vleak) # eq 21
    Isigk = -JSigK * Castr * gamma # eq 20
    if trpv:
        Vk_dt = 1/Castr * (-Isigk - Ibk - Ileak - Itrpv) # eq 20
    else:
        Vk_dt = 1/Castr * (-Isigk - Ibk - Ileak)
    return ip3_dt, calcium_a_dt, h_dt, ss_dt, eet_dt, nbk_dt, Vk_dt, Ibk, Jtrpv


def perivascular_space(potassium_p, k, Vm, calcium_p, Ibk, Jtrpv, Castr,
                       gamma, gkir0, mM, vkir1, vkir2, Csmc, VRpa, VRps,
                       Rdecay, potassium_p_min, dp, mmHg, mV, v2, gca, vca,
                       Cadecay, calcium_p_min, trpv=True, **kwargs):
    """Implements the equations for chemical transport in the perivascular
    space

    Parameters
    --------------
    potassium_p : float
        Value of perivascular K+ at current time step.
    k : float
        Value of gating variable k at current time step.
    Vm : float
        Value of SMC membrane potential at current time step.
    calcium_p : float
        Value of perivascular Ca2+ at current time step.
    Ibk : float
        Value of Ibk (eq 15) at current time step.
    Jtrpv : float
        Value of Jtrpv (eq 10) at current time step.

    Returns
    ---------
    float
        Value of ODE for perivascular K+ at time t.
    float
        Value of ODE for perivascular Ca2+ at time t.
    float
        Value of vkir (eq 26) at current time step.
    float
        Value of Ikir (eq 24) at current time step.
    float
        Value of Ica (eq 36) at current time step.
    """
    Jbk = Ibk/(Castr*gamma) # eq 22
    gkir = gkir0 * np.sqrt(potassium_p/mM) # eq 25
    vkir = vkir1 * np.log10(potassium_p/mM) - vkir2 # eq 26
    Ikir = gkir * k * (Vm - vkir) # eq 24
    Jkir = Ikir/(Csmc*gamma) # 22
    potassium_p_dt = Jbk/VRpa + Jkir/VRps - Rdecay * (potassium_p -\
        potassium_p_min) # eq 22
    v1 = (-17.4-12*(dp/mmHg)/200)*mV # eq 38
    minf = 0.5 * (1 + np.tanh((Vm-v1)/v2)) # eq 37
    Ica = gca * minf * (Vm - vca) # eq 36
    Jca = -Ica/(Csmc*gamma) # not defined in manuscript
    if trpv:
        calcium_p_dt = -Jtrpv/VRpa - Jca/VRps - Cadecay * (calcium_p - calcium_p_min) # eq 23
    else:
        calcium_p_dt = -Jca/VRps - Cadecay * (calcium_p - calcium_p_min) # eq 23
    return potassium_p_dt, calcium_p_dt, vkir, Ikir, Ica


def ion_currents(Vm, k, calcium_smc, n, vkir, Ica, Ikir, alphakir, av1, av2,
                 betakir, bv2, bv1, mV, gl, vl, gk, vk, Csmc, v5, Ca3, Ca4, v6,
                 phin, v4, **kwargs):
    """Implements the equations for ion currents in the neurovascular unit.

    Parameters
    --------------
    Vm : float
        Value of SMC membrane potential at current time step.
    k : float
        Value of gating variable k at current time step.
    calcium_smc : float
        Value of SMC Ca2+ at current time step.
    n : float
        Value of gating variable n at current time step.
    vkir : float
        Value of vkir (eq 26) at current time step.
    Ica : float
        Value of Ica (eq 36) at current time step.
    Ikir : float
        Value of Jkir (eq 24) at current time step.

    Returns
    ---------
    float
        Value of ODE for gating variable k at current time step.
    float
        Value of ODE for SMC membrane potential at current time step.
    float
        Value of ODE for gating variable n at current time step.
    """
    alphak = alphakir / (1 + np.exp((Vm - vkir + av1)/av2)) # eq 28
    betak = betakir * np.exp(bv2/mV * (Vm - vkir + bv1)/mV) # eq 29
    tauk = 1/(alphak+betak) # eq 27
    kinf = alphak/(alphak+betak) # eq 27
    k_dt = 1/tauk * (kinf - k) # eq 27
    Il = gl * (Vm - vl) # eq 30
    Ik = gk * n * (Vm - vk) # eq 31
    Vm_dt = (1/Csmc) * (-Il - Ik - Ica - Ikir) # eq 30
    v3 = -(v5/2) * np.tanh((calcium_smc-Ca3)/Ca4) + v6 # eq 35
    lamn = phin * np.cosh((Vm-v3)/(2*v4)) # eq 34
    ninf = 0.5 * (1 + np.tanh((Vm-v3)/v4)) # eq 33
    n_dt = lamn * (ninf - n) # eq 32
    return k_dt, Vm_dt, n_dt


def vessel_mechanics(t, calcium_smc, x, yy, omega, Ica, Kd, Bt, alpha, kca, dp,
                     Ax, um, x0, x3, x1, x2, x4, x5, x8, x6, x7, x9, we, Sx,
                     sigma0h, u2, u1, u3, wm, tau, Cam, q, kpsi, psim,
                     Caref, sigmay0h, y1, y2, y4, y0, y3, vref, ad, cd, bd, dd,
                     **kwargs):
    """Implements the equations for blood vessel mechanics in the neurovascular
    unit.

    Parameters
    --------------
    t : float
        Time.
    calcium_smc : float
        Value of SMC Ca2+ at time t.
    x : float
        Value of vessel circumference at time t.
    yy : float
        Value of length of the contractile component at time t.
    omega : float
        Value of myogenic stress at time t.
    Ica : float
        Value of Ica (eq 36) at time t.

    Returns
    ---------
    float
        Value of ODE for the vessel circumference at time t.
    float
        Value of ODE for SMC Ca2+ at time t.
    float
        Value of ODE for myogenic stress at time t.
    float
        Value of ODE for the length of the contractile component at time t.
    """
    # SMC calcium
    rho_smc = (Kd+calcium_smc)**2/((Kd+calcium_smc)**2 + Kd*Bt) # eq A.2
    calcium_smc_dt = -rho_smc * (alpha*Ica + kca*calcium_smc) # eq A.1

    # Vessel mechanics
    fdp = 0.5 * dp * (x/np.pi - Ax/x) * um # eq A.3
    xd = x/x0
    sigmax = x3*(1 + np.tanh((xd-x1)/x2)) + x4*(xd-x5) - x8*(x6/(xd-x7))**2 -\
                x9 # eq A.4
    fx = we*Sx*sigmax*sigma0h # eq A.11
    yd = yy/x0
    ud = xd-yd
    sigmau = u2 * np.exp(u1*ud) - u3 # eq A.5
    fu = wm*Sx*sigmau*sigma0h # eq A.11
    x_dt = 1/tau * (fdp - fx - fu) # eq A.12
    psi = calcium_smc**q/(Cam**q+calcium_smc**q) # eq A.7
    omega_dt = kpsi * (psi/(psim+psi) - omega) # eq A.8
    psiref = Caref**q/(Cam**q+Caref**q)
    omega_ref = psiref/(psim + psiref) # eq A.9
    sigmay0 = sigmay0h * omega/omega_ref # eq A.9
    sy = (y1/(yd+y2))**y4
    # eq A.10
    sigmay = sigmay0/sigma0h * (np.exp(-(yd-y0)**2/(2*sy**2)) - y3)/(1-y3)
    cond = sigmau/sigmay
    if cond < 1:
        ycond = -vref * (psi/psiref) * ad * (1-cond)/(ad+cond)
    else:
        ycond = cd*(np.exp(bd * (cond-dd)) - np.exp(bd * (1-dd)))
    yy_dt = x0*ycond
    
    return x_dt, calcium_smc_dt, omega_dt, yy_dt


def nvu(t, y, Jrho_IN, x_rel, units, param):
    """Implements the equations describing chemical transport and blood vessel
    mechanics in the neurovascular unit

    Parameters
    --------------
    t : float
        Time.
    y : array
        Array containing current solution values.
    Jrho_IN : array
        Array containing model input values.
    x_rel : float
        Value of relaxed vessel circumference at time t.
    units : dict
        Dictionary containing values for units.
    param : dict
        Dictionary containing values for parameters.

    Returns
    ---------
    array
        Array containing values of ODEs for all variables in the neurovascular
        unit at time t.
    """
    potassium_s = y[0]
    ip3 = y[1]
    calcium_a = y[2]
    h = y[3]
    ss = y[4]
    eet = y[5]
    nbk = y[6]
    Vk = y[7]
    potassium_p = y[8]
    calcium_p = y[9]
    k = y[10]
    Vm = y[11]
    n = y[12]
    x = y[13]
    calcium_smc = y[14]
    omega = y[15]
    yy = y[16]
    
    # Synaptic space
    potassium_s_dt, JSigK = synapse(t, potassium_s, Jrho_IN, **param)    
    
    # Astrocytic space
    ip3_dt, calcium_a_dt, h_dt, ss_dt, eet_dt, nbk_dt, Vk_dt, Ibk, Jtrpv =\
        astrocyte(t, ip3, calcium_a, h, ss, Vk, calcium_p, x, eet, nbk,
                  Jrho_IN, x_rel, JSigK, **units, **param)    
    
    # Perivascular space
    potassium_p_dt, calcium_p_dt, vkir, Ikir, Ica = perivascular_space(
            potassium_p, k, Vm, calcium_p, Ibk, Jtrpv, **units, **param)   
    
    # Ion currents
    k_dt, Vm_dt, n_dt = ion_currents(Vm, k, calcium_smc, n, vkir, Ica, Ikir,
                                     **units, **param)
    
    # Vessel mechanics
    x_dt, calcium_smc_dt, omega_dt, yy_dt = vessel_mechanics(t, calcium_smc, x,
                                            yy, omega, Ica, **units, **param)
    
    return [potassium_s_dt, ip3_dt, calcium_a_dt, h_dt, ss_dt, eet_dt, nbk_dt,
            Vk_dt, potassium_p_dt, calcium_p_dt, k_dt, Vm_dt, n_dt, x_dt,
            calcium_smc_dt, omega_dt, yy_dt]


def nvu_Vk(t, y, Jrho_IN, x_rel, Vk_f, units, param):
    """Same as nvu(), but Vk is fixed to values provided by Vk_f

    Parameters
    --------------
    t : float
        Time.
    y : array
        Array containing current solution values.
    Jrho_IN : array
        Array containing model input values.
    x_rel : float
        Value of relaxed vessel circumference at time t.
    Vk_f : function
        Function returning values for astrocyte membrane potential at time t.
    units : dict
        Dictionary containing values for units.
    param : dict
        Dictionary containing values for parameters.

    Returns
    ---------
    array
        Array containing values of ODEs for all variables in the neurovascular
        unit at time t.
    """
    potassium_s = y[0]
    ip3 = y[1]
    calcium_a = y[2]
    h = y[3]
    ss = y[4]
    eet = y[5]
    nbk = y[6]
    Vk = Vk_f(t)
    potassium_p = y[8]
    calcium_p = y[9]
    k = y[10]
    Vm = y[11]
    n = y[12]
    x = y[13]
    calcium_smc = y[14]
    omega = y[15]
    yy = y[16]
    
    # Synaptic space
    potassium_s_dt, JSigK = synapse(t, potassium_s, Jrho_IN, **param)    
    
    # Astrocytic space
    ip3_dt, calcium_a_dt, h_dt, ss_dt, eet_dt, nbk_dt, Vk_dt, Ibk, Jtrpv =\
        astrocyte(t, ip3, calcium_a, h, ss, Vk, calcium_p, x, eet, nbk,
                  Jrho_IN, x_rel, JSigK, **units, **param)    
    
    # Perivascular space
    potassium_p_dt, calcium_p_dt, vkir, Ikir, Ica = perivascular_space(
            potassium_p, k, Vm, calcium_p, Ibk, Jtrpv, **units, **param)   
    
    # Ion currents
    k_dt, Vm_dt, n_dt = ion_currents(Vm, k, calcium_smc, n, vkir, Ica, Ikir,
                                     **units, **param)
    
    # Vessel mechanics
    x_dt, calcium_smc_dt, omega_dt, yy_dt = vessel_mechanics(t, calcium_smc, x,
                                            yy, omega, Ica, **units, **param)
    
    return [potassium_s_dt, ip3_dt, calcium_a_dt, h_dt, ss_dt, eet_dt, nbk_dt,
            Vk_dt, potassium_p_dt, calcium_p_dt, k_dt, Vm_dt, n_dt, x_dt,
            calcium_smc_dt, omega_dt, yy_dt]


def nvu_k(t, y, Jrho_IN, x_rel, k_f, units, param):
    """Same as nvu(), but k is fixed to values provided by k_f

    Parameters
    --------------
    t : float
        Time.
    y : array
        Array containing current solution values.
    Jrho_IN : array
        Array containing model input values.
    x_rel : float
        Value of relaxed vessel circumference at time t.
    k_f : function
        Function returning values for gating variable k at time t.
    units : dict
        Dictionary containing values for units.
    param : dict
        Dictionary containing values for parameters.

    Returns
    ---------
    array
        Array containing values of ODEs for all variables in the neurovascular
        unit at time t.
    """
    potassium_s = y[0]
    ip3 = y[1]
    calcium_a = y[2]
    h = y[3]
    ss = y[4]
    eet = y[5]
    nbk = y[6]
    Vk = y[7]
    potassium_p = y[8]
    calcium_p = y[9]
    k = k_f(t)
    Vm = y[11]
    n = y[12]
    x = y[13]
    calcium_smc = y[14]
    omega = y[15]
    yy = y[16]
    
    # Synaptic space
    potassium_s_dt, JSigK = synapse(t, potassium_s, Jrho_IN, **param)    
    
    # Astrocytic space
    ip3_dt, calcium_a_dt, h_dt, ss_dt, eet_dt, nbk_dt, Vk_dt, Ibk, Jtrpv =\
        astrocyte(t, ip3, calcium_a, h, ss, Vk, calcium_p, x, eet, nbk,
                  Jrho_IN, x_rel, JSigK, **units, **param)    
    
    # Perivascular space
    potassium_p_dt, calcium_p_dt, vkir, Ikir, Ica = perivascular_space(
            potassium_p, k, Vm, calcium_p, Ibk, Jtrpv, **units, **param)   
    
    # Ion currents
    k_dt, Vm_dt, n_dt = ion_currents(Vm, k, calcium_smc, n, vkir, Ica, Ikir,
                                     **units, **param)
    
    # Vessel mechanics
    x_dt, calcium_smc_dt, omega_dt, yy_dt = vessel_mechanics(t, calcium_smc, x,
                                            yy, omega, Ica, **units, **param)
    
    return [potassium_s_dt, ip3_dt, calcium_a_dt, h_dt, ss_dt, eet_dt, nbk_dt,
            Vk_dt, potassium_p_dt, calcium_p_dt, k_dt, Vm_dt, n_dt, x_dt,
            calcium_smc_dt, omega_dt, yy_dt]



def nvu_trpv(t, y, Jrho_IN, x_rel, r_f, units, param):
    """Same as nvu(), but vessel radius is fixed to values provided by r_f

    Parameters
    --------------
    t : float
        Time.
    y : array
        Array containing current solution values.
    Jrho_IN : array
        Array containing model input values.
    x_rel : float
        Value of relaxed vessel circumference at time t.
    r_f : function
        Function returning values for vessel radius at time t.
    units : dict
        Dictionary containing values for units.
    param : dict
        Dictionary containing values for parameters.

    Returns
    ---------
    array
        Array containing values of ODEs for all variables in the neurovascular
        unit at time t.
    """
    potassium_s = y[0]
    ip3 = y[1]
    calcium_a = y[2]
    h = y[3]
    ss = y[4]
    eet = y[5]
    nbk = y[6]
    Vk = y[7]
    potassium_p = y[8]
    calcium_p = y[9]
    k = y[10]
    Vm = y[11]
    n = y[12]
    x = 2*np.pi*r_f(t)
    calcium_smc = y[14]
    omega = y[15]
    yy = y[16]
    
    # Synaptic space
    potassium_s_dt, JSigK = synapse(t, potassium_s, Jrho_IN, **param)    
    
    # Astrocytic space
    ip3_dt, calcium_a_dt, h_dt, ss_dt, eet_dt, nbk_dt, Vk_dt, Ibk, Jtrpv =\
        astrocyte(t, ip3, calcium_a, h, ss, Vk, calcium_p, x, eet, nbk,
                  Jrho_IN, x_rel, JSigK, **units, **param)    
    
    # Perivascular space
    potassium_p_dt, calcium_p_dt, vkir, Ikir, Ica = perivascular_space(
            potassium_p, k, Vm, calcium_p, Ibk, Jtrpv, **units, **param)   
    
    # Ion currents
    k_dt, Vm_dt, n_dt = ion_currents(Vm, k, calcium_smc, n, vkir, Ica, Ikir,
                                     **units, **param)
    
    # Vessel mechanics
    x_dt, calcium_smc_dt, omega_dt, yy_dt = vessel_mechanics(t, calcium_smc, x,
                                            yy, omega, Ica, **units, **param)
    
    return [potassium_s_dt, ip3_dt, calcium_a_dt, h_dt, ss_dt, eet_dt, nbk_dt,
            Vk_dt, potassium_p_dt, calcium_p_dt, k_dt, Vm_dt, n_dt, x_dt,
            calcium_smc_dt, omega_dt, yy_dt]
    
    
def nvu_notrpv(t, y, Jrho_IN, x_rel, units, param):
    """Same as nvu(), but trpv channels are deleted

    Parameters
    --------------
    t : float
        Time.
    y : array
        Array containing current solution values.
    Jrho_IN : array
        Array containing model input values.
    x_rel : float
        Value of relaxed vessel circumference at time t.
    r_f : function
        Function returning values for vessel radius at time t.
    units : dict
        Dictionary containing values for units.
    param : dict
        Dictionary containing values for parameters.

    Returns
    ---------
    array
        Array containing values of ODEs for all variables in the neurovascular
        unit at time t.
    """
    potassium_s = y[0]
    ip3 = y[1]
    calcium_a = y[2]
    h = y[3]
    ss = y[4]
    eet = y[5]
    nbk = y[6]
    Vk = y[7]
    potassium_p = y[8]
    calcium_p = y[9]
    k = y[10]
    Vm = y[11]
    n = y[12]
    x = y[13]
    calcium_smc = y[14]
    omega = y[15]
    yy = y[16]
    
    # Synaptic space
    potassium_s_dt, JSigK = synapse(t, potassium_s, Jrho_IN, **param)    
    
    # Astrocytic space
    ip3_dt, calcium_a_dt, h_dt, ss_dt, eet_dt, nbk_dt, Vk_dt, Ibk, Jtrpv =\
        astrocyte(t, ip3, calcium_a, h, ss, Vk, calcium_p, x, eet, nbk,
                  Jrho_IN, x_rel, JSigK, trpv=False, **units, **param)    
    
    # Perivascular space
    potassium_p_dt, calcium_p_dt, vkir, Ikir, Ica = perivascular_space(
            potassium_p, k, Vm, calcium_p, Ibk, Jtrpv, trpv=False, **units,
            **param)   
    
    # Ion currents
    k_dt, Vm_dt, n_dt = ion_currents(Vm, k, calcium_smc, n, vkir, Ica, Ikir,
                                     **units, **param)
    
    # Vessel mechanics
    x_dt, calcium_smc_dt, omega_dt, yy_dt = vessel_mechanics(t, calcium_smc, x,
                                            yy, omega, Ica, **units, **param)
    
    return [potassium_s_dt, ip3_dt, calcium_a_dt, h_dt, ss_dt, eet_dt, nbk_dt,
            Vk_dt, potassium_p_dt, calcium_p_dt, k_dt, Vm_dt, n_dt, x_dt,
            calcium_smc_dt, omega_dt, yy_dt]


def run_simulation(time, y0, *args, atol=1e-7, rtol=1e-7, mode=''):
    """Sets up the simulation using ode() from scipy.integrate.

    Parameters
    --------------
    time : array
        Simulation time.
    y0 : array
        Array containing initial conditions for solution values.
    args :
        Optional parameters.
    atol : float
        Absolute tolerance value for ode().
    rtol : function
        Relative tolerance value for ode().
    mode : string
        Simulation type selector, can be empty, 'Vk', 'k' or 'trpv'.

    Returns
    ---------
    array
        Array containing full solution of the neurovascular unit model.
    """
    integrator = "lsoda"
    ode15s = ode(nvu)
    if mode == 'Vk':
        ode15s = ode(nvu_Vk)
    elif mode == 'k':
        ode15s = ode(nvu_k)
    elif mode == 'trpv':
        ode15s = ode(nvu_trpv)
    elif mode == 'notrpv':
        ode15s = ode(nvu_notrpv)
    ode15s.set_f_params(*args)
    ode15s.set_integrator(integrator, atol=atol, rtol=rtol)
    ode15s.set_initial_value(y0, t=time[0])
    nt = len(time)
    sol = np.zeros([nt, len(y0)])
    sol[0,:] = y0
    for i in range(1, nt):
        if ode15s.successful():
            sol[i,:] = ode15s.integrate(time[i])
    return sol


def plot_solution(t, sol, fig_dims, uM, mV, mM, um, fname='', **kwargs):
    """Plot solution.

    Parameters
    --------------
    t : array
        Time.
    sol : array
        Array containing model solution.
    fig_dims : tuple
        Figure dimensions.
    """
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.rcParams['legend.fontsize'] = 9
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.serif'] = ['Arial']
    
    f, axarr = plt.subplots(4, 2)
    f.set_size_inches(fig_dims[0], h=fig_dims[1])
    
    # left side
    axarr[0, 0].plot(t, sol[:,0]/uM, label="", lw=2)
    axarr[0, 0].set_ylabel("K+ syn (uM)")
    axarr[1, 0].plot(t, sol[:,1]/uM, label="", lw=2)
    axarr[1, 0].set_ylabel("IP3 (uM)")
    axarr[2, 0].plot(t, sol[:,2]/uM, label="", lw=2)
    axarr[2, 0].set_ylabel("Ca2+ ast (uM)")
    axarr[3, 0].plot(t, sol[:,5]/uM, label="", lw=2)
    axarr[3, 0].set_ylabel("EET (uM)")

    # right side
    axarr[0, 1].plot(t, sol[:,7]/mV, label="", lw=2)
    axarr[0, 1].set_ylabel("Vk (mV)")
    axarr[1, 1].plot(t, sol[:,8]/mM, label="", lw=2)
    axarr[1, 1].set_ylabel("K+ pvs (mM)")
    axarr[2, 1].plot(t, sol[:,14]/uM, label="", lw=2)
    axarr[2, 1].set_ylabel("Ca2+ smc (uM)")
    axarr[3, 1].plot(t, sol[:,13]/(2*np.pi*um), label="", lw=2)
    axarr[3, 1].set_ylabel("r (um)")
    
    f.suptitle("time (s)", y=0.05)
    
    # Fine-tune figure; hide x ticks for top plots
    plt.setp([a.get_xticklabels() for a in axarr[0,:]], visible=False)
    plt.setp([a.get_xticklabels() for a in axarr[1,:]], visible=False)
    plt.setp([a.get_xticklabels() for a in axarr[2,:]], visible=False)

    # Fine-tune figure; make subplots farther from each other.
    f.subplots_adjust(wspace=0.3, hspace=0.2)
    if fname == '':
        plt.show()
    else:
        plt.savefig(fname, dpi=600, bbox_inches='tight')
        
        
def plot_solutions(t, sol1, sol2, fig_dims, uM, mV, mM, um, fname='', **kwargs):
    """Plot two solutions together.

    Parameters
    --------------
    t : array
        Time.
    sol : array
        Array containing model solution.
    fig_dims : tuple
        Figure dimensions.
    """
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.rcParams['legend.fontsize'] = 9
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.serif'] = ['Arial']
    
    f, axarr = plt.subplots(4, 2)
    f.set_size_inches(fig_dims[0], h=fig_dims[1])
    
    # left side
    axarr[0, 0].plot(t, sol2[:,0]/uM, label="", lw=2, ls='--', color='k')
    axarr[0, 0].plot(t, sol1[:,0]/uM, label="", lw=2)
    axarr[0, 0].set_ylabel("K+ syn (uM)")
    axarr[1, 0].plot(t, sol2[:,1]/uM, label="", lw=2, ls='--', color='k')
    axarr[1, 0].plot(t, sol1[:,1]/uM, label="", lw=2)
    axarr[1, 0].set_ylabel("IP3 (uM)")
    axarr[2, 0].plot(t, sol2[:,2]/uM, label="", lw=2, ls='--', color='k')
    axarr[2, 0].plot(t, sol1[:,2]/uM, label="", lw=2)
    axarr[2, 0].set_ylabel("Ca2+ ast (uM)")
    axarr[3, 0].plot(t, sol2[:,5]/uM, label="", lw=2, ls='--', color='k')
    axarr[3, 0].plot(t, sol1[:,5]/uM, label="", lw=2)
    axarr[3, 0].set_ylabel("EET (uM)")

    # right side
    axarr[0, 1].plot(t, sol2[:,7]/mV, label="", lw=2, ls='--', color='k')
    axarr[0, 1].plot(t, sol1[:,7]/mV, label="", lw=2)
    axarr[0, 1].set_ylabel("Vk (mV)")
    axarr[1, 1].plot(t, sol2[:,8]/mM, label="", lw=2, ls='--', color='k')
    axarr[1, 1].plot(t, sol1[:,8]/mM, label="", lw=2)
    axarr[1, 1].set_ylabel("K+ pvs (mM)")
    axarr[2, 1].plot(t, sol2[:,14]/uM, label="", lw=2, ls='--', color='k')
    axarr[2, 1].plot(t, sol1[:,14]/uM, label="", lw=2)
    axarr[2, 1].set_ylabel("Ca2+ smc (uM)")
    axarr[3, 1].plot(t, sol2[:,13]/(2*np.pi*um), label="", lw=2, ls='--', color='k')
    axarr[3, 1].plot(t, sol1[:,13]/(2*np.pi*um), label="", lw=2)
    axarr[3, 1].set_ylabel("r (um)")
    
    f.suptitle("time (s)", y=0.05)
    
    # Fine-tune figure; hide x ticks for top plots
    plt.setp([a.get_xticklabels() for a in axarr[0,:]], visible=False)
    plt.setp([a.get_xticklabels() for a in axarr[1,:]], visible=False)
    plt.setp([a.get_xticklabels() for a in axarr[2,:]], visible=False)

    # Fine-tune figure; make subplots farther from each other.
    f.subplots_adjust(wspace=0.3, hspace=0.2)
    if fname == '':
        plt.show()
    else:
        plt.savefig(fname, dpi=600, bbox_inches='tight')


def plot_vasodilation(t, sol, fig_dims, x_rel, um, uM, fname='', **kwargs):
    """Plot artery radial strain and perivascular Ca2+.

    Parameters
    --------------
    t : array
        Time.
    sol : array
        Array containing model solution.
    fig_dims : tuple
        Figure dimensions.
    x_rel : float
        Value of the relaxed vessel circumference.
    """
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.rcParams['legend.fontsize'] = 9
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.serif'] = ['Arial']

    r = sol[:,13]/(2*np.pi*um)
    r_rel = x_rel/(2*np.pi*um)
    strain = (r-r_rel)/r_rel
    
    f, axarr = plt.subplots(2, 1)
    f.set_size_inches(fig_dims[0], h=fig_dims[1])
    
    axarr[0].plot(t, strain, label="", lw=2)
    axarr[0].set_ylabel("Radial strain")
    axarr[1].plot(t, sol[:,9]/uM, label="", lw=2)
    axarr[1].set_ylabel("Ca2+ pvs (uM)")
    
    f.suptitle("time (s)", y=0.05)
    
    # Fine-tune figure; hide x ticks for top plots
    plt.setp(axarr[0].get_xticklabels(), visible=False)

    # Fine-tune figure; make subplots farther from each other.
    f.subplots_adjust(wspace=0.3, hspace=0.2)
    if fname == '':
        plt.show()
    else:
        plt.savefig(fname, dpi=600, bbox_inches='tight')
    
    
def plot_input(Jrho_IN, fig_dims, uM, s, fname='', **kwargs):
    """Plots the input function.

    Parameters
    --------------
    Jrho_IN : array
        Array containing model input values.
    fig_dims : tuple
        Figure dimensions.
    """
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.rcParams['legend.fontsize'] = 9
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.serif'] = ['Arial']
    
    f, ax1 = plt.subplots()
    f.set_size_inches(fig_dims[0], h=fig_dims[1]/3)

    ax1.plot(Jrho_IN[:,0], Jrho_IN[:,1]/(uM/s), label="K+", lw=2)
    ax1.plot(Jrho_IN[:,0], Jrho_IN[:,2], label="Glu", lw=2)
    ax1.set_ylabel("Glu (1) / K+ (uM/s)")
    ax1.set_xlabel("time (s)")
    ax1.legend()
    if fname == '':
        plt.show()
    else:
        plt.savefig(fname, dpi=600, bbox_inches='tight')
    
    
def plot_Vk(t, sol, fig_dims, r0, um, mV, fname='', **kwargs):
    """Plot astrocyte membrane potential and vessel dilation.

    Parameters
    --------------
    t : array
        Time.
    sol : array
        Array containing model solution.
    fig_dims : tuple
        Figure dimensions.
    """
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.rcParams['legend.fontsize'] = 9
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.serif'] = ['Arial']

    r = (sol[:,13]/(2*np.pi)-r0)/r0 * 100
    
    f, ax1 = plt.subplots()
    f.set_size_inches(fig_dims[0], h=fig_dims[1])

    ax1.plot(t, sol[:,7]/mV, label="", lw=2, color="C1")
    ax1.set_ylabel("Vk (mV)")
    ax1.set_xlabel("time (s)")

    ax2 = ax1.twinx()
    ax2.plot(t, r, lw=2)
    ax2.set_ylim([-10, 100])
    ax2.set_ylabel("vessel dilation (%)")
    
    if fname == '':
        plt.show()
    else:
        plt.savefig(fname, dpi=600, bbox_inches='tight')
        
        
def plot_currents(t, sol, fig_dims, pA, JSigKkNa, KKoa, Castr, gamma, gbk, vbk,
                  gtrpv, vtrpv, fname='', **kwargs):
    """Plot artery radial strain and perivascular Ca2+.

    Parameters
    --------------
    t : array
        Time.
    sol : array
        Array containing model solution.
    fig_dims : tuple
        Figure dimensions.
    x_rel : float
        Value of the relaxed vessel circumference.
    """
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.rcParams['legend.fontsize'] = 9
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.serif'] = ['Arial']
    
    potassium_s = sol[:,1]
    ss = sol[:,4]
    nbk = sol[:,6]
    Vk = sol[:,7]
    JSigK = JSigKkNa * potassium_s/(potassium_s + KKoa)
    Isigk = -JSigK * Castr * gamma
    Ibk = gbk * nbk * (Vk - vbk)
    Itrpv = gtrpv * ss * (Vk - vtrpv)
    
    f, ax1 = plt.subplots()
    f.set_size_inches(fig_dims[0], h=fig_dims[1]/3)
    
    ax1.plot(t, Isigk/pA, label="Isigk", lw=2)
    ax1.plot(t, Ibk/pA, label="Ibk", lw=2)
    ax1.plot(t, Itrpv/pA, label="Itrpv", lw=2)
    ax1.legend()
    ax1.set_ylabel("current (pA)")
    ax1.set_xlabel("time (s)")
    
    if fname == '':
        plt.show()
    else:
        plt.savefig(fname, dpi=600, bbox_inches='tight')
