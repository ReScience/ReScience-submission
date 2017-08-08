# -*- coding: utf-8 -*-

"""
Code for Diem AK (2017) Chemical Signalling in the Neurovascular Unit, The
ReScience Journal

https://github.com/akdiem/ReScience-submission/tree/Diem-2017

This is a reimplementation of Witthoft A, Karniadakis GE (2012) A bidirectional
model for communication in the neurovascular unit. Journal of Theoretical
Biology 311: 80-93.

This file initiates and runs a simulation of the neurovascular unit.
"""

from nvu import nvu, utils

import numpy as np
from sys import argv
import matplotlib.pylab as plt 


def init(r0):
    """Returns the initial conditions for the NVU model.

    Parameters
    --------------
    r0 : float
        Blood vessel radius at rest.

    Returns
    ---------
    list
        List of initial values for the variables of nvu()
    """
    potassium_s = 2.92655044308714e-8
    ip3 = 5.37611796987610e-10
    calcium_a = 1.47220569018281e-07
    h = 0.404507631346124
    ss = 0.0161921297424289
    eet = 4.78801348065449e-07
    nbk = 6.24930194169376e-5
    Vk = -0.0814061063457068
    potassium_p = 0.00353809145071707
    calcium_p = 4.60269585230500e-06
    k = 8.01125818473096e-09
    Vm = 8.33004194103223e-05
    n = 0.283859572906570
    x = 2*np.pi*r0
    calcium_smc = 3.41385670857693e-07
    omega = 0.536911672725179
    yy = 0.000115089683436595
    return [potassium_s, ip3, calcium_a, h, ss, eet, nbk, Vk, potassium_p,
            calcium_p, k, Vm, n, x, calcium_smc, omega, yy]
    
    
def K_glut_release(t1, t2, uM, s, **kwargs):
    """Returns the input conditions for the NVU model.
    
    Input variables are the K+ (uM/s) concentration and ratio of bound/unbound
    Glu receptors. They are realised as step functions with a steep linear increase
    to a maximum sustained value and steep linear decrease back to 0 once the
    stimulus is removed.

    Parameters
    --------------
    t1 : float
        Time of the beginning of the simulation.
    t2 : float
        Time of the end of the simulation.
    uM : float
        Value for unit uM.
    s : float
        Value for unit s.

    Returns
    ---------
    array
        Array containing input conditions over simulation time.
    """
    sizeJrho = 1600
    sec = sizeJrho/(t2-t1)
    Max_neural_Kplus = 0.55*uM/s
    Max_neural_glut = 0.5
    Jrho_IN = np.zeros((sizeJrho,3))
    Jrho_IN[:,0] = np.linspace(t1, t2, sizeJrho)
    pulse = 20
    it1 = int(0*sec)
    it2 = int(pulse*sec)
    it3 = int((25-pulse)*sec)
    it4 = int(1*sec)
    pos = it1
    t = np.linspace(0, pulse, it2)
    Jrho_IN[pos+1:pos+it2+1,1] = Max_neural_Kplus * 0.5*(1 + np.tanh((t-9)/3))
    Jrho_IN[pos+1:pos+it2+1,2] = Max_neural_glut * 0.5*(1 + np.tanh((t-9)/3))
    pos += it2
    Jrho_IN[pos+1:pos+it3+1,1] = Max_neural_Kplus * np.ones(it3)
    Jrho_IN[pos+1:pos+it3+1,2] = Max_neural_glut * np.ones(it3)
    pos += it3
    Jrho_IN[pos+1:pos+it4+1,1] = Max_neural_Kplus * np.linspace(1, 0, it4)
    Jrho_IN[pos+1:pos+it4+1,2] = Max_neural_glut * np.linspace(1, 0, it4)
    return Jrho_IN


def main(fparam, fig_dims):
    # read units and parameters from parameter file
    units, param = utils.read_config(fparam)
    um = units['um']
    
    # set up simulation variables
    r0 = 20 * um
    y0 = init(r0)
    x_rel = y0[13]
    sol = np.zeros(len(y0))

    # simulation tolerance values
    atol = 1e-7
    rtol = 1e-7

    # Equilibration
    t1 = -20
    t2 = 0
    nt = 100
    Jrho_IN = np.zeros((nt,3))
    Jrho_IN[:,0] = np.linspace(t1, t2, nt)
    t = np.linspace(t1, t2, nt)
    sol = nvu.run_simulation(t, y0, Jrho_IN, x_rel, units, param, atol=atol, rtol=rtol)
    y0 = sol[-1,:]
    
    # Simulation
    t1 = 0
    t2 = 50 
    nt = 200
    Jrho_IN = K_glut_release(t1, t2, **units)
    t = np.linspace(t1, t2, nt)    
    sol = nvu.run_simulation(t, y0, Jrho_IN, x_rel, units, param, atol=atol, rtol=rtol)
    
    #nvu.plot_input(Jrho_IN, fig_dims, fname='../article/figures/input.png', **units)
#    nvu.plot_input(Jrho_IN, fig_dims, **units)
    
    # Plot solution
    #nvu.plot_solution(t, sol, fig_dims, fname='../article/figures/fig1.png', **units)
    nvu.plot_solution(t, sol, fig_dims, **units)
    

if __name__ == "__main__":
    script, param = argv

    np.seterr(over='ignore')
    
    WIDTH = 510  # the number latex spits out
    FACTOR = 1.0  # the fraction of the width you'd like the figure to occupy
    fig_width_pt  = WIDTH * FACTOR
    inches_per_pt = 1.0 / 72.27
    golden_ratio  = (np.sqrt(5) - 1.0) / 2  # because it looks good
    fig_width_in  = fig_width_pt * inches_per_pt  # figure width in inches
    fig_height_in = fig_width_in * golden_ratio   # figure height in inches
    fig_dims    = [fig_width_in, fig_height_in] # fig dims as a list
    
    main(param, fig_dims)
