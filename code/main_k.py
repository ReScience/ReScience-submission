# -*- coding: utf-8 -*-

"""
Code for Diem AK (2017) Chemical Signalling in the Neurovascular Unit, The ReScience Journal
https://github.com/akdiem/ReScience-submission/tree/Diem-2017
"""

from nvu import nvu, utils

import numpy as np
from sys import argv
from scipy.interpolate import interp1d


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


def k(t):
    # force k to fit equation (40)
    return 0.01 * (1 + np.tanh((t-270)/60))


def main(fparam, fig_dims):
    units, param = utils.read_config(fparam)
    um = units['um']
    mV = units['mV']
    
    r0 = 20 * um
    y0 = init(r0)
    x_rel = y0[13]
    sol = np.zeros(len(y0))

    atol = 1e-7
    rtol = 1e-7

    # Equilibration
    t1 = -20
    t2 = 0
    nt = 100
    Jrho_IN = np.zeros((nt,3))
    Jrho_IN[:,0] = np.linspace(t1, t2, nt)
    t = np.linspace(t1, t2, nt)
    sol = nvu.run_simulation(t, y0, Jrho_IN, x_rel, k, units, param, atol=atol, rtol=rtol, mode='k')
    y0 = sol[-1,:]
    
    # Simulation
    t1 = 0
    t2 = 350 
    nt = 2000
    t = np.linspace(t1, t2, nt)    
    sol = nvu.run_simulation(t, y0, Jrho_IN, x_rel, k, units, param, atol=atol, rtol=rtol, mode='k')

    # Plot solution
    nvu.plot_vasodilation(t, sol, fig_dims, x_rel, **units)
    

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
