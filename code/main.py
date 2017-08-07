# -*- coding: utf-8 -*-

"""
Code for Diem AK (2017) Chemical Signalling in the Neurovascular Unit, The ReScience Journal
https://github.com/akdiem/ReScience-submission/tree/Diem-2017
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
    it1 = int(5.0*sec)
    it2 = int(sec)
    it3 = int(18*sec)
    it4 = int(sec)
    pos = it1
    Jrho_IN[pos+1:pos+it2+1,1] = Max_neural_Kplus * np.linspace(0, 1, it2)
    Jrho_IN[pos+1:pos+it2+1,2] = Max_neural_glut * np.linspace(0, 1, it2)
    pos += it2
    Jrho_IN[pos+1:pos+it3+1,1] = Max_neural_Kplus * np.ones(it3)
    Jrho_IN[pos+1:pos+it3+1,2] = Max_neural_glut * np.ones(it3)
    pos += it3
    Jrho_IN[pos+1:pos+it4+1,1] = Max_neural_Kplus * np.linspace(1, 0, it4)
    Jrho_IN[pos+1:pos+it4+1,2] = Max_neural_glut * np.linspace(1, 0, it4)
    return Jrho_IN


def plot_input(Jrho_IN, fig_dims, uM, s, **kwargs):
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.rcParams['legend.fontsize'] = 9
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.serif'] = ['Arial']
    
    f, ax1 = plt.subplots()
    f.set_size_inches(fig_dims[0], h=fig_dims[1])

    ax1.plot(Jrho_IN[:,0], Jrho_IN[:,1]/(uM/s), label="K+", lw=2)
    ax1.plot(Jrho_IN[:,0], Jrho_IN[:,2], label="Glu", lw=2)
    ax1.set_ylabel("Glu (1) / K+ (uM/s)")
    ax1.set_xlabel("time (s)")
    ax1.legend()
#    plt.savefig('../article/figures/input.png', dpi=600, bbox_inches='tight')
    plt.show()
    

def main(fparam, fig_dims):
    units, param = utils.read_config(fparam)
    um = units['um']
    
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
    sol = nvu.run_simulation(t, y0, Jrho_IN, x_rel, units, param, atol=atol, rtol=rtol)
    y0 = sol[-1,:]
    
    # Plot solution
    #nvu.plot_solution(t, sol, fig_dims, **units)
    
    # Simulation
    t1 = 0
    t2 = 50 
    nt = 200
    Jrho_IN = K_glut_release(t1, t2, **units)
    t = np.linspace(t1, t2, nt)    
    sol = nvu.run_simulation(t, y0, Jrho_IN, x_rel, units, param, atol=atol, rtol=rtol)
    
#    plot_input(Jrho_IN, fig_dims, **units)
    
    # Plot solution
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
