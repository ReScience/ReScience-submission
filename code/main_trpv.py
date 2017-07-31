# -*- coding: utf-8 -*-

from nvu import nvu, utils

import numpy as np
from sys import argv
from scipy.interpolate import interp1d
import matplotlib.pylab as plt


def init(r0):
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


def plot_Vk(t, sol, fig_dims, r0, um=0, mV=0, **kwargs):
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
    
#    print(sol[0,7]/mV)
    
    plt.savefig('../article/figures/Vk_inflation.png', dpi=600, bbox_inches='tight')
#    plt.show()


def main(fparam, fig_dims):
    units, param = utils.read_config(fparam)
    um = units['um']
    mV = units['mV']
    
    r0 = 20 * um
    y0 = init(r0)
    x_rel = y0[13]
    sol = np.zeros(len(y0))

    atol = 1e-7
    rtol = 1e-8

    # manual artery inflation according to Cao2011
    cao = np.loadtxt("../data/cao2011.csv", delimiter=',')
    x0 = 2*r0*np.pi
    r_dil = r0 + r0*cao[:,1]/100
    r_f = interp1d(cao[:,0], r_dil, bounds_error=False, fill_value=r0)
    y0[13] = 2*np.pi*r_f(0)

    # Equilibration
    t1 = -20
    t2 = 0
    nt = 100
    Jrho_IN = np.zeros((nt,3))
    Jrho_IN[:,0] = np.linspace(t1, t2, nt)
    t = np.linspace(t1, t2, nt)
    sol = nvu.run_simulation(t, y0, Jrho_IN, x_rel, r_f, units, param, atol=atol, rtol=rtol, mode='trpv')
    y0 = sol[-1,:]
    
    # Simulation
    t1 = 0
    t2 = 240 
    nt = 2000
    Jrho_IN = np.zeros((nt,3))
    Jrho_IN[:,0] = np.linspace(t1, t2, nt)
    t = np.linspace(t1, t2, nt)    
    sol = nvu.run_simulation(t, y0, Jrho_IN, x_rel, r_f, units, param, atol=atol, rtol=rtol, mode='trpv')
    sol[:,13] = 2*np.pi*r_f(t)

    # Plot solution
    plot_Vk(t, sol, fig_dims, r0, **units)
    

if __name__ == "__main__":
    script, param = argv

    np.seterr(over='ignore')
    
    WIDTH = 510  # the number latex spits out
    FACTOR = 1.0  # the fraction of the width you'd like the figure to occupy
    fig_width_pt  = WIDTH * FACTOR
    inches_per_pt = 1.0 / 72.27
    golden_ratio  = (np.sqrt(5) - 1.0) / 3  # because it looks good
    fig_width_in  = fig_width_pt * inches_per_pt  # figure width in inches
    fig_height_in = fig_width_in * golden_ratio   # figure height in inches
    fig_dims    = [fig_width_in, fig_height_in] # fig dims as a list
    
    main(param, fig_dims)