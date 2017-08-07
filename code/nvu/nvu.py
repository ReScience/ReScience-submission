# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from scipy.integrate import ode
from scipy.interpolate import interp1d
import matplotlib.pylab as plt


def synapse(t, potassium_s, Jrho_IN, JSigKkNa, KKoa, **kwargs):
    JSigK = JSigKkNa * potassium_s/(potassium_s + KKoa)
    JKss = interp1d(Jrho_IN[:,0], Jrho_IN[:,1], bounds_error=False, fill_value=0)
    potassium_s_dt = JKss(t) - JSigK
    return potassium_s_dt, JSigK


def astrocyte(t, ip3, calcium_a, h, ss, Vk, calcium_p, x, eet, nbk, Jrho_IN,
              x_rel, JSigK, delta, KG, rh, kdeg, Jmax, Ki, Kact, calcium_er,
              Vmax, Kp, Pl, gtrpv, vtrpv, Castr, gamma, beta, kon, Kinh, 
              tautrpv, uM, gammacai, gammacae, eps12, kappa, v1trpv, v2trpv,
              Veet, calcium_a_min, keet, v5bk, Ca3bk, Ca4bk, v6bk, psibk,
              v4bk, eetshift, gbk, vbk, gleak, vleak, **kwargs):
    rhos = interp1d(Jrho_IN[:,0], Jrho_IN[:,2], bounds_error=False,
                    fill_value=0)
    G = (rhos(t) + delta)/(KG + rhos(t) + delta)
    ip3_dt = rh*G - kdeg*ip3
    Jip3 = Jmax * (ip3/(ip3+Ki) * calcium_a/(calcium_a+Kact) * h)**3 *\
        (1 - calcium_a/calcium_er)
    Jpump = Vmax * calcium_a**2 / (calcium_a**2 + Kp**2)
    Jleak = Pl * (1 - calcium_a/calcium_er)
    Itrpv = gtrpv * ss * (Vk - vtrpv)
    Jtrpv = -Itrpv/(Castr*gamma)
    calcium_a_dt = beta * (Jip3 - Jpump + Jleak)
    h_dt = kon * (Kinh - (calcium_a + Kinh) * h)
    tauca = tautrpv / (calcium_p/uM)
    eps = (x - x_rel)/x_rel
    Hca = calcium_a/gammacai + calcium_p/gammacae
    sinf = (1/(1 + np.exp(-(eps-eps12)/kappa))) * (1/(1+Hca) *\
        (Hca + np.tanh((Vk - v1trpv)/v2trpv)))
    ss_dt = 1/tauca * (sinf - ss)
    eet_dt = Veet * (calcium_a - calcium_a_min) - keet*eet
    v3bk = -(v5bk/2) * np.tanh((calcium_a-Ca3bk)/Ca4bk) + v6bk
    phibk = psibk * np.cosh((Vk-v3bk)/(2*v4bk))
    ninf = 0.5 * (1 + np.tanh((Vk + eetshift*eet - v3bk)/v4bk))
    nbk_dt = phibk * (ninf - nbk)
    Ibk = gbk * nbk * (Vk - vbk)
    Ileak = gleak * (Vk - vleak)
    Isigk = -JSigK * Castr * gamma
    Vk_dt = 1/Castr * (-Isigk - Ibk - Ileak)
    return ip3_dt, calcium_a_dt, h_dt, ss_dt, eet_dt, nbk_dt, Vk_dt, Ibk, Jtrpv


def perivascular_space(potassium_p, k, Vm, calcium_p, Ibk, Jtrpv, Castr,
                       gamma, gkir0, mM, vkir1, vkir2, Csmc, VRpa, VRps,
                       Rdecay, potassium_p_min, dp, mmHg, mV, v2, gca, vca,
                       Cadecay, calcium_p_min, **kwargs):
    Jbk = Ibk/(Castr*gamma)
    gkir = gkir0 * np.sqrt(potassium_p/mM)
    vkir = vkir1 * np.log10(potassium_p/mM) - vkir2
    Ikir = gkir * k * (Vm - vkir)
    Jkir = Ikir/(Csmc*gamma)
    potassium_p_dt = Jbk/VRpa + Jkir/VRps - Rdecay * (potassium_p -\
        potassium_p_min)
    v1 = (-17.4-(12*(dp/mmHg)/200))*mV
    minf = 0.5 * (1 + np.tanh((Vm-v1)/v2))
    Ica = gca * minf * (Vm - vca)
    Jca = -Ica/(Csmc*gamma)
    calcium_p_dt = - Jca/VRps - Cadecay * (calcium_p - calcium_p_min)
    return potassium_p_dt, calcium_p_dt, vkir, Ikir, Ica


def ion_currents(Vm, k, calcium_smc, n, vkir, Ica, Ikir, alphakir, av1, av2,
                 betakir, bv2, bv1, mV, gl, vl, gk, vk, Csmc, v5, Ca3, Ca4, v6,
                 phin, v4, **kwargs):
    alphak = alphakir / (1 + np.exp((Vm - vkir + av1)/av2))
    betak = betakir * np.exp(bv2/mV * (Vm - vkir + bv1)/mV)
    tauk = 1/(alphak+betak)
    kinf = alphak/(alphak+betak)
    k_dt = 1/tauk * (kinf - k)
    Il = gl * (Vm - vl)
    Ik = gk * n * (Vm - vk)
    Vm_dt = (1/Csmc) * (-Il - Ik - Ica - Ikir)
    v3 = -(v5/2) * np.tanh((calcium_smc-Ca3)/Ca4) + v6
    lamn = phin * np.cosh((Vm-v3)/(2*v4))
    ninf = 0.5 * (1 + np.tanh((Vm-v3)/v4))
    n_dt = lamn * (ninf - n)
    return k_dt, Vm_dt, n_dt


def vessel_mechanics(t, calcium_smc, x, yy, omega, Ica, Kd, Bt, alpha, kca, dp,
                     Ax, um, x0, x3, x1, x2, x4, x5, x8, x6, x7, x9, we, Sx,
                     sigma0h, u2, u1, u3, wm, tau, Cam, q, kpsi, psim,
                     Caref, sigmay0h, y1, y2, y4, y0, y3, vref, ad, cd, bd, dd,
                     **kwargs):
    # SMC calcium
    rho_smc = (Kd+calcium_smc)**2/((Kd+calcium_smc)**2 + Kd*Bt)
    calcium_smc_dt = -rho_smc * (alpha*Ica + kca*calcium_smc)

    # Vessel mechanics
    fdp = 0.5 * dp * (x/np.pi - Ax/x) * um
    xd = x/x0
    sigmax = x3*(1 + np.tanh((xd-x1)/x2)) + x4*(xd-x5) - x8*(x6/(xd-x7))**2 - x9
    fx = we*Sx*sigmax*sigma0h
    yd = yy/x0
    ud = xd-yd
    sigmau = u2 * np.exp(u1*ud) - u3
    fu = wm*Sx*sigmau*sigma0h
    x_dt = 1/tau * (fdp - fx - fu)
    psi = calcium_smc**q/(Cam**q+calcium_smc**q)
    omega_dt = kpsi * (psi/(psim+psi) - omega)
    psiref = Caref**q/(Cam**q+Caref**q)
    omega_ref = psiref/(psim + psiref)
    sigmay0 = sigmay0h * omega/omega_ref
    sy = (y1/(yd+y2))**y4
    sigmay = sigmay0/sigma0h * (np.exp(-(yd-y0)**2/(2*sy**2)) - y3)/(1-y3)
    cond = sigmau/sigmay
    if cond < 1:
        ycond = -vref * (psi/psiref) * ad * (1-cond)/(ad+cond)
    else:
        ycond = cd*(np.exp(bd * (cond-dd)) - np.exp(bd * (1-dd)))
    yy_dt = x0*ycond
    
    return x_dt, calcium_smc_dt, omega_dt, yy_dt


def nvu(t, y, Jrho_IN, x_rel, units, param):
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


def run_simulation(time, y0, *args, atol=1e-7, rtol=1e-7, mode=''):
    integrator = "lsoda"
    ode15s = ode(nvu)
    if mode == 'Vk':
        ode15s = ode(nvu_Vk)
    elif mode == 'k':
        ode15s = ode(nvu_k)
    elif mode == 'trpv':
        ode15s = ode(nvu_trpv)
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


def plot_solution(t, sol, fig_dims, uM=0, mV=0, mM=0, um=0, **kwargs):
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
#    plt.savefig('../article/figures/no_trpv.png', dpi=600, bbox_inches='tight')
    plt.show()


def plot_vasodilation(t, sol, fig_dims, x_rel, um=0, uM=0, **kwargs):
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
    plt.savefig('../article/figures/pinacidil.png', dpi=600, bbox_inches='tight')
#    plt.show()
    
    
def plot_input(Jrho_IN, fig_dims, uM, s, **kwargs):
    """Plots the input function.

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
    
    
def plot_Vk(t, sol, fig_dims, r0, um, mV, **kwargs):
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
    
#    plt.savefig('../article/figures/Vk_inflation.png', dpi=600, bbox_inches='tight')
    plt.show()
