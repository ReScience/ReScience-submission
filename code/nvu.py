# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from scipy.integrate import odeint, ode
from scipy.interpolate import interp1d
import matplotlib.pylab as plt
import sys


# Units
s = 1
cm = 1e-2
um = 1e-4*cm
mM = 1e-3
uM = 1e-3 * mM
nM = 1e-6 * mM
C = 1
A = C/s
V = 1
mV = 1e-3 * V
F = C/V
pF = 1e-12 * F
S = A/V
pS = 1e-12 * S
dyn = 1
mmHg = 1333.22 * dyn/cm**2


# Parameter
JSigKkNa = 0.6 * mM/s
KKoa = 1.5 * mM
rh = 4.8 * uM
kdeg = 1.25 * 1/s
KG = 8.82
delta = 0.001235
beta = 0.0244
Jmax = 2880 * uM/s
Ki = 0.03 * uM
Kact = 0.17 * uM
Vmax = 20 * uM/s
Kp = 0.24 * uM
Pl = 5.2 * uM/s
calcium_er = 400 * uM # value from Bennet et al.
Castr = 40 * pF
gamma = 1970 * mV/uM
gtrpv = 200 * pS
vtrpv = 6 * mV
kon = 2 * 1/(uM*s)
Kinh = 0.1 * uM
tautrpv = 0.9 * 1/s
eps12 = 0.16
kappa = 0.04
gammacai = 0.2 * uM
gammacae = 0.2 * mM
v1trpv = 120 * mV
v2trpv = 13 * mV
Veet = 72 * 1/s
calcium_a_min = 0.1 * uM
keet = 7.1 * 1/s
psibk = 2.664 * 1/s
v4bk = 14.5 * mV
v5bk = 8 * mV
v6bk = -15 * mV
eetshift = 2 * mV/uM
Ca3bk = 400 * nM
Ca4bk = 150 * nM
gbk = 225.6 * pS
vbk = -95 * mV
gleak = 78.54 * pS
vleak = -70 * mV
Csmc = 19.635 * pF
gkir0 = 145 * pS
vkir1 = 57 * mV
vkir2 = 130 * mV
VRpa = 3.2e-5
VRps = 0.1
gca = 157 * pS
vca = 80 * mV
v2 = 25 * mV
dp = 60 * mmHg
Rdecay = 1 * 1/s
potassium_p_min = 3 * mM
Cadecay = 0.5 * 1/s
calcium_p_min = 2000 * uM
alphakir = 1020 * s
av1 = 18 * mV
av2 = 10.8 * mV
betakir = 26.9 * s
bv1 = 18 * mV
bv2 = 0.06 * mV
gl = 62.832 * pS
vl = -70 * mV
gk = 251.33 * pS
vk = -80 * mV
phin = 2.664
Ca3 = 400 * nM
Ca4 = 150 * nM
v4 = 14.5 * mV
v5 = 8 * mV
v6 = -15 * mV
tau = 0.2 * dyn/cm
Kd = 1000 * nM
Bt = 10000 * nM
alpha = 4.3987e15 * nM/C
kca = 1.3568e11 * nM/C
we = 0.9
x0 = 188.5 * um
x1 = 1.2
x2 = 0.13
x3 = 2.2443
x4 = 0.71182
x5 = 0.8
x6 = 0.01
x7 = 0.32134
x8 = 0.88977
x9 = 0.0090463
sigma0h = 3e6 * dyn/cm**2
u1 = 41.76
u2 = 0.047396
u3 = 0.0584
wm = 0.7
kpsi = 3.3
Cam = 500 * nM
psim = 0.3
q = 3
y0 = 0.928
y1 = 0.639
y2 = 0.35
y3 = 0.78847
y4 = 0.8
sigmay0h = 2.6e6 * dyn/cm**2
vref = 0.24
Caref = 510 * nM
ad = 0.28125
bd = 5
cd = 0.03
dd = 1.3
Sx = 40000 * um**2
Ax = 50.265 * um**2


def nvu(t, y, Jrho_IN, x_rel):
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
    JSigK = JSigKkNa * potassium_s/(potassium_s + KKoa)
    JKss = interp1d(Jrho_IN[:,0], Jrho_IN[:,1], bounds_error=False, fill_value=0)
    potassium_s_dt = JKss(t) - JSigK
    
    
    # Astrocytic space
    rhos = interp1d(Jrho_IN[:,0], Jrho_IN[:,2], bounds_error=False, fill_value=0)
    G = (rhos(t) + delta)/(KG + rhos(t) + delta)
    ip3_dt = rh*G - kdeg*ip3
    Jip3 = Jmax * ((ip3/(ip3+Ki)) * (calcium_a/(calcium_a+Kact)) * h)**3 *\
        (1 - calcium_a/calcium_er)
    Jpump = Vmax * calcium_a**2 / (calcium_a**2 + Kp**2)
    Jleak = Pl * (1 - calcium_a/calcium_er)
    Itrpv = gtrpv * ss * (Vk - vtrpv)
    Jtrpv = -Itrpv/(Castr*gamma)
    calcium_a_dt = beta * (Jip3 - Jpump + Jleak + Jtrpv)
    h_dt = kon * (Kinh - (calcium_a + Kinh) * h)
    tauca = tautrpv / (calcium_p/uM)
    eps = (x - x_rel)/x_rel
    Hca = calcium_a/gammacai + calcium_p/gammacae
    # sinf is inconsistent between code and equations
    # gammacai is off by order of magnitude
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
    Vk_dt = 1/Castr * (-Isigk - Ibk - Ileak - Itrpv)
    
    
    # Perivascular space
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
    Jca = Ica/(Csmc*gamma)
    calcium_p_dt = -Jtrpv - Jca - Cadecay * (calcium_p - calcium_p_min)
    
    
    # Ion currents
    alphak = alphakir / (1 + np.exp((Vm - vkir + av1)/av2))
    betak = betakir * np.exp(bv2/mV * (Vm - vkir + bv1)/mV)
    tauk = 1/(alphak+betak)
    kinf = alphak/(alphak+betak)
    k_dt = 1/tauk * (kinf - k)
    Il = gl * (Vm - vl)
    Ik = gk * n * (Vm - vk)
    Vm_dt = (1/Csmc) * (-Il - Ik - Ica - Ikir)
    v3 = -(v5/2) * np.tanh((calcium_smc-Ca3)/Ca4) + v6
    lamn = phin * np.cosh(0.5*(Vm-v3)/v4)
    ninf = 0.5 * (1 + np.tanh((Vm-v3)/v4))
    n_dt = lamn * (ninf - n)
    
    # Vessel SMC calcium
    rho_smc = (Kd+calcium_smc)**2/((Kd+calcium_smc)**2 + Kd*Bt)
    calcium_smc_dt = -rho_smc * (alpha*Ica + kca*calcium_smc)
    
    
    # Vessel mechanics
    fdp = 0.5 * dp * (x/np.pi - Ax/x) * cm # why *cm here? Is that the artery length?
    xd = x/x0
    sigmax = x3*(1 + np.tanh((xd-x1)/x2)) + x4*(xd-x5) - x8*(x6/(xd-x7))**2 - x9
    fx = we*Sx*sigmax*sigma0h
    u = x-yy
    ud = u/x0
    sigmau = u2 * np.exp(u1*ud) - u3
    fu = wm*Sx*sigmau*sigma0h
    x_dt = 1/tau * (fdp - fx - fu)
    psi = calcium_smc**q/(Cam**q+calcium_smc**q)
    omega_dt = kpsi * (psi/(psim+psi) - omega)
    yd = yy/x0
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
    
    return [potassium_s_dt, ip3_dt, calcium_a_dt, h_dt, ss_dt, eet_dt, nbk_dt,
            Vk_dt, potassium_p_dt, calcium_p_dt, k_dt, Vm_dt, n_dt, x_dt,
            calcium_smc_dt, omega_dt, yy_dt]


def init():
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
    x = 2*np.pi*20*um
    calcium_smc = 3.41385670857693e-07
    omega = 0.536911672725179
    yy = 0.000115089683436595
    return [potassium_s, ip3, calcium_a, h, ss, eet, nbk, Vk, potassium_p,
            calcium_p, k, Vm, n, x, calcium_smc, omega, yy]
    
    
def K_glut_release(t1, t2):
    sizeJrho = 1600
    sec = sizeJrho/(t2-t1)
    Max_neural_Kplus = 0.48*uM/s
    Max_neural_glut = 0.5
    Jrho_IN = np.zeros((sizeJrho,3))
    Jrho_IN[:,0] = np.linspace(t1, t2, sizeJrho)
    it1 = int(2.5*sec)
    it2 = int(sec)
    it3 = int(15*sec)
    it4 = int(0.25*sec)
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
    
    
def main(fig_dims):
    y0 = init()
    x_rel = y0[13]
    
#    # Equilibration
#    t1 = -20
#    t2 = 0
#    n = 100
#    Jrho_IN = np.zeros((n,3))
#    Jrho_IN[:,0] = np.linspace(t1, t2, n)
#    t = np.linspace(t1, t2, 200)    
#    ode15s = ode(nvu)
#    ode15s.set_f_params(Jrho_IN, x_rel)
#    ode15s.set_integrator('lsoda')
#    ode15s.set_initial_value(y0, t=t1)
#    nt = len(t)
#    sol = np.zeros([nt, len(y0)])
#    sol[0,:] = y0
#    for i in range(1,nt):
#        if ode15s.successful():
#            sol[i,:] = ode15s.integrate(t[i])
#    y0 = sol[-1,:]
#    
#    plt.figure(figsize=fig_dims)
#    plt.plot(t, sol[:,7]/mV, label="", lw=2)
#    plt.xlabel("time")
#    plt.ylabel("")
#    plt.show()
#    
#    sys.exit()
    

    
    # Simulation
    t1 = 0
    t2 = 50 
    Jrho_IN = K_glut_release(t1, t2)
    t = np.linspace(t1, t2, 200)    
    ode15s = ode(nvu)
    ode15s.set_f_params(Jrho_IN, x_rel)
    ode15s.set_integrator('lsoda')
    ode15s.set_initial_value(y0, t=t1)
    nt = len(t)
    sol = np.zeros([nt, len(y0)])
    sol[0,:] = y0
    for i in range(1,nt):
        if ode15s.successful():
            sol[i,:] = ode15s.integrate(t[i])    
    
    plt.figure(figsize=fig_dims)
    plt.plot(t, sol[:,0]/uM, label="", lw=2)
    plt.xlabel("time")
    plt.ylabel("")
    plt.show()
    
    plt.figure(figsize=fig_dims)
    plt.plot(t, sol[:,1]/uM, label="", lw=2)
    plt.xlabel("time")
    plt.ylabel("")
    plt.show()
    
    plt.figure(figsize=fig_dims)
    plt.plot(t, sol[:,2]/uM, label="", lw=2)
    plt.xlabel("time")
    plt.ylabel("")
    plt.show()
    
    plt.figure(figsize=fig_dims)
    plt.plot(t, sol[:,5]/uM, label="", lw=2)
    plt.xlabel("time")
    plt.ylabel("")
    plt.show()
    
    plt.figure(figsize=fig_dims)
    plt.plot(t, sol[:,7]/mV, label="", lw=2)
    plt.xlabel("time")
    plt.ylabel("")
    plt.show()
    
    plt.figure(figsize=fig_dims)
    plt.plot(t, sol[:,8]/mM, label="", lw=2)
    plt.xlabel("time")
    plt.ylabel("")
    plt.show()
    
    plt.figure(figsize=fig_dims)
    plt.plot(t, sol[:,14]/uM, label="", lw=2)
    plt.xlabel("time")
    plt.ylabel("")
    plt.show()
    
    plt.figure(figsize=fig_dims)
    plt.plot(t, sol[:,13]/(2*np.pi*um), label="", lw=2)
    plt.xlabel("time")
    plt.ylabel("")
    plt.show()
    
    
    
    
if __name__ == "__main__":
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.rcParams['legend.fontsize'] = 9
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.serif'] = ['Arial']
    
    WIDTH = 510  # the number latex snp.pits out
    FACTOR = 0.5  # the fraction of the width you'd like the figure to occupy
    fig_width_pt  = WIDTH * FACTOR
    inches_per_pt = 1.0 / 72.27
    golden_ratio  = (np.sqrt(5) - 1.0) / 2.0  # because it looks good
    fig_width_in  = fig_width_pt * inches_per_pt  # figure width in inches
    fig_height_in = fig_width_in * golden_ratio   # figure height in inches
    fig_dims    = [fig_width_in, fig_height_in] # fig dims as a list
    
    main(fig_dims)