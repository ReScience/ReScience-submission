# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Andrei Maksimov (maksimov.andrei7@gmail.com)
# -----------------------------------------------------------------------------
# References:
#
# *Cellular and network mechanisms of slow oscillatory activity (<1 Hz)
# and wave propagations in a cortical network model*, A. Compte,
# M.V. Sanchez-Vives, D.A. McCormick, X.-J. Wang,
# Journal of Neurophysiology, 2707--2725, 2003"
# -----------------------------------------------------------------------------
# File description:
#
# Simulate response of pyramidal neurons to DC injection (Fig. 8)
# -----------------------------------------------------------------------------

import numpy as np
import awesomizer as aw
import matplotlib.pyplot as plt
import copy
import nest

nest.ResetKernel()

style = aw.style_Plos()
aw.style_apply(style)


###################################################
###############  PARAMETERS  ######################
###################################################

sim_time = 25000.  # [ms] simulation time
dt = 0.1  # [ms] simulation time step

# parameters for dc current injection
Param_dc = {
    'start': 2000.,  # [ms] start time of DC injection
    'stop': 22000.,  # [ms] stop time of DC injection
}

# [pA] range of amplitudes for DC injection for e and i neurons
temp = np.append(np.arange(-10., 1.8, 0.5), np.arange(1.8, 2., 0.02))
Amplitude_dc = np.append(temp, np.arange(2., 3.6, 0.2))


syns = nest.GetDefaults('compte2003_ex')['receptor_types']  # receptor IDs
C_m = 1E6  # [pF/cm^2] unit membrane capacitance

# parameters for exc neurons
A_s = 0.015E-2  # [cm^2] somatic area
A_d = 0.035E-2  # [cm^2] dendritic area
Params_neuron = {'PY':
                 {
                     't_ref': 2.,    # [ms] absolute refractory time 
                     'E_Ca':  120.,  # [mV] Ca2+ reversal potential
                     'E_K': -100.,   # [mV] K+ reversal potential
                     'E_L': -60.95,  # [mV] resting potential
                     'E_Na':   55.,  # [mV] Na+ reversal potential
                     'E_ex':    0.,  # [mV] AMPA, NMDA reversal potential
                     'E_in': -70.,   # [mV] GABA reversal potential
                     'g_conn': 1.75E3,  # [nS] axial dendritic conductance
                     'tau_syn_ampa': 2.,  # [ms] AMPA time constant
                     'tau_syn_nmda_fast': 2.,  # [ms] fast NMDA time constant
                     'tau_syn_nmda_slow': 100.,  # [ms] slow NMDA time constant
                     'tau_syn_gaba': 10.,  # [ms] GABA time constant

                     # parameters for somatic compartment
                     'C_m_s': C_m * A_s,  # [pF] capacitance
                     'V_m_s': -75.,  # [mV] initial membrane potential
                      # [nS] conductance of fast Na+ channel
                      'g_Na': 50E6 * A_s,
                      # [nS] conductance of fast K+ channel
                      'g_K': 10.5E6 * A_s,
                      'g_L': 0.067E6 * A_s,  # [nS] leak conductance
                      # [nS] conductance of fast A-type K+ channel
                      'g_K_A': 1E6 * A_s,
                      # [nS] conductance of Na+ dependent K+ channel
                      'g_K_Na': 1.33E6 * A_s,
                      # [nS] conductance of non-inactivating K+ channel
                      'g_K_s': 0.576E6 * A_s,
                      
                      # parameters for dendritic compartment
                      'C_m_d': C_m * A_d,
                      'V_m_d': -75.,
                      # [nS] conductance of Ca+ dependent K+ channel
                      'g_K_Ca': 0.57E6 * A_d,
                      # [nS] conductance of inwardly rectifying channel
                      'g_K_AR': 0.0257E6 * A_d,
                      # [nS] conductance of high-threshold Ca+ channel
                      'g_Ca': 0.43E6 * A_d,
                      # [nS] conductance of Na-persistent channel
                      'g_Na_p': 0.0686E6 * A_d,
    
                      'K_A_h': 0.31,  # initial inactivation value of g_K_A
                      'K_m': 0.015,  # initial activation value of g_K
                      'K_s_m': 0.002,  # initial activation value of g_K_s
                      'Na_h': 1.,  # initial inactivation value of g_Na
                      'n_Ca': 0.,  # [mM] initial concentration of Ca2+
                      'n_Na': 9.5  # [mM] initial concentration of Na+

                 }
               }


# parameters for inh neurons
A_s = 0.02E-2  # [cm^2] somatic area
Params_neuron['FS'] = {
    't_ref': 2.,
    'E_K': -90.,
    'E_L': -63.8,
    'E_Na':  55.,
    'E_ex':   0.,
    'E_in': -70.,
    'tau_syn_ampa': 2.,
    'tau_syn_nmda_fast': 2.,
    'tau_syn_nmda_slow': 100.,
    'tau_syn_gaba': 10.,

    # somatic compartment
    'C_m_s': C_m * A_s,
    'V_m_s': -61.,
    'g_Na': 35E6 * A_s,
    'g_K': 9E6 * A_s,
    'g_L': 0.1025E6 * A_s
}


#####################################################
#########   SIMULATE: CURRENT INJECTION     #########
#####################################################
# For pyramidal neuron type - simulate the current injection
neuron_type = 'PY'

nest.ResetKernel()
nest.SetKernelStatus({"resolution": dt,             
                      "local_num_threads": 4
                     })


# set defaults for nest objects
nest.SetDefaults('dc_generator', Param_dc)

nest.SetDefaults('spike_detector', {
    'start': 0.,
    'stop': sim_time,
    'withgid': True,
    'withtime': True,
    'to_file': False})

nest.SetDefaults('multimeter', {
    'interval': dt,
    'record_from': nest.GetDefaults('compte2003_ex')['recordables'],
    'start': Param_dc['start'] - 100.,  # 100ms before start of DC injection
    'stop': sim_time
})

# create neurons
neuron = nest.Create('compte2003_ex', n=len(Amplitude_dc),
                     params=Params_neuron[neuron_type])


mm = []
sd = []
# for injection amplitude create and connect generator, multimeter and
# spike detector
for dc_id in range(len(Amplitude_dc)):

    # create and connect DC input
    gen = nest.Create('dc_generator',
                      params={'amplitude': Amplitude_dc[dc_id]})
    nest.Connect(gen, [neuron[dc_id]], {'rule': 'one_to_one'},
                 {'receptor_type': syns['curr']})

    # create and connect multimeters
    mm += [nest.Create('multimeter')]
    nest.Connect(mm[-1], [neuron[dc_id]], {'rule': 'one_to_one'})

    # create and connect spike detector
    sd += [nest.Create('spike_detector')]
    nest.Connect([neuron[dc_id]], sd[-1], {'rule': 'one_to_one'})

# simulate
nest.Simulate(sim_time)


#######################################################
########   ANALYSIS     ###############################
#######################################################

Mm_raw = []  # here relevant recordings are stored for each current injection
V_max = []  # [mV] maximal depolarization
I_dc = []   # injected current amplitude
G_origin = []  # [MOhm] membrane conductance, calculated according to
# original model as inverse sum of open conductances
G_th = []  # theoretical estimation of conductance based on classical
# definition

for dc_id in range(len(Amplitude_dc)):

    data = nest.GetStatus(mm[dc_id])[0]['events']

    # many channels reside on dendritic compartment for excitatory neurons,
    # while for inhibitory neurons all channels are in the soma. Here we test
    # excitatory neurons
    comp = '1'

    # reconstruct all neuronal conductances from simulation

    temp = {'V_m': data['V_m_s']}
    temp['V_m_d'] = data['V_m_d']
    temp['times'] = data['times']
    temp['g_ex'] = (data['g_ampa'] + data['g_nmda_slow'] -
                    data['g_nmda_fast'])
    temp['g_in'] = data['g_gaba']
    temp['senders'] = data['senders']

    a_m = 0.1 * (data['V_m_s'] + 33.) / \
        (1. - np.exp(-(data['V_m_s'] + 33.) / 10.))
    b_m = 4. * np.exp(-(data['V_m_s'] + 53.7) / 12.)
    m_inf = a_m / (a_m + b_m)
    temp['g_Na'] = (Params_neuron[neuron_type]['g_Na'] *
                    np.power(m_inf, 3) * data['Na_h'])

    temp['g_K'] = (Params_neuron[neuron_type]['g_K'] *
                   np.power(data['K_m'], 4))

    m_inf = 1. / (1. + np.exp(-(data['V_m_s'] + 50.) / 20.))
    temp['g_K_A'] = (Params_neuron[neuron_type]['g_K_A'] *
                     np.power(m_inf, 3) * data['K_A_h'])

    temp['g_K_s'] = (Params_neuron[neuron_type]['g_K_s'] *
                     data['K_s_m'])

    Na_p_m = 1. / (1. + np.exp(-(data['V_m_d'] + 55.7) / 7.7))
    temp['g_Na_p'] = (Params_neuron[neuron_type]['g_Na_p'] *
                      np.power(Na_p_m, 3))

    K_ar_h = 1. / (1. + np.exp((data['V_m_d'] + 75.) / 4.))
    temp['g_K_ar'] = (Params_neuron[neuron_type]['g_K_AR'] *
                      K_ar_h)

    m_inf = 1. / (1. + np.exp(-(data['V_m_d'] + 20.) / 9.))
    temp['g_Ca'] = (Params_neuron[neuron_type]['g_Ca'] *
                    np.power(m_inf, 2))

    temp['n_Ca'] = data['n_Ca']
    temp['g_K_Ca'] = (1. / (1. + 30. / data['n_Ca']) *
                      Params_neuron[neuron_type]['g_K_Ca'])

    temp['n_Na'] = data['n_Na']
    temp['g_K_Na'] = (0.37 / (1. + np.power(38.7 / data['n_Na'], 3.5)) *
                      Params_neuron[neuron_type]['g_K_Na'])

    Mm_raw += [temp]

    # if dc injection is sub-threshold - extract maximal depolarization 10 ms
    # before the end of dc injection and calculate conductance
    data_spikes = nest.GetStatus(sd[dc_id])[0]['events']
    if len(data_spikes['times']) == 0:
        mask = Mm_raw[-1]['times'] > (Param_dc['stop'] - 10.)
        V_max += [Mm_raw[-1]['V_m'][mask][0] - Mm_raw[-1]['V_m'][0]]
        I_dc += [Amplitude_dc[dc_id]]

        # calculate membrane conductance according to original definition
        # in the model 10ms before the end of dc injection
        g_temp = (Params_neuron[neuron_type]['g_L'] + temp['g_Na']
                  + temp['g_K'] + temp['g_K_A'] +
                  temp['g_K_s'] + temp['g_Na_p']
                  + temp['g_K_ar'] + temp['g_Ca'] + temp['g_K_Ca']
                  + temp['g_K_Na'] + temp['g_ex'] + temp['g_in'])
        G_origin += [g_temp[mask][0]]  # [nS]

        # define the steady state activation of all channels at a given
        # membrane potential. Implemented only for excitatory neurons!

        data = copy.deepcopy(temp)

        # fast Na+ current
        def g_Na(V):
            a = 0.1 * (V + 33.) / (1. - np.exp(-(V + 33.) / 10.))
            b = 4. * np.exp(-(V + 53.7) / 12.)
            m = a / (a + b)

            c = 0.07 * np.exp(-(V + 50.) / 10.)
            d = 1. / (1. + np.exp(-(V + 20.) / 10.))
            h = c / (c + d)
            return (np.power(m, 3) * h *
                    Params_neuron[neuron_type]['g_Na'])

        # fast K+ current
        def g_K(V):
            a = 0.01 * (V + 34.) / (1 - np.exp(-(V + 34.) / 10.))
            b = 0.125 * np.exp(-(V + 44.) / 25.)
            m = a / (a + b)
            return (np.power(m, 4) *Params_neuron[neuron_type]['g_K'])

        # K+ A-current
        def g_K_A(V):
            m = 1. / (1. + np.exp(-(V + 50.) / 20.))
            h = 1. / (1. + np.exp((V + 80.) / 6.))
            return (np.power(m, 3) * h *
                    Params_neuron[neuron_type]['g_K_A'])

        # Na+ persistent current
        def g_Na_p(V):
            m = 1. / (1 + np.exp(-(V + 55.7) / 7.7))
            return (np.power(m, 3) *
                    Params_neuron[neuron_type]['g_Na_p'])

        # K+ slow non-inactivating current
        def g_K_s(V):
            m = 1. / (1. + np.exp(-(V + 34.) / 6.5))
            return m * Params_neuron[neuron_type]['g_K_s']

        # K+ inward-rectifying current
        def g_K_AR(V):
            m = 1. / (1 + np.exp((V + 75.) / 4.))
            return m * Params_neuron[neuron_type]['g_K_AR']

        # Ca2+ current
        def g_Ca(V):
            m = 1. / (1. + np.exp(-(V + 20.) / 9.))
            return m**2 * Params_neuron[neuron_type]['g_Ca']

        #  K+ Ca-dependent current
        def g_K_Ca(V):
            alpha_Ca = (0.005) * 1E-3  # [uM/pA/ms] - 1E-3 is a transition
            # factor from original nA to present pA
            tau_Ca = 150.     # [ms] Ca concentration decay time
            n_Ca = (tau_Ca * alpha_Ca * g_Ca(V) *
                    (V - Params_neuron[neuron_type]['E_Ca']))

            m = n_Ca / (n_Ca + 30.)
            return m * Params_neuron[neuron_type]['g_K_Ca']

        # K+ Na-dependent current
        def g_K_Na(V):
            alpha_Na = (0.01) * 1E-3  # [mM/pA/ms] 1E-3 is a transition factor
            # from original nA to present pA
            R_pump = 0.018         # [mM/ms]
            Na_eq = 9.5

            I = (g_Na(V) + g_Na_p(V)) * \
                (V - Params_neuron[neuron_type]['E_Na'])

            koef = I * alpha_Na / R_pump + Na_eq**3 / (Na_eq**3 + 15**3)
            n_Na = np.power(koef * 15**3 / (1. - koef), 1. / 3)

            m = 0.37 / (1. + np.power(38.7 / n_Na, 3.5))

            return m * Params_neuron[neuron_type]['g_K_Na']

        # total current
        def I_m(V):
            I = ((g_Na(V)+g_Na_p(V)) * (V-Params_neuron[neuron_type]['E_Na']) +
                 
                 (g_K(V) + g_K_A(V) + g_K_s(V) + g_K_AR(V) + g_K_Ca(V) + 
                 g_K_Na(V)) * (V - Params_neuron[neuron_type]['E_K']) +
                 
                 g_Ca(V) * (V - Params_neuron[neuron_type]['E_Ca']) +
                 
                 Params_neuron[neuron_type]['g_L'] * 
                 (V - Params_neuron[neuron_type]['E_L']) +
                 
                 data['g_ex'] * (V - Params_neuron[neuron_type]['E_ex']) +
                 data['g_in'] * (V - Params_neuron[neuron_type]['E_in']))

            return I

        # select the time point just 10 ms before the end of DC injection
        mask = Mm_raw[-1]['times'] > (Param_dc['stop'] - 10.)
        V = Mm_raw[-1]['V_m']
        V_d = Mm_raw[-1]['V_m_d']

        # assume hyperpolarization by 0.5 mV
        dV = -0.5  # [mV[
        G = (I_m(V + dV) - I_m(V)) / (dV)
        G_th += [G[mask][0]]  # [nS]


# calculate differential membrane conductance as a slope of I-V curve.
# Factor 1000 transforms mV/pA into MOhm.
I_dc = np.array(I_dc)
V_max = np.array(V_max)
G_ohm = (I_dc[2:] - I_dc[:-2]) / (V_max[2:] - V_max[:-2])
G_th = np.array(G_th)


##########################################################
############   PLOT   ####################################
##########################################################

# initialize figure
fig = plt.figure()
fig_size = [1, 3]
ax_id = 1
axes = []

# plot membrane potential traces
ax = fig.add_subplot(fig_size[0], fig_size[1], ax_id)
ax_id += 1

for dc_id in range(0, len(Amplitude_dc), 4):
    ax.plot(Mm_raw[dc_id]['times'] / 1000., Mm_raw[dc_id]['V_m'], 'k',
            label='%.1f' % Amplitude_dc[dc_id])
ax.set_xlabel('$\mathrm{time \, (s)}$')
ax.set_ylabel('$V_{\mathrm{m}} \, \mathrm{(mV)}$')
ax.set_ylim(top=-66.)
ax.set_xticks(np.arange(0, 26, 5))
ax.set_xlim([0., sim_time / 1000.])
ax.set_yticks(np.arange(-78., -64., 4.))
ax.vlines([Param_dc['start'] / 1000., Param_dc['stop'] / 1000.],
          ax.get_ylim()[0], ax.get_ylim()[1], color='0.5',
          linestyle='dashed')

axes += [{'ax': ax, 'cb': None}]

# plot I-V curve
ax = fig.add_subplot(fig_size[0], fig_size[1], ax_id)
ax_id += 1
ax.plot(I_dc, V_max, 'k')
ax.set_xlabel('$\mathrm{I \, (pA)}$')
ax.set_xticks(np.linspace(min(I_dc),max(I_dc),3,dtype=int))
ax.set_ylabel('$\Delta V_{\mathrm{m}} \, \mathrm{(mV)}$')
ax.set_yticks([-2.,4.,10.])
axes += [{'ax': ax, 'cb': None}]


# plot differential, Compte's and reconstructed Ohm's conductances vs voltage
ax = fig.add_subplot(fig_size[0], fig_size[1], ax_id)
ax_id += 1
ax.plot(V_max, G_origin, 'k--')
ax.plot(V_max[1:-1], G_ohm, 'k-')
ax.plot(V_max, G_th, 'o', color='0.5', alpha=0.5, markeredgecolor='0.5')
ax.set_ylabel('$G_{\mathrm{m}} \, \mathrm{(nS)}$')
ax.set_xlabel('$\Delta V_{\mathrm{m}} \, \mathrm{(mV)}$')
ax.set_yticks(np.linspace(-5,15,3))
ax.set_xticks(np.linspace(-3,6,4))
# ax.set_yscale('log')
axes += [{'ax': ax, 'cb': None}]


# Finalization phase: adjust subplot dimensions
aw.tight_layout(axes, style, './fig8.png', fig_width='medium')

plt.show()
