# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Andrei Maksimov (maksimov.andrei7@gmail.com)
# -----------------------------------------------------------------------------
# Reference:
#
# Cellular and network mechanisms of slow oscillatory activity (<1 Hz)
# and wave propagations in a cortical network model*, A. Compte,
# M.V. Sanchez-Vives, D.A. McCormick, X.-J. Wang,
# Journal of Neurophysiology, 2707--2725, 2003"
# -----------------------------------------------------------------------------
# File description:
#
# Simulate response of neurons to DC injection and generate Fig. 1
# -----------------------------------------------------------------------------

import numpy as np
import awesomizer as aw
import matplotlib.pyplot as plt
import nest


# initialize NEST kernel
nest.ResetKernel()

# initialize matplotlib defaults
style = aw.style_Plos()
aw.style_apply(style)


###################################################
###############  PARAMETERS  ######################
###################################################

sim_time = 5000.  # [ms] simulation time
dt = 0.1  # [ms] simulation time step


# parameters for dc current injection
Param_dc = {
    'start': 3000.,  # [ms] start time of DC injection
    'stop': 3500.,   # [ms] stop time of DC injection
    'amplitude': 250.  # [pA] current amplitude
}


syns = nest.GetDefaults('compte2003_ex')['receptor_types']  # receptor IDs

C_m = 1E6  # [pF/cm^2] unit membrane capacitance

# parameters for exc neurons
A_s = 0.015E-2  # [cm^2] somatic area
A_d = 0.035E-2  # [cm^2] dendritic area
Params_neuron = {'PY':
                 {
                     'E_Ca':  120.,  # [mV] Ca2+ reversal potential
                     'E_K': -100.,   # [mV] K+ reversal potential
                     'E_L': -60.95,  # [mV] resting potential
                     'E_Na':   55.,  # [mV] Na+ reversal potential
                     'E_ex':    0.,  # [mV] AMPA, NMDA reversal potential
                     'E_in': -70.,   # [mV] GABA reversal potential
                     'g_0': 1.75E3,  # [nS] axial dendritic conductance
                     'tau_syn_ampa': 2.,  # [ms] AMPA time constant
                     'tau_syn_nmda_fast': 2.,  # [ms] fast NMDA time constant
                     'tau_syn_nmda_slow': 100.,  # [ms] slow NMDA time constant
                     'tau_syn_gaba': 10.,  # [ms] GABA time constant

                     # parameters for somatic compartment
                     'comp_0': {
                         'C_m': C_m * A_s,  # [pF] capacitance
                         'V_m': -75.,  # [mV] initial membrane potential
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

                         # [nS] conductance of Ca+ dependent K+ channel
                         'g_K_Ca': 0.,
                         # [nS] conductance of inwardly rectifying channel
                         'g_K_AR': 0.,
                         # [nS] conductance of high-threshold Ca+ channel
                         'g_Ca': 0.,  
                         # [nS] conductance of Na-persistent channel
                         'g_Na_p': 0.,

                         'K_A_h': 0.31,  # initial inactivation value of g_K_A
                         'K_m': 0.015,  # initial activation value of g_K
                         'K_s_m': 0.002,  # initial activation value of g_K_s
                         'Na_h': 1.,  # initial inactivation value of g_Na
                         'n_Ca': 0.,  # [mM] initial concentration of Ca2+
                         'n_Na': 9.5  # [mM] initial concentration of Na+

                     },

                     # parameters for dendritic compartment
                     'comp_1': {
                         'C_m': C_m * A_d,
                         'V_m': -75.,
                         'g_Na': 0.,
                         'g_L': 0.,
                         'g_K': 0.,
                         'g_K_A': 0.,
                         'g_K_Na': 0.,
                         'g_K_s': 0.,

                         'g_K_Ca': 0.57E6 * A_d,
                         'g_K_AR': 0.0257E6 * A_d,
                         'g_Ca': 0.43E6 * A_d,
                         'g_Na_p': 0.0686E6 * A_d,

                         'K_A_h': 0.31,
                         'K_m': 0.015,
                         'K_s_m': 0.002,
                         'Na_h':  1.,
                         'n_Ca':  0.,
                         'n_Na': 9.5
                     }
                 }
                 }


# parameters for inh neurons
# axial conductance g_0=0. mimics 1-compartment neuron
A_s = 0.02E-2  # [cm^2] somatic area
Params_neuron['FS'] = {
    'E_Ca': 120.,
    'E_K': -90.,
    'E_L': -63.8,
    'E_Na':  55.,
    'E_ex':   0.,
    'E_in': -70.,
    'g_0':   0.,
    'tau_syn_ampa': 2.,
    'tau_syn_nmda_fast': 2.,
    'tau_syn_nmda_slow': 100.,
    'tau_syn_gaba': 10.,

    # somatic compartment
    'comp_0': {
        'C_m': C_m * A_s,
        'V_m': -61.,
        'g_Na': 35E6 * A_s,
        'g_K': 9E6 * A_s,
        'g_L': 0.1025E6 * A_s,
        'g_K_A': 0.,
        'g_K_Na': 0.,
        'g_K_s': 0.,

        'g_K_Ca': 0.,
        'g_K_AR': 0.,
        'g_Ca': 0.,
        'g_Na_p': 0.
    },

    # dendritic compartment is irrelevant here
    'comp_1': {
        'C_m': 100.,
        'V_m': -61.,
        'g_Na': 0.,
        'g_L': 0.,
        'g_K': 0.,
        'g_K_A': 0.,
        'g_K_Na': 0.,
        'g_K_s': 0.,

        'g_K_Ca': 0.,
        'g_K_AR': 0.,
        'g_Ca': 0.,
        'g_Na_p': 0.,
    }

}


# initialize figure for DC injection results
fig = plt.figure(0)
fig_size = [1, 2]
axes = []

# for each neuron type - simulate the current injection
for type_id in range(2):

    neuron_type = ['PY', 'FS'][type_id]

    #####################################################
    #########   SIMULATE: CURRENT INJECTION     #########
    #####################################################

    # initialize NEST kernel
    nest.ResetKernel()

    # create neuron
    neuron_name = 'compte2003_' + ['ex', 'in'][type_id]
    neuron = nest.Create(neuron_name, params=Params_neuron[neuron_type])

    # create and connect DC input
    gen = nest.Create('dc_generator', params=Param_dc)
    nest.Connect(gen, neuron, {'rule': 'one_to_one'},
                 {'receptor_type': syns['C0_curr']})

    # create and connect multimeter
    params = {
        'interval': dt,
        'record_from': nest.GetDefaults('compte2003_ex')['recordables'],
    }
    mm = nest.Create('multimeter', params=params)
    nest.Connect(mm, neuron, {'rule': 'one_to_one'})

    # create and connect spike detector
    sd = nest.Create('spike_detector')
    nest.Connect(neuron, sd)

    # simulate
    nest.Simulate(sim_time)

    #######################################################
    ########   ANALYSIS: CURRENT INJECTION     ############
    #######################################################

    # calculate firing rate during DC injection.
    spikes = nest.GetStatus(sd)[0]['events']['times']
    mask = (spikes > Param_dc['start']) * (spikes < Param_dc['stop'])
    # Factor 1000 converts 1/ms to 1/s
    rate = np.sum(mask) / (Param_dc['stop'] - Param_dc['start']) * 1000.

    print '%s: spike rate=%.2f' % (neuron_type, rate)

    ##########################################################
    ##############    PLOT: CURRENT INJECTION    #############
    ##########################################################

    ax = fig.add_subplot(fig_size[0], fig_size[1], type_id + 1)

    # get data from multimeter
    data = nest.GetStatus(mm)[0]['events']

    # plot membrane potential trace
    ax.plot(data['times'], data['V_m.0'], 'k', label='C0')
    ax.set_xlabel('$\mathrm{time} \, \mathrm{(ms)}$')
    ax.set_ylabel(r'$V_{\mathrm{m}} \, \mathrm{(mV)}$')
    ax.set_xlim([Param_dc['start'] - 100., Param_dc['stop'] + 100.])
    ax.set_xticks([Param_dc['start'], Param_dc['stop']])
    ax.set_ylim([-80.,-40.])
    ax.set_yticks(np.linspace(-80.,-40.,3))
    axes += [{'ax': ax, 'cb': None}]

# Finalization phase: adjust subplot dimensions
aw.tight_layout(axes, style, './fig1.png', fig_width='medium')

plt.show()
