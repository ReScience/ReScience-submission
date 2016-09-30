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
# *) Reproduce the AMPA, NMDA, GABA post-synaptic response to pre-synaptic spiking
#    based on the original implementation
# *) Simulate post-synaptic response to the same input with
#    the simplified synaptic model
# *) Compare the performance of the two synaptic implementations
# *) Generate Fig. 2
# -----------------------------------------------------------------------------

import numpy as np
import awesomizer as aw
import matplotlib.pyplot as plt
import copy
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


# parameters for Tsodyks synapse
tsodyks_tau_rec = 130.  # [ms] recovery time constant
tsodyks_u = 0.5  # initial utilization of synaptic vesicles
tsodyks_tau_fac = 0.  # [ms] facilitation time constant


# parameters for DC injection to evoke pre-synaptic spiking
Param_dc = {
    'start': 3000.,  # [ms] start time of DC injection
    'stop': 3500.,  # [ms] stop time of DC injection
    'amplitude': {
        'PY': [50., 150., 250.],
        'FS': [30., 90., 150.]
    }  # [pA] current amplitudes
}


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


#############################################
####    SUPPORT FUNCTIONS  ##################
#############################################


# first-order activation scheme
def f_gate_1order(alpha, tau, f):
    '''
    implements first-order activation scheme with forward Euler for synaptic 
    gating variable s

    alpha, tau - parameters
    f - filtered shape of pre-synaptic membrane potential
    '''
    gate = np.zeros(len(Vm))  # here gating variable is stored
    for time_id in range(1, len(Vm)):

        gate[time_id] = (gate[time_id - 1] +
                         dt * (alpha * f[time_id - 1] - gate[time_id - 1]/tau))

    return gate


# second-order activation scheme

def f_gate_2order(alpha, tau, alpha_x, tau_x, f):
    '''
    implements second-order activation scheme with forward Euler for synaptic 
    gating variable s

    alpha, tau, alpha_x, tau_x - parameters
    f - filtered shape of pre-synaptic membrane potential
    '''

    gate = np.zeros(len(Vm))  # here gating variable is stored
    x = np.zeros(len(Vm))  # here the hidden variable x is stored
    for time_id in range(1, len(Vm)):

        x[time_id] = (x[time_id - 1] +
                      dt * (alpha_x * f[time_id - 1] - x[time_id - 1] / tau_x))

        gate[time_id] = (gate[time_id - 1] + dt *
                         (alpha * (1 - gate[time_id - 1]) * x[time_id - 1] -
                          gate[time_id - 1] / tau))

    return gate, x


###############################################
######   MAIN SIMULATION  #####################
###############################################

# time traces for gating variables
gate_ampa = []
gate_nmda = []
gate_gaba = []

# spike times for each current injection are stored here
spike_times = {'PY': [], 'FS': []}

# amplitudes of [AMPA, NMDA, GABA] gating variable activation
# during first spike in a DC injection window for all DC current amplitudes
gate_ampl = np.zeros([len(Param_dc['amplitude']['PY']), 3])

# for each neuron type, simulate the current injection, extract pre-synaptic
# membrane potential and calculate gating variables. AMPA and NMDA kinetics is
# obtained from pre-synaptic PY neuron, GABA from pre-synaptic FS neuron

for type_id in range(2):

    neuron_type = ['PY', 'FS'][type_id]

    #####################################################
    #########   SIMULATE: CURRENT INJECTION     #########
    #####################################################

    nest.ResetKernel()
    nest.SetKernelStatus({"resolution": dt,             
                          "local_num_threads": 5
                        })
    # create neurons
    neuron_name = 'compte2003_' + ['ex', 'in'][type_id]
    neuron = nest.Create(neuron_name, n=len(Param_dc['amplitude']['PY']),
                         params=Params_neuron[neuron_type])

    # create and connect DC input
    gen = []
    for gen_id in range(len(Param_dc['amplitude']['PY'])):

        params = {'start': Param_dc['start'],
                  'stop': Param_dc['stop'],
                  'amplitude': Param_dc['amplitude'][neuron_type][gen_id]
                  }

        gen += [nest.Create('dc_generator', params=params)]

        nest.Connect(gen[-1], [neuron[gen_id]], {'rule': 'one_to_one'},
                     {'receptor_type': syns['curr']})

    # create and connect multimeter
    params = {
        'interval': dt,
        'record_from': nest.GetDefaults(neuron_name)['recordables'],
    }

    mm = nest.Create('multimeter', n=len(gen), params=params)
    nest.Connect(mm, neuron, {'rule': 'one_to_one'})

    # create and connect spike detector
    sd = nest.Create('spike_detector', n=len(gen))
    nest.Connect(neuron, sd, {'rule': 'one_to_one'})

    # simulate
    nest.Simulate(sim_time)

    #######################################################
    ########   EXTRACT GATING VARIABLES        ############
    #######################################################

    rate = []
    dc_ampl = []
    for gen_id in range(len(gen)):

        spikes = nest.GetStatus([sd[gen_id]])[0]['events']['times']

        if len(spikes) > 0:   # perform analysis only if neuron emitted spikes

            dc_ampl += [Param_dc['amplitude'][neuron_type][gen_id]]
            spike_times[neuron_type] += [spikes]

        # calculate steady firing rate: calculate number of spikes
        # inside DC injection time window. Factor 1000 converts 1/ms to 1/s

        mask = (spikes > Param_dc['start']) * (spikes < Param_dc['stop'])
        rate += [np.sum(mask) / (Param_dc['stop'] - Param_dc['start']) * 1000.]

        # extract the presynaptic spike shape related to first spike in
        # DC injection time window

        # [ms] first 2 spike times
        t0, t1 = spikes[spikes > Param_dc['start']][:2]  
        
        data = nest.GetStatus([mm[gen_id]])[0]['events']

        # mask for data that is around first spike in time
        mask_0 = (data['times'] > (t0 - 3.)) * (data['times'] < t1)

        gate_time = data['times']
        Vm = data['V_m_s']
        
        # filtered pre-synaptic membrane potential
        f = 1. / (1. + np.exp(-(Vm - 20.) / 2.))

        if neuron_type == 'PY':

            # gating variable for AMPA channel
            alpha = 3.48
            tau = Params_neuron[neuron_type]['tau_syn_ampa']
            gate_ampa += [f_gate_1order(alpha, tau, f)]
            gate_ampl[gen_id][0] = np.max(gate_ampa[-1][mask_0])

            # gating variable for NMDA channel
            alpha = 0.5
            tau = Params_neuron[neuron_type]['tau_syn_nmda_slow']
            alpha_x = 3.48
            tau_x = Params_neuron[neuron_type]['tau_syn_nmda_fast']
            gate_nmda += [f_gate_2order(alpha, tau, alpha_x, tau_x, f)[0]]
            gate_ampl[gen_id][1] = np.max(gate_nmda[-1][mask_0])

        else:
            # gating variable for GABA channel
            alpha = 1.
            tau = Params_neuron[neuron_type]['tau_syn_gaba']
            gate_gaba += [f_gate_1order(alpha, tau, f)]
            gate_ampl[gen_id][2] = np.max(gate_gaba[-1][mask_0])

    print ('neuron type:%s ' % neuron_type)
    print ('mean rate (spikes/s): ', rate)
    print ('dc_ampl (pA): ', dc_ampl)


#############################################################
# SIMULATE ORIGINAL AND SIMPLIFIED SYNAPTIC INPUT ###########
# ONTO EXAMPLE PY NEURON                          ###########
#############################################################

# initialize figure
fig = plt.figure()
fig_size = [len(Param_dc['amplitude']['PY']), 3]
axes = []
time_gate = data['times']

# calibration coefficients to match original and simplified
# synaptic implementations
calibr = np.ones([3])


# First perform calibration trial: estimate post-synaptic conductance response
# to injection of input with weight=1 nS to original and simplified
# implementations.
# Second, compare post-synaptic responses in both implementations, with synaptic
# weights adjusted to match amplitude of response to first input.
# Note: while several pre-synaptic firing rates are tested, synaptic weight
# calibration is based on only one injection and then applied to all cases.


for flag_calibr in [True, False]:

    # amplitudes of [AMPA, NMDA, GABA] conductances during first spike
    # in a DC injection window for all DC amplitudes
    g_ampl = np.zeros([len(Param_dc['amplitude']['PY']), 3])

    # for each pre-synaptic firing rate, compare both implementations
    for gen_id in range(len(Param_dc['amplitude']['PY'])):

        nest.ResetKernel()
        nest.SetKernelStatus({"resolution": dt,             
                              "local_num_threads": 5
                            })
        # create 3 neurons for testing AMPA, NMDA and GABA inputs
        neurons = nest.Create('compte2003_ex', n=3, params=Params_neuron['PY'])

        # create 'parrot' neuron for testing Tsodyks synapse.
        # "parrot" neuron is required, because it is not possible to connect
        # spike generator by Tsodyks synapse in NEST
        parrot = nest.Create('parrot_neuron')

        # create input spike generator for AMPA and NMDA channels
        gen = nest.Create('spike_generator',
                          params={'spike_times': spike_times['PY'][gen_id]})

        # connect generator to AMPA channels
        nest.Connect(gen, [neurons[0]], {'rule': 'one_to_one'},
                     {'receptor_type': syns['ampa'],
                      'weight': calibr[0],
                      'delay': 0.1})

        # connect generator to parrot neuron
        nest.Connect(gen, parrot, {'rule': 'one_to_one'})

        # connect parrot neuron to NMDA channels
        nest.SetDefaults('tsodyks2_synapse', params={'tau_rec': tsodyks_tau_rec,
                                                     'U': tsodyks_u,
                                                     'tau_fac': tsodyks_tau_fac})

        nest.Connect(parrot, [neurons[1]], {'rule': 'one_to_one'},
                     {'model': 'tsodyks2_synapse',
                      'receptor_type': syns['nmda_fast'],
                      'weight': calibr[1], 'delay': 0.1})

        # create spike generator for GABA channels
        gen = nest.Create('spike_generator',
                          params={'spike_times': spike_times['FS'][gen_id]})

        # connect generator to GABA channels
        nest.Connect(gen, [neurons[2]], {'rule': 'one_to_one'},
                     {'receptor_type': syns['gaba'],
                      'weight': calibr[2], 'delay': 0.1})

        # create and connect multimeter
        params = {
            'interval': dt,
            'record_from': nest.GetDefaults('compte2003_ex')['recordables'],
        }
        mm = nest.Create('multimeter', n=3, params=params)
        nest.Connect(mm, neurons, {'rule': 'one_to_one'})

        # simulate
        nest.Simulate(sim_time)


        # perform analysis


        # extract amplitude of simulated AMPA conductance
        data = nest.GetStatus(mm)[0]['events']
        time = data['times']
        g_ampa = data['g_ampa']

        spikes = spike_times['PY'][gen_id]
        
        # [ms] first 2 spike times
        t0, t1 = spikes[spikes > Param_dc['start']][:2]  
        
        mask_0 = (time > t0 - 3) * (time < t1)
        g_ampl[gen_id, 0] = np.max(g_ampa[mask_0])

        # extract amplitude of simulated NMDA conductance
        data = nest.GetStatus(mm)[1]['events']
        g_nmda = data['g_nmda_slow'] - data['g_nmda_fast']
        g_ampl[gen_id, 1] = np.max(g_nmda[mask_0])

        # extract amplitude of simulated GABA conductance
        data = nest.GetStatus(mm)[2]['events']
        g_gaba = data['g_gaba']

        spikes = spike_times['FS'][gen_id]
        
        # [ms] first 2 spike times
        t0, t1 = spikes[spikes > Param_dc['start']][:2]  
        
        mask_0 = (time > t0 - 3) * (time < t1)
        g_ampl[gen_id, 2] = np.max(g_gaba[mask_0])

        # when calibration is complete, create figure with both
        # implementations
        if not flag_calibr:

            # plot results for AMPA

            ax = fig.add_subplot(fig_size[0], fig_size[1], 1 + 3 * gen_id)
            # plot gating variable
            ax.plot(time_gate, gate_ampa[gen_id], '0.5')
            # plot recorded conductance
            ax.plot(time, g_ampa, 'k--')
            ax.set_xlabel('$\mathrm{time \, (ms)}$')
            ax.set_ylabel('$g \, (\mathrm{nS})$')
            ax.set_xlim([Param_dc['start'], Param_dc['stop']])
            ax.set_xticks([Param_dc['start'], Param_dc['stop']])
            
            if gen_id == 0:
                ax.set_title('AMPA')

            ax.set_yticks([0., 1.])
            axes += [{'ax': ax, 'cb': None}]

            # plot results for NMDA

            ax = fig.add_subplot(fig_size[0], fig_size[1], 2 + 3 * gen_id)
            # plot gating variable
            ax.plot(time_gate, gate_nmda[gen_id], '0.5')
            # plot recorded conductance
            ax.plot(time, g_nmda, 'k--')
            ax.set_xlabel('$\mathrm{time \, (ms)}$')
            ax.set_ylabel('$g \, (\mathrm{nS})$')
            ax.set_xlim([Param_dc['start'], Param_dc['stop']])
            ax.set_xticks([Param_dc['start'], Param_dc['stop']])
            
            if gen_id == 0:
                ax.set_title('NMDA')

            y_top = np.round(max(np.max(gate_nmda[gen_id]), np.max(g_nmda)), 1)
            ax.set_yticks([0., y_top])
            axes += [{'ax': ax, 'cb': None}]

            # plot results for GABA

            ax = fig.add_subplot(fig_size[0], fig_size[1], 3 + 3 * gen_id)
            # plot gating variable
            ax.plot(time_gate, gate_gaba[gen_id], '0.5')
            # plot recorded conductance
            ax.plot(time, g_gaba, 'k-')
            ax.set_xlabel('$\mathrm{time \, (ms)}$')
            ax.set_ylabel('$g \, (\mathrm{nS})$')
            ax.set_xlim([Param_dc['start'], Param_dc['stop']])
            ax.set_xticks([Param_dc['start'], Param_dc['stop']])
            
            if gen_id == 0:
                ax.set_title('GABA')

            y_top = np.round(max(np.max(gate_gaba[gen_id]), np.max(g_gaba)), 1)
            ax.set_yticks([0., y_top])
            axes += [{'ax': ax, 'cb': None}]

    # calculate calibration coefficients, required to match amplitude of first-input 
    # response in simplified implementation to the one in the original implementation. 
    # Because post-synaptic response is tested on multiple frequencies, 'gen_id' 
    # variable defines the exact post-synaptic conductance trace that is used here.

    if flag_calibr:

        gen_id = 1
        calibr[0] = gate_ampl[gen_id, 0] / g_ampl[gen_id, 0]
        calibr[1] = gate_ampl[gen_id, 1] / g_ampl[gen_id, 1]
        calibr[2] = gate_ampl[gen_id, 2] / g_ampl[gen_id, 2]

# Finalization phase: adjust subplot dimensions
aw.tight_layout(axes, style, './fig2.png', fig_width='medium',
                label_order=['A', 'B', 'C', '', '', '', '', '', ''])

print ('Synaptic calibration:\n', calibr)

plt.show()
