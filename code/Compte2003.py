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
# Simulate reimplemented Compte et al. (2003) model.
# Generates Fig. 4-6, 8
#
# Note the simplified synaptic kinetics, which now uses exponential 
# post-synaptic currents that do not depend on the membrane potential of 
# pre-synaptic neurons. Also, Tsodyks synapses are used to simplify the 
# second-order activation scheme of NMDA channels. 
# -----------------------------------------------------------------------------

import numpy as np
import awesomizer as aw
import matplotlib.pyplot as plt
import copy
import functions as f
import nest


np.random.seed(20)  # set python seed for reproducibility

###########################################
####   PARAMETERS   #######################
###########################################

flag_sim = True   # False - previous simulation results are used and only 
                  # analysis is performed. True - also perform new simulation
sim_time = 9000.  # [ms] simulation time
dt = 0.1          # [ms] simulation time step

num_mm = 20  # number of neurons of each type to record intracellularly   


# network parameters


factor = 1  # coefficient used in simulation to downscale the network, while 
            # preserving local network dynamics. Use to reduce simulation time

N = np.array([1024, 256], dtype=int) / factor  # number of exc and inh neurons

chain_len = 5000./factor  # [um] geometrical length of chain of neurons

sigma_e = 250./chain_len * N  # [units of neuron id] Gaussian connectivity width
sigma_i = 125./chain_len * N  # [units of neuron id] Gaussian connectivity width
sigma = np.transpose(np.array([sigma_e, sigma_i]))

outdegree_mu = 20  # mean number of connections per neuron
outdegree_sd = 5   # standard deviation of number of connections per neuron

#  [nS] synaptic conductances to target populations from [AMPA, NMDA, GABA] inputs
#                 AMPA  NMDA   GABA 
G_syn = np.array([[1.6, 0.24 , 2.  ],  # e        
                  [0.1, 0.024, 0.05]]) # i


# parameters for tsodyks synapse             
tsodyks_tau_rec = 130.  # [ms] recombination time constant
tsodyks_u = 0.5         # initial utilization of synaptic vesicles
tsodyks_tau_fac = 0.    # [ms] facilitation time constant 

syn_delay = dt  #[ms] synaptic delay

# neuronal parameters
syns = nest.GetDefaults('compte2003_ex')['receptor_types']  # receptor IDs                                                                             
C_m = 1E6  # [pF/cm^2] unit membrane capacitance                                                             

# parameters for exc neurons
A_s = 0.015E-2  # [cm^2] somatic area
A_d = 0.035E-2  # [cm^2] dendritic area
Params_neuron = [
    {
      'E_Ca':  120.,  # [mV] Ca2+ reversal potential
      'E_K' : -100.,  # [mV] K+ reversal potential
      'E_L' : -60.95, # [mV] resting potential
      'E_Na':   55.,  # [mV] Na+ reversal potential
      'E_ex':    0.,  # [mV] AMPA, NMDA reversal potential
      'E_in':  -70.,  # [mV] GABA reversal potential
      'g_0' : 1.75E3, # [nS] axial dendritic conductance
      'tau_syn_ampa': 2.,        # [ms] AMPA time constant
      'tau_syn_nmda_fast': 2.,   # [ms] fast NMDA time constant
      'tau_syn_nmda_slow': 100., # [ms] slow NMDA time constant
      'tau_syn_gaba': 10.,       # [ms] GABA time constant
      
      # parameters for somatic compartment
      'comp_0': {   
        'C_m'   : C_m*A_s,    # [pF] capacitance
        'V_m'   : -75.,       # [mV] initial membrane potential   
        'g_Na'  : 50E6*A_s,   # [nS] conductance of fast Na+ channel
        'g_K'   : 10.5E6*A_s, # [nS] conductance of fast K+ channel
        'g_L'   : 0.067E6*A_s,# [nS] leak conductance
        'g_K_A' : 1E6*A_s,    # [nS] conductance of fast A-type K+ channel
        'g_K_Na': 1.33E6*A_s, # [nS] conductance of Na+ dependent K+ channel
        'g_K_s' : 0.576E6*A_s,# [nS] conductance of non-inactivating K+ channel

        'g_K_Ca': 0.,         # [nS] conductance of Ca+ dependent K+ channel
        'g_K_AR': 0.,         # [nS] conductance of inwardly rectifying channel 
        'g_Ca'  : 0.,         # [nS] conductance of high-threshold Ca+ channel
        'g_Na_p': 0.,         # [nS] conductance of Na-persistent channel
        
        'K_A_h': 0.31,        # initial inactivation value of g_K_A       
        'K_m'  : 0.015,       # initial activation value of g_K
        'K_s_m': 0.002,       # initial activation value of g_K_s
        'Na_h' : 1.,          # initial inactivation value of g_Na
        'n_Ca' : 0.,          # [mM] initial concentration of Ca2+
        'n_Na' : 9.5          # [mM] initial concentration of Na+
        
        },
        
      # parameters for dendritic compartment    
      'comp_1': {
        'C_m'   : C_m*A_d,   
        'V_m'   : -75.,
        'g_Na'  : 0.  ,  
        'g_L'   : 0.  ,
        'g_K'   : 0.  ,   
        'g_K_A' : 0.  ,
        'g_K_Na': 0.  ,
        'g_K_s' : 0.  ,   

        'g_K_Ca': 0.57E6*A_d,  
        'g_K_AR': 0.0257E6*A_d,
        'g_Ca'  : 0.43E6*A_d,  
        'g_Na_p': 0.0686E6*A_d,
    
        'K_A_h': 0.31,
        'K_m'  : 0.015,
        'K_s_m': 0.002,
        'Na_h' :  1.,
        'n_Ca' :  0.,    
        'n_Na' : 9.5     
      }
    }
  ]


# parameters for inh neurons 
# axial conductance g_0=0. mimics 1-compartment neuron
A_s = 0.02E-2    # [cm^2] somatic area
Params_neuron += [
      {
      'E_Ca': 120., 
      'E_K' : -90.,  
      'E_L' : -63.8,
      'E_Na':  55., 
      'E_ex':   0., 
      'E_in': -70., 
      'g_0' :   0.,  
      'tau_syn_ampa': 2.,
      'tau_syn_nmda_fast': 2.,
      'tau_syn_nmda_slow': 100.,
      'tau_syn_gaba': 10.,
      
      # somatic compartment
      'comp_0': {
        'C_m'   : C_m*A_s,  
        'V_m'   : -61.,
        'g_Na'  : 35E6*A_s, 
        'g_K'   : 9E6*A_s,  
        'g_L'   : 0.1025E6*A_s,
        'g_K_A' : 0.,
        'g_K_Na': 0.,
        'g_K_s' : 0.,   

        'g_K_Ca': 0.,
        'g_K_AR': 0.,
        'g_Ca'  : 0.,
        'g_Na_p': 0.
        },
    
      # dendritic compartment is irrelevant here
      'comp_1': {
        'C_m'   :100.   ,   
        'V_m'   :-61.   ,
        'g_Na'  :0.     ,  
        'g_L'   :0.     ,
        'g_K'   :0.     ,   
        'g_K_A' :0.     ,
        'g_K_Na':0.     ,
        'g_K_s' :0.     ,   

        'g_K_Ca':0.     ,
        'g_K_AR':0.     ,
        'g_Ca'  :0.     ,
        'g_Na_p':0.     ,
      }   
    
    }]


# some neuronal parameters are Gaussian distributed for exc and inh neurons    
cv_gL = [0.1,0.024]  # sd/mean for leak conductance
sd_EL=[0.3,0.15]     # [mV] sd of leak potential
sd_Gax=0.1E3         # [nS] sd of axial dendritic conductance for exc neurons
    
# Parameters for hyperpolarizing generator (to measure membrane resistance)
# During the whole simulation one neuron is constantly hyperpolarized by 
# injected dc current. With some frequency additional dc pulses are injected 
# to see membrane potential response 
params_hyper = {
    'dur_on': 100. ,  # [ms] duration of hyperpolarizing phase
    'dur_off': 500.,  # [ms] time between hyperpolarizing pulses
    'ampl_base': -250.,  # [pA] current amplitude to keep neuron hyperpolarized   
    'ampl_pulse': -300.  # [pA] current amplitude of injected pulses
    }
    
# spike-detector parameters
Sd_params = {
    'start': 0.,  # [ms] start time of recording
    'stop': sim_time,  # [ms] end time of recording
    'withgid': True, 
    'withtime': True,
    'to_file': False
    }

# multimeter parameters

# list of neuronal variables that can be recorded by multimeter
rqs = nest.GetDefaults('compte2003_ex')['recordables']  

Mm_params = {
    'start': 0.,
    'stop': sim_time,
    'withgid': True, 
    'withtime': True,
    'to_file': False,
    'interval': dt,
    'record_from': rqs
    }

# analysis parameters
bin_width_PH = 10. #[ms] bin width for construction of population time histogram
time_show_spike = [0.,sim_time]  # [ms] time window to construct time histogram
bin_width_syn = 20.  # [ms] width of smoothing rectangular window when comparing 
                     # exc and inh synaptic conductances                         

# if flag_sim==True - perform new simulation          
if flag_sim:
  
    ##################################################
    ########   CREATE CONNECTIVITY MAP   #############
    ##################################################

    def f_gaussian(x, mu, sigma):
      '''
      for geometrical distances x compute Gaussian connection probability 
      x, mu, sigma - in [um]
      '''
      return np.exp(-np.power(x - mu, 2.) / 2. / np.power(sigma, 2.))
    
    
    
    def f_con_1D(t_pop,s_pop):
        '''
        Description:
            -connect target and source populations according to Gaussian rule: 
            each neuron is connected to neighbor neurons to the left and right 
            by a number of synapses, derived from Gaussian law. 
    
        Inputs:	
            -t_pop,s_pop = int - indices of target and source populations

        Return:
            -Connect_map = 2Darray - 2D array with number of synapses 
                                     per each possible connection
        '''
        
        # initialize connection map
        Connect_map = np.zeros([N[t_pop], N[s_pop]], dtype=int)    
        
        # create Gaussian distributed out-degrees for source neurons 
        Outdegree = np.random.normal(outdegree_mu, outdegree_sd, N[s_pop])    
        Outdegree[Outdegree < 1] = 1
        Outdegree = map(int, Outdegree)
  
        # for each source neuron, randomly select target indices 
        target_list = np.arange(0, N[t_pop])
        for s_ind in range(N[s_pop]):   
            
            # create probability weights (normal distribution)
            target_position = 1. * s_ind/N[s_pop] * (N[t_pop]-1)
            Weight_list = f_gaussian(target_list, target_position, 
                                     sigma[t_pop][s_pop])
    
            # remove autapses
            if t_pop==s_pop:
                Weight_list[int(target_position)] = 0.
    
            # randomly select targets according to connection probability 
            target_choice = np.random.choice(target_list, 
                                            size = Outdegree[s_ind],
                                            p = Weight_list/np.sum(Weight_list),
                                            replace = True)
       
            # fill matrix with numbers of synapses per each possible connection
            for t_ind in target_choice:
                Connect_map[t_ind, s_ind]+=1
       
        return Connect_map




    # create connection maps with number of synapses per each possible 
    # connection
    ConMap=[]
    for t_pop in range(2):
        ConMap += [[]]
        for s_pop in range(2):
            ConMap[t_pop] += [f_con_1D(t_pop, s_pop)]


    ###########################################
    ####   CREATE OBJECTS  ####################
    ###########################################

    # initialize NEST kernel
    nest.ResetKernel()
    nest.SetKernelStatus({"resolution": dt,             
                          "local_num_threads": 5
                        })


    # create neurons
         
    Neurons=[]
    for t_pop in range(2):
        neuron_name = 'compte2003_' + ['ex','in'][t_pop]
        nest.SetDefaults(neuron_name,Params_neuron[t_pop])

        # distribute membrane E_L, g_L randomly    
        g_L = (np.random.normal(1., cv_gL[t_pop], N[t_pop]) *
                                Params_neuron[t_pop]['comp_0']['g_L'])
                                
        E_L = (np.random.normal(0., sd_EL[t_pop], N[t_pop]) +
                                Params_neuron[t_pop]['E_L'])
        
        # axial conductance is distributed only for exc neurons
        if t_pop==0:
            g_0 = (np.random.normal(0., sd_Gax, N[t_pop]) +
                   Params_neuron[t_pop]['g_0'])
                   
        else:
            g_0 = np.zeros(N[t_pop])
            
        # create list of dictionaries with distributed neuronal parameters
        params = [{'comp_0': {'g_L': g_L[i]},
                   'E_L'   : E_L[i],
                   'g_0'   : g_0[i]} for i in range(N[t_pop])]

        Neurons += [nest.Create(model=neuron_name, n=N[t_pop], params=params)] 

  
    # create spike detectors
    Sd  = nest.Create('spike_detector', params=Sd_params)   # for exc neurons
    Sd += nest.Create('spike_detector', params=Sd_params)   # for inh neurons

    # create multimeters
    Mm  =[nest.Create('multimeter', n=min(N[0],num_mm), params=Mm_params)]
    Mm +=[nest.Create('multimeter', n=min(N[1],num_mm), params=Mm_params)]
  
    # create basic hyper-polarizing generator 
    gen_basic = nest.Create('dc_generator', params={
                                         'start': 0. ,
                                         'stop' : sim_time,
                                         'amplitude':params_hyper['ampl_base']})     
    
    # create generators for each pulse injection
    t0 = params_hyper['dur_off']
    gen_pulse = []
    while t0 < (sim_time):
        # create generator for base current injection (off period)
        gen_pulse += nest.Create('dc_generator', params={
                                     'start': t0,
                                     'stop': t0+params_hyper['dur_on'],
                                     'amplitude': params_hyper['ampl_pulse']})
                                       

        t0 += params_hyper['dur_off'] + params_hyper['dur_on'] 
      
    ###########################################
    ####   CONNECT OBJECTS  ###################
    ###########################################
                
    # create temporal copies of synaptic models   
    nest.CopyModel('static_synapse','static_temp')
    nest.CopyModel('tsodyks2_synapse','tsodyks_temp')


    for t_pop in xrange(2):      # for each target population
        
        for s_pop in xrange(2):    # for each source population
  
            for s_ind in xrange(N[s_pop]): # for each source neuron   
                
                # connect to those target neurons, where number of synapses per
                # connection is>0
                mask = ConMap[t_pop][s_pop][:,s_ind]>0
                target_gid = list(np.array(Neurons[t_pop])[mask])
              
                if s_pop == 0:  # if source population is excitatory
                    
                    # excitatory inputs project to dendrites for exc, and 
                    # to soma for inh neurons
                    if t_pop == 0:
                        comp = 'C1_'
                    else :
                        comp = 'C0_'

                    # connect AMPA
                    weights = np.transpose( np.array([ G_syn[t_pop, 0] * 
                                          ConMap[t_pop][s_pop][:, s_ind][mask]]))
                    
                    nest.SetDefaults('static_temp', params={
                                     'receptor_type': syns[comp + 'ampa']})
                    
                    nest.Connect([Neurons[s_pop][s_ind]], target_gid,
                                                    {'rule': 'all_to_all'},
                                                    {'model': 'static_temp', 
                                                     'weight': weights, 
                                                     'delay': syn_delay})
          
          
                    # connect NMDA
                    weights = np.transpose( np.array([G_syn[t_pop, 1] * 
                                          ConMap[t_pop][s_pop][:, s_ind][mask]]))
                    
                    nest.SetDefaults('tsodyks_temp', params={
                                        'receptor_type': syns[comp + 'nmda_fast'],
                                        'tau_rec': tsodyks_tau_rec,
                                        'tau_fac': tsodyks_tau_fac,
                                        'U': tsodyks_u})
                    
                    nest.Connect([Neurons[s_pop][s_ind]], target_gid,
                                                    {'rule': 'all_to_all'},
                                                    {'model': 'tsodyks_temp', 
                                                    'weight': weights, 
                                                    'delay': syn_delay})
                else:
      
                    # connect GABA
                    weights = np.transpose( np.array([ G_syn[t_pop,2] * 
                                          ConMap[t_pop][s_pop][:,s_ind][mask]]))
                    
                    nest.SetDefaults('static_temp', params={
                                     'receptor_type': syns['C0_gaba']})
                    
                    nest.Connect([Neurons[s_pop][s_ind]], target_gid,
                                                        {'rule': 'all_to_all'},
                                                        {'model': 'static_temp',
                                                        'weight': weights,
                                                        'delay': syn_delay})
       
        # connect spike detectors
        nest.Connect(Neurons[t_pop], [Sd[t_pop]], {'rule': 'all_to_all'})

    
        # randomly connect multimeters
        n_gids = list(np.random.choice(Neurons[t_pop], len(Mm[t_pop]),
                                                       replace = False))
        
        nest.Connect(Mm[t_pop], n_gids, {'rule':'one_to_one'})
    
        # connect hyperpolarizing generators to first recorded PY neuron
        if t_pop == 0:
            nest.Connect(gen_basic, [n_gids[0]], {'rule':'all_to_all'},
                                      {'receptor_type':syns['C0_curr']})        
            
            nest.Connect(gen_pulse, [n_gids[0]], {'rule':'all_to_all'},
                                      {'receptor_type':syns['C0_curr']})      
  
    #############################################
    #####   SIMULATE  ###########################
    #############################################

    nest.Simulate(sim_time)

  
    ###########################################
    ####   SAVE DATA          #################
    ###########################################

    # transform data to standard form
    Mm_raw = []
    Sd_raw = []
    for t_pop in range(2):
        Mm_raw += [[]]
        mm_id = 0
        for mm in Mm[t_pop]:
            data = nest.GetStatus([mm])[0]['events']

            if t_pop == 0:
                comp = '1'
            else:
                comp = '0'

            temp = {'V_m': data['V_m.0']}  
            temp['times'] = data['times']
            temp['g_ex'] = (data['g_ampa.'+comp] + data['g_nmda_slow.'+comp] -
                            data['g_nmda_fast.'+comp])
            temp['g_in'] = data['g_gaba.0']
            temp['senders'] = data['senders']
      
            a_m = 0.1 *(data['V_m.0']+33.)/(1.-np.exp(-(data['V_m.0']+33.)/10.))
            b_m = 4. * np.exp( -(data['V_m.0']+53.7) / 12.)
            m_inf = a_m / (a_m+b_m)
            temp['g_Na'] = (Params_neuron[t_pop]['comp_0']['g_Na'] * 
                            np.power(m_inf,3)*data['Na_h.0'] )
      
            temp['g_K']  = (Params_neuron[t_pop]['comp_0']['g_K'] * 
                            np.power(data['K_m.0'],4) )
      
            m_inf = 1. / (1. + np.exp( -(data['V_m.0']+50.) / 20.))
            temp['g_K_A'] = (Params_neuron[t_pop]['comp_0']['g_K_A'] *      
                             np.power(m_inf,3) * data['K_A_h.0'])
            
            temp['g_K_s'] = (Params_neuron[t_pop]['comp_0']['g_K_s'] *
                             data['K_s_m.0'] )
      
            Na_p_m = 1. / (1. + np.exp( -(data['V_m.'+comp]+ 55.7) / 7.7))
            temp['g_Na_p'] = (Params_neuron[t_pop]['comp_'+comp]['g_Na_p'] *
                              np.power(Na_p_m,3) )
      
            K_ar_h = 1. / (1. + np.exp((data['V_m.'+comp] + 75.) / 4.))
            temp['g_K_ar'] = K_ar_h*Params_neuron[t_pop]['comp_'+comp]['g_K_AR']      
      
            m_inf = 1. / (1. + np.exp( -(data['V_m.'+comp]+20.) / 9.))
            temp['g_Ca'] = (Params_neuron[t_pop]['comp_'+comp]['g_Ca'] *
                            np.power(m_inf,2) )
      
            temp['n_Ca'] = data['n_Ca.'+comp]
            temp['g_K_Ca'] = (Params_neuron[t_pop]['comp_'+comp]['g_K_Ca'] * 1./
                              (1. + 30./data['n_Ca.'+comp]) )           

            temp['n_Na'] = data['n_Na.0']
            temp['g_K_Na'] = (Params_neuron[t_pop]['comp_0']['g_K_Na'] * 0.37 /
                              (1.+np.power(38.7/data['n_Na.0'],3.5)) )
      
  
            Mm_raw[t_pop] += [temp]
      
            mm_id += 1

        Sd_raw += [nest.GetStatus([Sd[t_pop]])[0]['events']]

    # save data
    np.savez('./Compte2003_mm',Mm_raw)
    np.savez('./Compte2003_sd',Sd_raw)
  
else:
    # if flag_sim==False - no simulation is performed. Read previous results. 
    Mm_raw = np.load('./Compte2003_mm.npz').items()[0][1]
    Sd_raw = np.load('./Compte2003_sd.npz').items()[0][1]










#######################################
#   ANALYSIS ##########################
#######################################


# analysis of spiking activity
 
 
Sd_sorted = []
PH = []
for t_pop in range(2):
    
    # sort population spikes to spike trains of individual neurons
    temp = f.f_Sd_sort(Sd_raw[t_pop], N[t_pop])
    for id_temp in reversed(range(len(temp['ids']))):
        
        # if neuron did not spike, exclude it from consideration
        if temp['ids'][id_temp] == None:
            temp['ids'].pop(id_temp)
            temp['times'].pop(id_temp)

    Sd_sorted += [temp]
    
    # construct population time histogram
    t_range = time_show_spike
    PH_temp = f.f_psth(Sd_sorted[t_pop]['times'], bin_width=bin_width_PH,
                       t_min=t_range[0], t_max=t_range[1])
    
    PH += [{'times': PH_temp['times'], 
            'rates': np.mean(PH_temp['rates'], axis=0)}]
  
  
  
  
  
  

# reconstruct membrane resistance according to virtual current injection method
# for PY neuron. Both instantaneous and steady-state [Na], [Ca] are considered.

t_pop = 0 # further analysis is implemented only for excitatory population
n_id = 2  # ID of neuron to analyze
data = Mm_raw[t_pop][n_id]  # recorded intracellular data

# define the steady state activation of all channels at a given 
# membrane potential
          
                        
# fast Na+ current
def g_Na(V):
     a = 0.1 * (V+33.) / (1.-np.exp(-(V+33.)/10.))
     b = 4. * np.exp(-(V+53.7)/12.)
     m = a/(a+b)
  
     c = 0.07 * np.exp(-(V+50.)/10.);
     d = 1. / (1.+np.exp(-(V+20.)/10.));
     h = c/(c+d)
     return (np.power(m, 3) * h * Params_neuron[t_pop]['comp_0']['g_Na'])

# fast K+ current
def g_K(V):
     a = 0.01*(V+34.)/(1-np.exp(-(V+34.)/10.))
     b = 0.125*np.exp(-(V+44.)/25.)
     m = a/(a+b)
     return (np.power(m, 4) * Params_neuron[t_pop]['comp_0']['g_K'])

# K+ A-current
def g_K_A(V):
     m = 1. / (1.+np.exp(-(V+50.)/20.))
     h = 1. / (1.+np.exp((V+80.)/6.))
     return (np.power(m, 3) * h * Params_neuron[t_pop]['comp_0']['g_K_A'])

# Na+ persistent current
def g_Na_p(V):
     m = 1. / (1+np.exp(-(V+55.7)/7.7))
     return (np.power(m,3) * Params_neuron[t_pop]['comp_1']['g_Na_p'])
  
# K+ slow non-inactivating current
def g_K_s(V):
     m = 1. / (1.+np.exp(-(V+34.)/6.5))
     return m * Params_neuron[t_pop]['comp_0']['g_K_s']
  
# K+ inward-rectifying current
def g_K_AR(V):
     m = 1. / (1+np.exp((V+75.)/4.))
     return m * Params_neuron[t_pop]['comp_1']['g_K_AR']
  
# Ca2+ current 
def g_Ca(V):
     m = 1. / (1.+np.exp(-(V+20.)/9.))
     return m**2 * Params_neuron[t_pop]['comp_1']['g_Ca']


#  K+ Ca-dependent current
def g_K_Ca(V):
    alpha_Ca = (0.005)*1E-3 # [uM/pA/ms] - 1E-3 is a conversion 
                            # factor from original nA to present pA
    tau_Ca = 150.     # [ms] Ca concentration decay time
    n_Ca = (tau_Ca * alpha_Ca*g_Ca(V) * 
                              (V-Params_neuron[t_pop]['E_Ca']) )
                
    m = n_Ca / (n_Ca + 30.)
    return m * Params_neuron[t_pop]['comp_1']['g_K_Ca'] 
            
            
# K+ Na-dependent current
def g_K_Na(V):
    alpha_Na = (0.01)*1E-3 # [mM/pA/ms] 1E-3 is a conversion factor 
                           # from original nA to present pA
    R_pump = 0.018         # [mM/ms]
    Na_eq=9.5
                
    I = (g_Na(V) + g_Na_p(V)) * (V - Params_neuron[t_pop]['E_Na'])
               
    koef = I * alpha_Na / R_pump + Na_eq**3 / (Na_eq**3 + 15**3)
    n_Na = np.power(koef * 15**3 / (1. - koef), 1./3)
               
    m = 0.37 / (1. + np.power(38.7 / n_Na, 3.5))

    return m * Params_neuron[t_pop]['comp_0']['g_K_Na']    

            
# total current with instantaneous [Na] and [Ca]
def I_m_inst(V):
     I = ( (g_Na(V) + g_Na_p(V)) * (V-Params_neuron[t_pop]['E_Na']) +
           
           (g_K(V) + g_K_A(V) + g_K_s(V) + g_K_AR(V) + data['g_K_Ca'] + 
           data['g_K_Na']) * (V-Params_neuron[t_pop]['E_K']) +
           
           g_Ca(V) * (V-Params_neuron[t_pop]['E_Ca']) + 
           
           Params_neuron[t_pop]['comp_0']['g_L'] * 
           (V-Params_neuron[t_pop]['E_L']) +
           
           data['g_ex'] * (V-Params_neuron[t_pop]['E_ex']) + 
           data['g_in'] * (V-Params_neuron[t_pop]['E_in'])  )
  
     return I
     
# total current with steady-state [Na] and [Ca]
def I_m_steady(V):
     I = ( (g_Na(V) + g_Na_p(V)) * (V-Params_neuron[t_pop]['E_Na']) +
           
           (g_K(V) + g_K_A(V) + g_K_s(V) + g_K_AR(V) + g_K_Ca(V) + g_K_Na(V)) * 
           (V-Params_neuron[t_pop]['E_K']) +
           
           g_Ca(V) * (V-Params_neuron[t_pop]['E_Ca']) + 
           
           Params_neuron[t_pop]['comp_0']['g_L'] * 
           (V-Params_neuron[t_pop]['E_L']) +
           
           data['g_ex'] * (V-Params_neuron[t_pop]['E_ex']) + 
           data['g_in'] * (V-Params_neuron[t_pop]['E_in'])  )
  
     return I


# assume hyperpolarization by 0.1 mV
dV = -0.1  # [mV]
V = data['V_m']

G_steady = (I_m_steady(V+dV) - I_m_steady(V)) / dV        
G_inst   = (I_m_inst(V+dV)   - I_m_inst(V)) / dV       
  
# mask out spike events
halfwidth = 20.  # [ms] this time window to left and right from each spike is 
                 # not considered

n_gid = data['senders'][0]
spike_times = Sd_sorted[t_pop]['times'][
                                np.where(Sd_sorted[t_pop]['ids']==n_gid)[0][0] ]




mask = np.ones(len(data['times']), dtype=bool)  
for spike in spike_times:  
  mask_temp = (data['times'] > spike-halfwidth) * \
              (data['times'] < spike+halfwidth)  
  mask[mask_temp] = False
  


time_virtual_masked =  data['times'][mask]  # spike-free time array
G_steady_masked = G_steady[mask]      # spike-free array of virtual conductances
G_inst_masked = G_inst[mask]      # spike-free array of virtual conductances
time_virtual = copy.deepcopy(data['times']) # not masked time array
V_virtual = copy.deepcopy(V)    # not masked V_m array   
spike_times_virtual=copy.deepcopy(spike_times)

#######################################
#   PLOT  SPIKE DATA  #################
####################################### 

# initialization: assign default style format
style = aw.style_Plos()
aw.style_apply(style)

fig = plt.figure(1)
ax_id = 1
fig_size = [3,2]
axes = []


# plot raster plot
ax = fig.add_subplot(fig_size[0], fig_size[1], ax_id)
ax_id += 1

# plot exc and inh populations mixed
for t_pop in [1,0]:
    
    # shift neuronal IDs in order to start from zero
    ids_temp = np.array(Sd_sorted[t_pop]['ids'])
    ids_temp -= min(ids_temp)
    for id_local in xrange(len(ids_temp)):
        # match ID to geometrical position 
        geom_loc = 1. * chain_len * ids_temp[id_local]//N[t_pop]
        
        ax.plot(Sd_sorted[t_pop]['times'][id_local]/1000., geom_loc / 1000. * 
                np.ones(len(Sd_sorted[t_pop]['times'][id_local])),
                'o'+['r','b'][t_pop], markersize=1., markeredgecolor=['r','b'][t_pop])

ax.set_xlabel('time (s)')
ax.set_ylabel('chain position (mm)')
axes += [{'ax':ax, 'cb':None}]

# plot time histogram
ax = fig.add_subplot(fig_size[0], fig_size[1], ax_id)
ax_id += 1

for t_pop in [0,1]:
    ax.plot(PH[t_pop]['times']/1000., PH[t_pop]['rates'], ['r','b'][t_pop])
ax.set_xlabel('time (s)')
ax.set_ylabel('rate (spikes/s)')
axes += [{'ax':ax, 'cb':None}]



for mm_id in [8,2]:

    # plot membrane potential for two neurons: with and without spontaneous 
    # depolarization - choose IDs of multimeters accordingly

    ax = fig.add_subplot(fig_size[0], fig_size[1], ax_id)

  
    data = Mm_raw[0][mm_id]
    ax.plot(data['times']/1000., data['V_m'], 'k')
    ax.set_xlabel('time (s)')
    ax.set_ylim(top = -40.)
    ax.set_ylabel('$V_{\mathrm{m}}$ (mV)')
    axes += [{'ax':ax, 'cb':None}]  
  
    # plot [Na+] for same two neurons
    ax = fig.add_subplot(fig_size[0], fig_size[1], ax_id+2)
    ax_id += 1 
  
    data = Mm_raw[0][mm_id]
    ax.plot(data['times']/1000., data['n_Na'], 'k')
    ax.set_xlabel('time (s)')
    ax.set_ylabel('$[Na^+]$ (mM)')
    axes += [{'ax':ax, 'cb':None}]  

#Finalization phase: adjust subplot dimensions
aw.tight_layout(axes,style, './fig4.png', fig_width='medium',
               label_order=['A','B','C','D','',''])



#######################################
#   PLOT INTRACELLULAR DATA   #########
####################################### 


mm_id = 0
for mm in np.array(Mm_raw[0])[[2]]: # choose range of mm IDs to plot
    
    fig = plt.figure()
    ax_id = 0
    fig_size = [5,2]
    axes = []

    
    
    # plot membrane potential
    
    
    ax = plt.subplot2grid(fig_size, (np.mod(ax_id,5),ax_id/5), colspan=1)
    ax_id += 1

    ax.plot(mm['times']/1000., mm['V_m'], 'k', linewidth=0.5)
    ax.set_xlabel('$\mathrm{time \, (s)}$')
    ax.set_ylabel('$V_{\mathrm{m}}$ $\mathrm{(mV)}$')
    ax.set_ylim(top =- 40.)
    ax.set_yticks([-70., -40.])
    axes += [{'ax':ax, 'cb':None}]  

    
    
    # plot membrane resistance, calculated as inverse sum of conductances
    
    
    ax = plt.subplot2grid(fig_size, (np.mod(ax_id,5),ax_id/5), colspan=1)
    ax_id += 1

    Rm_model = 1000./(mm['g_Na'] + mm['g_K'] + Params_neuron[0]['comp_0']['g_L']  
                     + mm['g_K_A'] + mm['g_K_s'] + mm['g_Na_p'] + mm['g_K_ar'] + 
                     mm['g_Ca'] + mm['g_K_Ca'] + mm['g_K_Na'] + mm['g_ex'] + 
                     mm['g_in'] )
    ax.plot(mm['times']/1000., Rm_model, 'k', linewidth=0.5)
    ax.set_xlabel('$\mathrm{time \, (s)}$')
    ax.set_ylabel('$R_{\mathrm{m}}$ $\mathrm{(M\Omega)}$')
    ax.set_yticks([0., 40., 80.])
    axes += [{'ax':ax, 'cb':None}]  

    
    
    # plot synaptic conductances
    
    
    for syn_id in range(2):
        ax = plt.subplot2grid(fig_size, (np.mod(ax_id,5),ax_id/5), colspan=1)
        ax_id += 1
  
        cond_type = ['ex', 'in'][syn_id] 
        ax.plot(mm['times']/1000., mm['g_'+cond_type], 'k', linewidth=0.5)
        ax.set_xlabel('$\mathrm{time \, (s)}$')
        ax.set_ylabel('$g_{\mathrm{%s}}$ $\mathrm{(nS)}$' %cond_type)
        axes += [{'ax':ax, 'cb':None}]  
    
    
    
    # plot ratio of synaptic conductances
    
    
    ax = plt.subplot2grid(fig_size, (np.mod(ax_id,5),ax_id/5), colspan=1)
    ax_id += 1
  
    # find periods with no spikes
    mask = mm['V_m'] < -60.
    g_in = mm['g_in'][mask]
    g_ex = mm['g_ex'][mask]
  
    # smooth data with rectangular kernel
    num_per_bin = int(bin_width_syn / dt)
    g_in_smooth = [np.mean(g_in[num_per_bin*step : num_per_bin*(step+1)]) 
                  for step in xrange(len(g_in)/num_per_bin)]

    g_ex_smooth = [np.mean(g_ex[num_per_bin*step : num_per_bin*(step+1)]) 
                  for step in xrange(len(g_ex)/num_per_bin)]  

    # maximal value of smoothed conductances
    max_ex = max(g_ex_smooth)
    max_in = max(g_in_smooth)
  
    ax.plot(g_in_smooth, g_ex_smooth, 'o', markersize=1.)

    ax.set_xlabel('$g_{\mathrm{in}} \, \mathrm{(nS)}$')
    ax.set_ylabel('$g_{\mathrm{ex}} \, \mathrm{(nS)}$')
    ax.set_xlim([0., max_in*1.3])
    ax.set_ylim([0., max_ex*1.3])
    axes += [{'ax':ax, 'cb':None}]  
      
 
    # plot conductances
    g_id=0
    for g_type in ['Na_p', 'K_s', 'K_ar', 'K_Ca', 'K_Na']:
        ax=plt.subplot2grid(fig_size, (np.mod(ax_id,5),ax_id/5),colspan=1)
        ax_id+=1
  
        label=['NaP', 'KS', 'AR', 'KCa', 'KNa'][g_id]
        ax.plot(mm['times']/1000., mm['g_'+g_type], 'k', linewidth=0.5)
        ax.set_xlabel('$\mathrm{time \, (s)}$')
        ax.set_ylabel(r'$g_{_{\mathrm{%s}}}$ $\mathrm{(nS)}$' %label)
        ax.set_ylim([0., [2.,2.,7.,1.5,1.5][g_id]])
        axes+=[{'ax':ax, 'cb':None}]  
        g_id+=1

    # Finalization phase: adjust subplot dimensions
    aw.tight_layout(axes,style,'./fig5.png', fig_width='medium',
                    label_order=['A','F','B','G','C','H','D','I','E','J'])
    mm_id += 1

       
    
###################################################
#   COMPARISON OF RESISTANCE MEASUREMENTS  ########
###################################################    
       
mm_id = 0
data = Mm_raw[0][mm_id]

# calculate resistance based on membrane potential deflection in response to 
# hyperpolarizing pulses

Vm = data['V_m']
time = data['times']
Rm_pulse = {'time':[], 'R':[]}
t0 = params_hyper['dur_off']   
while t0 < sim_time:
    ind_base = np.argmin(abs(time-t0))
    V_base = Vm[ind_base]
    t0 += params_hyper['dur_on']   
    Rm_pulse['time'] += [t0]
    
    ind_hyper = np.argmin(abs(time-t0))
    V_hyper = Vm[ind_hyper]
    t0 += params_hyper['dur_off']    
  
    # factor 1000 to transform mV/pA to MOhm
    Rm_pulse['R'] += [(V_hyper-V_base) / params_hyper['ampl_pulse'] *1000.]

Rm_pulse = {'time': np.array(Rm_pulse['time']), 'R': np.array(Rm_pulse['R'])}

# calculate resistance as inverse sum of conductances. 
# Factor 1000 transforms 1/nS to MOhm
Rm_model = 1000. / (Params_neuron[0]['comp_0']['g_L'] + data['g_Na'] +   
           data['g_K'] + data['g_K_A'] + data['g_K_s'] + data['g_Na_p'] + 
           data['g_K_ar'] + data['g_Ca'] + data['g_K_Ca'] + data['g_K_Na'] + 
           data['g_ex'] + data['g_in'])




fig = plt.figure()
ax_id = 1
fig_size = [1,2]
axes = []


# plot membrane potential trace


ax = fig.add_subplot(fig_size[0],fig_size[1],ax_id)
ax_id += 1

ax.plot(data['times']/1000., data['V_m'], 'k', linewidth=0.5) 
ax.set_xlabel('$\mathrm{time \, (s)}$')
ax.set_ylabel('$V_{\mathrm{m}}$ $\mathrm{(mV)}$')
ax.set_ylim(top=-70.)
ax.set_yticks([-100., -70.])
ax.set_xlim(right=9.1)
axes += [{'ax':ax, 'cb':None}]  



# plot measured resistance


ax = fig.add_subplot(fig_size[0], fig_size[1], ax_id)
ax_id += 1

ax.plot(data['times']/1000., Rm_model, '0.5')   # 1000 transforms ms->s
ax.plot(Rm_pulse['time']/1000., Rm_pulse['R'], 'ok', markersize=5) 
ax.set_xlabel('$\mathrm{time \, (s)}$')
ax.set_ylabel('$R_{\mathrm{m}} \, \mathrm{(M\Omega)}$')
ax.set_xlim(right=9.1)
axes += [{'ax':ax, 'cb':None}] 

# Finalization phase: adjust subplot dimensions
aw.tight_layout(axes, style, './fig6.png',
                fig_width='medium', label_order=['A','B'])
mm_id+=1




############################################################################
#   PLOT CONDUCTANCE MEASURED WITH VIRTUAL HYPERPOLARIZATION METHOD    #####
############################################################################


fig=plt.figure()
fig_size=[1,3]
ax_id=1
axes=[]


# plot membrane potential


ax = fig.add_subplot(fig_size[0],fig_size[1],ax_id)
ax_id += 1

ax.plot(time_virtual/1000., V_virtual, 'k',)
ax.set_xlabel('$\mathrm{time \, (s)}$')
ax.set_ylabel('$V_{\mathrm{m}}$ \, $\mathrm{(mV)}$')
ax.set_yticks(np.linspace(-80.,40.,4))
ax.set_xlim([0., 6.])
axes += [{'ax':ax, 'cb':None}]


# plot conductances for steady-state and instantaneous [Na] and [Ca]
ind = 0
for  g in [G_inst_masked,G_steady_masked]:
    ax = fig.add_subplot(fig_size[0], fig_size[1], ax_id)
    ax_id += 1

    ax.plot(time_virtual_masked/1000., g, 'ko', markersize=0.8)

    # add info about spike times
    for time in spike_times_virtual:
        ax.bar(left = (time-halfwidth)/1000.,
               width = halfwidth*2/1000., 
               bottom = min(g), 
               height = max(g) - min(g),
               linewidth=0., color='0.5', alpha=0.3)
                         
                         
    # add zero line
    ax.hlines(0., ax.get_ylim()[0], ax.get_ylim()[1],
                  linestyle = 'dashed', color='0.2', alpha=0.5)

    ax.set_xlabel('$\mathrm{time \, (s)}$')
    ax.set_ylabel('$G_{\mathrm{m}} \, \mathrm{(nS)}$' )
    ax.set_yticks(np.linspace(-30.,15.,4))
    ax.set_xlim([0., 6.])
    axes+=[{'ax':ax,'cb':None}]
    ind+=1


# Finalization phase: adjust subplot dimensions
aw.tight_layout(axes, style, './fig8.png',
                fig_width='medium', label_order=['A','B','C'])



plt.show()

