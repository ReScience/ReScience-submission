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
# Generates Fig. 5
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

np.random.seed(65)  # set python seed for reproducibility

###########################################
####   PARAMETERS   #######################
###########################################

flag_sim = True # False - previous simulation results are used and only 
                # analysis is performed. True - also perform new simulation
sim_time = 6200. # [ms] simulation time
dt = 0.1         # [ms] simulation time step

num_mm = 20  # number of neurons of each type to record intracellularly   


# network parameters


factor =  1 # coefficient used in simulation to downscale the network, while 
            # preserving local network dynamics. Use to reduce simulation time

N = np.array([1024/factor, 256/factor], dtype=int)   # number of exc and 
                                                     # inh neurons

chain_len = 5000./factor  # [um] geometrical length of chain of neurons

sigma_e = 250./chain_len * N  # [units of neuron id] Gaussian connectivity width
sigma_i = 125./chain_len * N  # [units of neuron id] Gaussian connectivity width
sigma = np.transpose(np.array([sigma_e, sigma_i]))

outdegree_mu = 20  # mean number of connections per neuron
outdegree_sd = 5   # standard deviation of number of connections per neuron

#  [nS] synaptic conductances to target populations from [AMPA, NMDA, GABA] inputs
#                 AMPA  NMDA   GABA                                                   
G_syn = np.array([[7., 0.15 , 16.  ],  # e        
                  [3., 0.  ,  2.   ]])  # i
                                   

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
                     't_ref': 1.,    # [ms] absolute refractory time 
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
               ]
               
               



# parameters for inh neurons
A_s = 0.02E-2  # [cm^2] somatic area
Params_neuron += [{
    't_ref': 1.,
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
}]


# some neuronal parameters are Gaussian distributed for exc and inh neurons    
cv_gL = [0.1,0.024]  # sd/mean for leak conductance
sd_EL = [0.3,0.15]   # [mV] sd of leak potential
sd_Gax= 0.1E3        # [nS] sd of axial dendritic conductance for exc neurons


# parameters for generators used to damp spontaneous activity and then to 
# initiate network activity propagation. Hyperpolarization is applied to all 
# excitatory neurons. Further stimulation is applied to two groups of neurons at
# both sides of the chain.  
params_stim = {
    'start_hyper': 0000. ,  # [ms] start time of hyperpolarization
    'stop_hyper':  18000.,  # [ms] stop time of hyperpolarization
    'ampl_hyper': -6.,  # [pA] current amplitude to keep neuron 
                               # hyperpolarized   
    'start_pulse': [2000.,6000.] ,  # [ms] start times of stimulation periods
    'duration':  50.,  # [ms] duration of stimulation
    'ampl_stim': 200., # [pA] current amplitude
    'num_neuron': [40,0] # number of neurons to stimulate for exc and inh 
                          # populations, starting from neuronal index N/2
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
time_show_PH = [2000.,5000.] # [ms] time window used to construct population 
                             # time histogram                    

# if flag_sim==True - perform new simulation          
if flag_sim:
  
    ##################################################
    ########   CREATE CONNECTIVITY MAPS  #############
    ##################################################

    def f_gaussian(x, mu, sigma):
        '''
        for geometrical distances x compute Gaussian connection probability 
        x, mu, and sigma should have the same units
        '''
        return np.exp(-np.power(x - mu, 2.) / 2. / np.power(sigma, 2.))

 
    
    # initialize connection maps
    ConMap=[]
    for t_pop in range(2):
        ConMap += [[]]
        for s_pop in range(2):
            ConMap[t_pop] += [np.zeros([N[t_pop], N[s_pop]])]


    ind_to_con_ex = np.arange(0, N[0], 1) # indices of neurons in exc population
    ind_to_con_in = np.arange(0, N[1], 1) # indices of neurons in inh population
    
    for s_pop in range(2):
            
        # create Gaussian distributed out-degrees for source neurons 
        Outdegree = np.random.normal(outdegree_mu, outdegree_sd, N[s_pop])    
        Outdegree[Outdegree < 1] = 1
        Outdegree = np.array(Outdegree, dtype=int)
        
        for s_ind in range(N[s_pop]):   
            
            # create probability weights for target excitatory population
            target_position = 1. * s_ind / (N[s_pop]-1) * (N[0]-1)
            Weight_list_ex = f_gaussian(ind_to_con_ex, target_position, 
                                       sigma[0][s_pop])
                                     
                                     
            # create probability weights for target inhibitory population
            target_position = 1. * s_ind/N[s_pop] * (N[1]-1)
            Weight_list_in = f_gaussian(ind_to_con_in, target_position, 
                                       sigma[1][s_pop])
            
            # create lists of weights and ids for merged populations                           
            Weight_list = np.append(Weight_list_ex, Weight_list_in)
            ind_to_con = np.arange(0,np.sum(N),1)
                        
            # remove autapses
            s_gid = s_ind + N[0]*s_pop
            Weight_list = np.delete(Weight_list, s_gid)
            ind_to_con = np.delete(ind_to_con, s_gid)
                
            # randomly select targets from both populations according to 
            # connection probability 
            target_choice = np.random.choice(ind_to_con, 
                                            size = Outdegree[s_ind],
                                            p = Weight_list/np.sum(Weight_list),
                                            replace = True)
            
            # split targets to individual populations
            target_ids = [target_choice[target_choice < N[0]]]  
            target_ids += [target_choice[target_choice >= N[0]]]
            
            target_ids[1] -= N[0] # to make ids to start from zero           
            
            # fill connectivity maps with numbers of synapses for chosen 
            # connections
            for t_pop in range(2):
              
                Con_local=ConMap[t_pop][s_pop]
                for t_ind in target_ids[t_pop]:
                    Con_local[t_ind, s_ind] += 1
                
                ConMap[t_pop][s_pop] = Con_local
   
      
        
    ###########################################
    ####   CREATE OBJECTS  ####################
    ###########################################

    # initialize NEST kernel
    nest.ResetKernel()
    nest.SetKernelStatus({"resolution": dt,             
                          "local_num_threads": 4
                        })


    # create neurons
         
    Neurons=[]
    for t_pop in range(2):
        neuron_name = 'compte2003_' + ['ex','in'][t_pop]
        nest.SetDefaults(neuron_name,Params_neuron[t_pop])

        # distribute membrane E_L, g_L randomly    
        g_L = (np.random.normal(1., cv_gL[t_pop], N[t_pop]) *
                                Params_neuron[t_pop]['g_L'])
                                
        E_L = (np.random.normal(0., sd_EL[t_pop], N[t_pop]) +
                                Params_neuron[t_pop]['E_L'])
                               
        
        # axial conductance is distributed only for exc neurons
        if t_pop==0:
            g_conn = (np.random.normal(0., sd_Gax, N[t_pop]) +
                   Params_neuron[t_pop]['g_conn'])
                               
            # create list of dictionaries with distributed neuronal parameters
            params = [{'g_L': g_L[i],
                       'E_L': E_L[i],
                       #'V_m_s': E_L[i],
                       'g_conn'   : g_conn[i]} for i in range(N[t_pop])]
        else:
                               
            # create list of dictionaries with distributed neuronal parameters
            params = [{'g_L': g_L[i],
                       'E_L'   : E_L[i]} for i in range(N[t_pop])]          

        Neurons += [nest.Create(model=neuron_name, n=N[t_pop], params=params)] 

  
    # create spike detectors
    Sd  = nest.Create('spike_detector', params=Sd_params)   # for exc neurons
    Sd += nest.Create('spike_detector', params=Sd_params)   # for inh neurons

    # create multimeters
    Mm  =[nest.Create('multimeter', n=min(N[0],num_mm), params=Mm_params)]  
 
 
    # create hyper-polarizing generator to damp spontaneous activity
    gen_stim_basic = nest.Create('dc_generator', params={
                                   'start': params_stim['start_hyper'] ,
                                   'stop' : params_stim['stop_hyper'],
                                   'amplitude': params_stim['ampl_hyper']})     
    
    # create generators for each pulse injection while stimulating network
    gen_stim_pulse = []
    for t_stim in params_stim['start_pulse']:
    
        gen_stim_pulse += nest.Create('dc_generator', params={
                                     'start': t_stim,
                                     'stop': t_stim + params_stim['duration'],
                                     'amplitude': params_stim['ampl_stim']})  
        
    ###########################################
    ####   CONNECT OBJECTS  ###################
    ###########################################
                
    # create temporal copies of synaptic models   
    nest.CopyModel('static_synapse','static_temp')
    nest.CopyModel('tsodyks2_synapse','tsodyks_temp')


    for t_pop in range(2):      # for each target population
        
        for s_pop in range(2):    # for each source population
  
            for s_ind in range(N[s_pop]): # for each source neuron   
                
                # connect to those target neurons, where number of synapses per
                # connection is>0
                mask = ConMap[t_pop][s_pop][:,s_ind]>0
                target_gid = list(np.array(Neurons[t_pop])[mask])
              
                if s_pop == 0:  # if source population is excitatory
                    

                    # connect AMPA
                    weights = np.transpose( np.array([ G_syn[t_pop, 0] * 
                                          ConMap[t_pop][s_pop][:, s_ind][mask]]))
                    
                    nest.SetDefaults('static_temp', params={
                                     'receptor_type': syns['ampa']})
                    
                    nest.Connect([Neurons[s_pop][s_ind]], target_gid,
                                                    {'rule': 'all_to_all'},
                                                    {'model': 'static_temp', 
                                                     'weight': weights, 
                                                     'delay': syn_delay})
          
          
                    # connect NMDA
                    weights = np.transpose( np.array([G_syn[t_pop, 1] * 
                                          ConMap[t_pop][s_pop][:, s_ind][mask]]))
                    
                    nest.SetDefaults('tsodyks_temp', params={
                                        'receptor_type': syns['nmda_fast'],
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
                                     'receptor_type': syns['gaba']})
                    
                    nest.Connect([Neurons[s_pop][s_ind]], target_gid,
                                                        {'rule': 'all_to_all'},
                                                        {'model': 'static_temp',
                                                        'weight': weights,
                                                        'delay': syn_delay})
       
        # connect spike detectors
        nest.Connect(Neurons[t_pop], [Sd[t_pop]], {'rule': 'all_to_all'})

    
        # randomly connect multimeters
        if t_pop==0:
          n_gids = list(np.random.choice(Neurons[t_pop], len(Mm[t_pop]),
                                                       replace = False))
        
          nest.Connect(Mm[t_pop], n_gids, {'rule':'one_to_one'})


        # connect generators to damp spontaneous activity
        if t_pop == 0:
            nest.Connect(gen_stim_basic, Neurons[t_pop], 
                                         {'rule':'all_to_all'},
                                         {'receptor_type':syns['curr']})         
       
        # connect generators to stimulate network
        num_to_stim = params_stim['num_neuron']                              
        nest.Connect(gen_stim_pulse, 
                     Neurons[t_pop][N[t_pop]//2 : N[t_pop]//2+num_to_stim[t_pop]], 
                     {'rule':'all_to_all'},
                     {'receptor_type':syns['curr']})   
    

    #############################################
    #####   SIMULATE  ###########################
    #############################################

    nest.Simulate(sim_time)

    ###########################################
    ####   SAVE DATA          #################
    ###########################################

    # prepare multimeter data for saving
    Mm_raw = []
    for t_pop in [0]:
        Mm_raw += [[]]
        mm_id = 0
        for mm in Mm[t_pop]:
            data = nest.GetStatus([mm])[0]['events']

            temp = {'V_m': data['V_m_s']}  
            temp['times'] = data['times']
            temp['g_ex'] = (data['g_ampa'] + data['g_nmda_slow'] -
                            data['g_nmda_fast'])
            temp['g_in'] = data['g_gaba']
            temp['senders'] = data['senders']
      
            a_m = 0.1 *(data['V_m_s']+33.)/(1.-np.exp(-(data['V_m_s']+33.)/10.))
            b_m = 4. * np.exp( -(data['V_m_s']+53.7) / 12.)
            m_inf = a_m / (a_m+b_m)
            temp['g_Na'] = (Params_neuron[t_pop]['g_Na'] * 
                            np.power(m_inf,3)*data['Na_h'] )
      
            temp['g_K']  = (Params_neuron[t_pop]['g_K'] * 
                            np.power(data['K_m'],4) )
      
            m_inf = 1. / (1. + np.exp( -(data['V_m_s']+50.) / 20.))
            temp['g_K_A'] = (Params_neuron[t_pop]['g_K_A'] *      
                             np.power(m_inf,3) * data['K_A_h'])
            
            temp['g_K_s'] = (Params_neuron[t_pop]['g_K_s'] *
                             data['K_s_m'] )
      
            Na_p_m = 1. / (1. + np.exp( -(data['V_m_d']+ 55.7) / 7.7))
            temp['g_Na_p'] = (Params_neuron[t_pop]['g_Na_p'] *
                              np.power(Na_p_m,3) )
      
            K_ar_h = 1. / (1. + np.exp((data['V_m_d'] + 75.) / 4.))
            temp['g_K_ar'] = K_ar_h*Params_neuron[t_pop]['g_K_AR']      
      
            m_inf = 1. / (1. + np.exp( -(data['V_m_d']+20.) / 9.))
            temp['g_Ca'] = (Params_neuron[t_pop]['g_Ca'] *
                            np.power(m_inf,2) )
      
            temp['n_Ca'] = data['n_Ca']
            temp['g_K_Ca'] = (Params_neuron[t_pop]['g_K_Ca'] * 1./
                              (1. + 30./data['n_Ca']) )           

            temp['n_Na'] = data['n_Na']
            temp['g_K_Na'] = (Params_neuron[t_pop]['g_K_Na'] * 0.37 /
                              (1.+np.power(38.7/data['n_Na'],3.5)) )
      
  
            Mm_raw[t_pop] += [temp]
      
            mm_id += 1
    
    # prepare sd data for saving
    Sd_raw = []
    for t_pop in range(2):
        Sd_raw += [nest.GetStatus([Sd[t_pop]])[0]['events']]

    # save data
    np.savez('./Compte2003_stim_mm',Mm_raw)
    np.savez('./Compte2003_stim_sd',Sd_raw)
  
else:
    # if flag_sim==False - no simulation is performed. Read previous results. 
    Mm_raw = np.load('./Compte2003_stim_mm.npz').items()[0][1]
    Sd_raw = np.load('./Compte2003_stim_sd.npz').items()[0][1]


#######################################
#   ANALYSIS ##########################
#######################################

# analysis of spiking activity
  

# sort population spikes to spike trains of individual neurons
Sd_sorted = []
for t_pop in range(2):
    
    temp = f.f_Sd_sort(Sd_raw[t_pop], N[t_pop])
      
    for id_temp in reversed(range(len(temp['ids']))):
        
        # if neuron did not spike, exclude it from consideration
        if temp['ids'][id_temp] == None:
            temp['ids'].pop(id_temp)
            temp['times'].pop(id_temp)

    Sd_sorted += [temp]
    
    
# sort population spikes inside the given window to spike trains of 
# individual neurons
Sd_sorted_lim = []
for t_pop in range(2):
    # select spike data from the given time window
    mask = ((Sd_raw[t_pop]['times'] > time_show_PH[0]) * 
            (Sd_raw[t_pop]['times'] < time_show_PH[1]))
    Sd_lim = {'senders':Sd_raw[t_pop]['senders'][mask],
              'times':Sd_raw[t_pop]['times'][mask]}
    
    # sort data to spike trains
    temp = f.f_Sd_sort(Sd_lim, N[t_pop])
    
    # if neuron did not spike, exclude it from consideration
    for id_temp in reversed(range(len(temp['ids']))):
        
        if temp['ids'][id_temp] == None:
            temp['ids'].pop(id_temp)
            temp['times'].pop(id_temp)
    
    Sd_sorted_lim += [temp]




# construct population time histogram. Spike times of neurons are calculated 
# relative to the first spike time of their neighbours of opposite type

PH = []
Lag = []
gid_0 = [min(Sd_sorted_lim[0]['ids'])]  # initial GID of exc neurons
gid_0 += [min(Sd_sorted_lim[1]['ids'])] # initial GID of inh neurons
for t_pop in range(2):
    
    
    # correct spike trains for the adjacent neuron first spike time and 
    # construct population time histogram
   
    # calculate difference of the first neuronal spike time relative to first 
    # spike time of adjacent neuron of opposite type 
   
    trains_corrected = []
    Lag += [[]]
    t_min=1E6; t_max=0. # initiate minimal and maximal values of relative spike times
    for n_id in range(1,len(Sd_sorted_lim[t_pop]['ids'])-1):
        
        # calculate neuron index, if given neuron would belong to opposite type
        gid = Sd_sorted_lim[t_pop]['ids'][n_id] # GID of given neuron
        gid_equiv = (gid_0[1-t_pop] + 
                     int(1. * (gid - gid_0[t_pop]) / N[t_pop] * N[1-t_pop]))
        
        # Calculate minimal adjacent first spike time.
        # Exclude neurons, that are stimulated directly. 
        if (((gid < gid_0[t_pop] + N[t_pop]/2) or 
           (gid > gid_0[t_pop] + N[t_pop]/2 + params_stim['num_neuron'][t_pop]))
           and (gid_equiv in Sd_sorted_lim[1-t_pop]['ids'])):
        
            id_equiv = Sd_sorted_lim[1-t_pop]['ids'].index(gid_equiv)
            t_ref=Sd_sorted_lim[1-t_pop]['times'][id_equiv][0]     
        
            Lag[t_pop] += [Sd_sorted_lim[t_pop]['times'][n_id][0] - t_ref]
        
            trains_corrected += [Sd_sorted_lim[t_pop]['times'][n_id] - t_ref]
  
            t_min = min(t_min,trains_corrected[-1][0])
            t_max = max(t_max,trains_corrected[-1][-1])    
    
    PH_temp = f.f_psth(trains_corrected, bin_width=bin_width_PH,
                       t_min=t_min, t_max=t_max)
    
    PH += [{'times': PH_temp['times'], 
            'rates': np.mean(PH_temp['rates'], axis=0)}]
  

#######################################
#   PLOT  SPIKE DATA  #################
####################################### 

# initialization: assign default style format
style = aw.style_Plos()
aw.style_apply(style)

fig = plt.figure(1)
ax_id = 1
fig_size = [2,2]
axes = []


# plot raster plot
ax = plt.subplot2grid(fig_size,(0, 0),colspan=2)
ax_id += 2

# plot exc and inh populations mixed
for t_pop in [1,0]:
    
    # shift neuronal IDs in order to start from zero
    ids_temp = np.array(Sd_sorted[t_pop]['ids'])
    ids_temp -= min(ids_temp)
    for id_local in range(len(ids_temp)):
        # match ID to geometrical position 
        geom_loc = 1. * chain_len * ids_temp[id_local]//N[t_pop]
        
        ax.plot(Sd_sorted[t_pop]['times'][id_local]/1000., geom_loc / 1000. * 
                np.ones(len(Sd_sorted[t_pop]['times'][id_local])),
                'o'+['r','b'][t_pop], markersize=0.5, markeredgecolor=['r','b'][t_pop],alpha=0.5)

ax.set_xlabel('time (s)')
ax.set_ylabel('chain position (mm)')
axes += [{'ax':ax, 'cb':None}]


# plot exc-inh first spike lag
ax = fig.add_subplot(fig_size[0], fig_size[1], ax_id)
ax_id += 1

ax.hist(Lag[0], bins=100)
ax.set_xlabel('exc-inh lag (ms)')
ax.set_xlim([-500.,500.])
ax.set_yticks([])
ax.set_ylabel(' \n ')

axes += [{'ax':ax, 'cb':None}]


# plot time histogram
ax = fig.add_subplot(fig_size[0], fig_size[1], ax_id)
ax_id += 1

for t_pop in [0,1]:
    ax.plot(PH[t_pop]['times']/1000., PH[t_pop]['rates'], ['r','b'][t_pop])

# plot scaled excitatory trace
ax.plot(PH[0]['times']/1000., PH[0]['rates'] * np.max(PH[1]['rates']) / 
                              np.max(PH[0]['rates']), '0.5', linewidth=0.7)

ax.set_xlabel('time into Up state (s)')
ax.set_xticks([0.,1.5])
ax.set_xlim([-1.,1.5])
ax.set_ylabel('rate (spikes/s)')

axes += [{'ax':ax, 'cb':None}]


#Finalization phase: adjust subplot dimensions
aw.tight_layout(axes,style, './fig5.png', fig_width='medium',
               label_order=['A','B','C'])




plt.show()

