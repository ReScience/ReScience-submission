# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2017, Christoph Metzner
# Distributed under the (new) BSD License.
#
# Contributors: Christoph Metzner (c.metzner@herts.ac.uk)
# -----------------------------------------------------------------------------
# References:
#
# * Vierling-Claassen, D., Siekmeier, P., Stufflebeam, S., & Kopell, N. (2008). 
#   Modeling GABA alterations in schizophrenia: a link between impaired 
#   inhibition and altered gamma and beta range auditory entrainment. 
#   Journal of neurophysiology, 99(5), 2656-2671.
# -----------------------------------------------------------------------------
# A script that runs the main simulations of the replication study.
#
# -----------------------------------------------------------------------------

import time
import numpy as np

from simple_model_class import simpleModel


# Parameters
n_ex    = 20                  # number of excitatory cells
n_inh   = 10                 # number of inhibitory cells

eta     = 5.0                  # synaptic scaling parameter
tau_R   = 0.1                 # 'synaptic rise time'
tau_ex  = 2.0               # exc. 'synaptic decay time'
tau_inh_ctrl = 8.0               # inh. 'synaptic decay time' (control case)
tau_inh_schiz = 28.0               # inh. 'synaptic decay time' (schizophrenia case)
 
    

g_ee    = 0.015                # E-E weight
g_ei    = 0.025             # E-I weight
g_ie    = 0.015                # I-E weight
g_ii    = 0.02                 # I-I weight
g_de    = 0.3                     # Drive-E weight
g_di    = 0.08                 # Drive-I weight



                
s = 2**13
sim_time = 500                # simulation stime (in ms) 
dt = float(sim_time)/float(s)        # time step


b_ex    = -0.01     # applied current for exc. cells
b_inh   = -0.01     # applied current for inh. cells
frequencies = [40.0,30.0,20.0]
background_rate = 33.3
A = 0.5


seeds = np.load('../data/Seeds.npy')
directory = '../data'

start = time.time() 
for drive_frequency in frequencies:
    print('Drive Frequency:',drive_frequency)
    meg_ctrl_avg = np.zeros((8192,))
    meg_schiz_avg = np.zeros((8192,))
    for i,seed in enumerate(seeds):
        print(i)

        filename_ctrl  = '/Single-Trial-Data-'+str(int(drive_frequency))+\
            'Hz/sims_ctrl_'+str(int(drive_frequency))+'Hz'+str(i)
        filename_schiz = '/Single-Trial-Data-'+str(int(drive_frequency))+\
            'Hz/sims_schiz_'+str(int(drive_frequency))+'Hz'+str(i)
            
        model_ctrl = simpleModel(n_ex,n_inh,eta,tau_R,tau_ex,tau_inh_ctrl,
            g_ee,g_ei,g_ie,g_ii,g_de,g_di,dt,b_ex,b_inh,drive_frequency,background_rate,
            A,seed,filename_ctrl,directory)
        model_schiz = simpleModel(n_ex,n_inh,eta,tau_R,tau_ex,tau_inh_schiz,
            g_ee,g_ei,g_ie,g_ii,g_de,g_di,dt,b_ex,b_inh,drive_frequency,background_rate,
            A,seed,filename_schiz,directory)
            
        meg_ctrl,ex_ctrl,inh_ctrl = model_ctrl.run(sim_time,1,1,1)
        meg_schiz,ex_schiz,inh_schiz = model_schiz.run(sim_time,1,1,1)
        
        meg_ctrl_avg += meg_ctrl    
        meg_schiz_avg += meg_schiz
        
        


    meg_ctrl_avg = meg_ctrl_avg *0.05
    meg_schiz_avg = meg_schiz_avg *0.05

    filename_ctrl_avg = directory+'/Data-Average/sims_ctrl_avg_'+\
        str(int(drive_frequency))+'Hz'
    filename_schiz_avg = directory+'/Data-Average/sims_schiz_avg_'+\
        str(int(drive_frequency))+'Hz'

    np.save(filename_ctrl_avg,meg_ctrl_avg)
    np.save(filename_schiz_avg,meg_schiz_avg)


end = time.time()
print('Time elapsed:', end-start)
