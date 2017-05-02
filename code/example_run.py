# -*- coding: utf-8 -*-
"""
Script that runs a single instance of the simple model

@author: Christoph Metzner
"""

from simple_model_class import simpleModel


# Parameters
n_ex    = 20  				# number of excitatory cells
n_inh   = 10 				# number of inhibitory cells

eta     = 5.0  				# 
tau_R   = 0.1 				# 'synaptic rise time'
tau_ex  = 2.0   			# exc. 'synaptic decay time'
tau_inh = 28.0  	 			# inh. 'synaptic decay time' (tau_inh=28 for the 
    

g_ee    = 0.015				# E-E weight
g_ei    = 0.025 			# default=0.025 E-I weight
g_ie    = 0.015 			# default=0.015 I-E weight
g_ii    = 0.02 				# default=0.02 I-I weight
g_de    = 0.3             		# Drive-E weight
g_di    = 0.08 				# default=0.08 Drive-I weight



#dt = 0.05				# time step
s=2**13
time 	= 500				# simulation stime (in ms) 
dt=float(time)/float(s)
print('time step:',dt)

b_ex    = -0.01 			# applied current for exc. cells
b_inh   = -0.01 			# applied current for inh. cells
drive_frequency = 40.0			# frequency of the periodic drive input
background_rate = 33.3			# average frequency of the Poissonian noise input
A = 0.5					# strength factor of the Poissonian noise
seed = 1234567				# seed for the random generator

filename = 'example'			# filename, in case data is recorded
directory = ''				# directory where data will be stored

model = simpleModel(n_ex,n_inh,eta,tau_R,tau_ex,tau_inh,g_ee,g_ei,g_ie,g_ii,g_de,g_di,dt,b_ex,b_inh,drive_frequency,background_rate,A,seed,filename,directory)

meg,ex,inh = model.run(time,0,0,0)

pxx,freqs = model.calculatePSD(meg,time)
model.plotPSD(freqs*1000,pxx,50.0,0)	# note that since time is in ms, frequencies of the PSD are off by a factor 1000 



model.rasterPlot(ex,time,0,'Ex')
model.rasterPlot(inh,time,0,'Inh')
model.plotMEG(meg,time,0)
