# ----------------------------------------------------------------------------
# Contributors: Renan O. Shimoura
#               Nilton L. Kamiji
#               Rodrigo F. O. Pena
#               Vinicius L. Cordeiro
#               Cesar C. Ceballos
#               Cecilia Romaro
#               Antonio C. Roque
# ----------------------------------------------------------------------------
# References:
#
# *The Cell-Type Specific Cortical Microcircuit: Relating Structure and Activity
# in a Full-Scale Spiking Network Model*,
# Tobias C. Potjans and Markus Diesmann,
# Cerebral Cortex, 24(3):785-806, 2014.
# ----------------------------------------------------------------------------
# File description:
#
# Main parameters of simulation is defined in this code.
# ----------------------------------------------------------------------------

import matplotlib
matplotlib.use('Agg')

from brian2 import *
from netParams import *
import neuronModels as neuronMod
import netModels as netMod

# Seting the simulation to run with openmp
#set_device('cpp_standalone', directory='PD')
#prefs.devices.cpp_standalone.openmp_threads = 10

###############################################################################
# Simulation parameters
###############################################################################
defaultclock.dt = 0.1*ms    # timestep of numerical integration method
tsim = 1.0*second           # time of simulation

# background type
bg_type = 0                 # 0 = layer-specific
                            # 1 = layer-independent

# stimulation
stim = 0                    # 0 = turn on the background noise
                            # 1 = DC current experiment
                            # 2 = background noise + DC current

# neuron model
eqs = neuronMod.LIFmodel
reset = neuronMod.resetLIF
tau_m, tau_ref, Cm, v_r, v_th = neuronMod.LIFparams()

# synapse parameters
w_ex = 87.8*pA			   # excitatory weight
g = 4.0             	   # inhibitory weight balance

# fixed seed of pseudo random numbers to test reproducibility
seed(1)

###############################################################################
# Creating neurons
###############################################################################
neurons = NeuronGroup(N, eqs, threshold='v>v_th', reset=reset, \
                        method='linear', refractory=tau_ref)

# seting initial values for membrane potential and currents
neurons.v = 'v_r + 0.1*v_r*randn()'
neurons.I = 0.0*pA      # initial value for synaptic currents
neurons.Iext = 0.0*pA   # constant external current

pop, con, bg_in, smon_net = netMod.PDnet(neurons, stim, bg_type, w_ex, g)

###############################################################################
# Running the simulation
###############################################################################
net = Network(collect())
net.add(neurons,pop, con, bg_in)    # Adding objects to the simulation
net.run(tsim, report='stdout')

###############################################################################
# Saving raster plot
###############################################################################
savetxt('data_raster_g' + str(g) + '_w' + str(w_ex/pA) + '.dat',\
            c_[smon_net.i,smon_net.t/ms],fmt="%i %.2f")
