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

from brian2 import *
from netParams import *
import neuronModels as neuronMod
import netModels as netMod

def runParams(tsim=1.0, bg_type=0, stim=0, w_ex=87.8, g=4.0, bg_freq=8.0, filename=None):

    ###############################################################################
    # Simulation parameters
    ###############################################################################
    defaultclock.dt = 0.1*ms    # timestep of numerical integration method

    # neuron model
    eqs = neuronMod.LIFmodel
    reset = neuronMod.resetLIF
    tau_m, tau_ref, Cm, v_r, v_th = neuronMod.LIFparams()

    ###############################################################################
    # Creating neurons
    ###############################################################################
    neurons = NeuronGroup(N, eqs, threshold='v>v_th', reset=reset, \
                            method='linear', refractory=tau_ref)

    # seting initial values for membrane potential and currents
    neurons.v = 'v_r + 0.1*v_r*randn()'
    neurons.I = 0.0*pA      # initial value for synaptic currents
    neurons.Iext = 0.0*pA   # constant external current

    pop, con, bg_in, smon_net = netMod.PDnet(neurons, stim, bg_type, w_ex, g, bg_freq)

    ###############################################################################
    # Running the simulation
    ###############################################################################
    net = Network(collect())
    net.add(neurons,pop, con, bg_in)    # Adding objects to the simulation
    net.run(tsim*second, report='stdout')

    ###############################################################################
    # Saving raster plot
    ###############################################################################
    savetxt(filename, c_[smon_net.i,smon_net.t/ms],fmt="%i %.2f")
