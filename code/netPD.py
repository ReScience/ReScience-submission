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

def runParams(tsim=1.0, bg_type=0, stim=0, w_ex=87.8, g=4.0, bg_freq=8.0, nsyn_type=0, thal='OFF', filename=None):

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
    neurons.v = '-58.0*mV + 10.0*mV*randn()'
    neurons.I = 0.0*pA      # initial value for synaptic currents
    neurons.Iext = 0.0*pA   # constant external current

    pop, con, bg_in, smon_net, thal_input, thal_con = netMod.PDnet(neurons, stim, \
                                            bg_type, w_ex, g, bg_freq, nsyn_type, thal)

    ###############################################################################
    # Running the simulation
    ###############################################################################
    net = Network(collect())

    if (thal == 'OFF'):
        net.add(neurons,pop, con, bg_in)    # Adding objects to the simulation
        net.run(tsim*second, report='stdout')

    elif (thal == 'ON'):
        w_thal = w_ex*pA            # excitatory synaptic weight from thalamus
        std_w_thal = w_thal*0.1     # standard deviation weigth
        net.add(neurons,pop, con, bg_in, thal_input, thal_con)    # Adding objects to the simulation

        for repeat in range(0,int(tsim)):
            net.run(0.7*second,report='stdout')
            gc.collect()    #garbage collector to clean memory

            # Adding thalamic input
            for r in range(0,8):
                thal_con[r].w = 'clip((w_thal + std_w_thal*randn()),w_thal*0.0, w_thal*inf)'
            net.run(0.01*second, report='stdout')
            gc.collect()    #garbage collector to clean memory

            # Removing thalamic input
            for r in range(0,8):
                thal_con[r].w = 0
            net.run(0.29*second, report='stdout')
            gc.collect()    #garbage collector to clean memory

    ###############################################################################
    # Saving raster plot
    ###############################################################################
    savetxt(filename, c_[smon_net.i,smon_net.t/ms],fmt="%i %.2f")
