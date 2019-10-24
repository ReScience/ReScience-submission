"""
Python script to reproduce the stable learning task
of the Clopath et al. 2010 publication.
See Fig.5 a in original publication.
Algorithm based on Matlab code of the Clopath et al. 2010 model.
Available on modelDB:
https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=144566
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from ANNarchy import *
setup(dt=1.0,seed=23456)

from network import *

# Presynaptic neuron model
"""
Because of the learning rule, we need an additional layer, that contains the
necessary variables for the learning. This population is one to one connected
with the Poisson input layer and spiked for ever corresponding neuron in the
Poisson layer. For a further description to define a neuron model, look at the
'net_fix.py' or 'net_homeostatic.py' file.
"""

params = """
EL = -70.4      :population
VTrest = -50.4  :population
taux = 15.0     :population  """

eqs = """
dg_vm/dt = EL/1000 : min = EL, init=-70.4
Spike = if state == 1: 1.0 else: 0.0
dReset/dt = if state == 1: +1 else: -Reset
dxtrace/dt = if state == 1: +1/taux else: -xtrace/taux : init = 0.0
state = if state >0: -1 else: 0"""

neuron = Neuron(parameters = params,
                equations = eqs,
                reset = """ g_vm = EL
                            state = 1""",
                spike = """g_vm > VTrest""")

# Global parameters
nb_pre = 500 # number of input neurons
nb_post= 1 # number of post synaptic neuron
duration = 100 #ms # number of time steps per epoch in ms
nb_epochs = 1000 # number of epochs per input pattern

# Population defintions
"""
Create the populations of presynaptic neurons and the population
of postsynaptic AdEx neurons.
"""
pre_pop = Population(geometry=nb_pre,neuron=neuron)
post_pop= Population(geometry=nb_post,neuron=AdExNeuron)
# Projection definitions

# Projection object to initialise the synapse with the learning rule
projInp_N = Projection(
    pre = pre_pop,
    post= post_pop,
    target='Exc',
    synapse = ffSyn
).connect_all_to_all(weights = Uniform(0.0,2.0))

projInp_N.set_fix = 0.0 # use the homeostatic mechanisms in the LTD term

# Define the input parameters as in Matlab source code
sigma = 10
in_max = 0.015
in_min = 0.0001
nb_pattern = 10

def run():
    # The input generating is taken from the original Matlab source code

    patterns = np.zeros((nb_epochs,duration))
    for i in range(nb_epochs):
        patterns[i,:] = np.floor(np.random.rand()*nb_pattern)
    patterns = np.reshape(patterns,nb_epochs*duration)

    # Initialise the gaussian input
    ind=np.linspace(0,nb_pre-1,nb_pre)
    gau= in_min + in_max*np.exp( - ( ind - nb_pre/2.)**2 / (2*sigma**2))
    gau = np.append(gau,gau)
    input_patterns = np.zeros((nb_pattern,nb_pre))
    for i in range(nb_pattern):
        mup = 1+(i)*nb_pre/nb_pattern;
        input_patterns[i,:] = gau[int(mup):int(mup+nb_pre)]


    compile()# Compile the network
    # Set parameters analoug to the parameters in the Matlab source code
    projInp_N.transmit = 4.0
    projInp_N.aLTP = 10*0.00008
    projInp_N.aLTD = 10*0.00014
    projInp_N.wMax = 3.0

    # Monitor object to save the weight, after each epoch to save memory
    monW = Monitor(projInp_N,'w',period=duration)

    # Start the simulation
    for t in range(1,duration*nb_epochs):
        inp = ((np.random.rand(nb_pre))< input_patterns[int(patterns[t])])*1

        # Set the membrane potential (vm) of the presynaptic neuron to emit a
        # spike depending on the input pattern
        pre_pop.g_vm = -60+(inp*30)
        simulate(1)

    # Get the weights from monitor
    w = monW.get('w')

    # Start plotting
    plt.figure()
    plt.imshow(np.squeeze(w).T, cmap='gray')
    plt.xlabel('Number of epoch')
    plt.ylabel('Synapse index')
    plt.savefig('Fig4_stable.png',bbox_inches='tight')
    plt.show()
    print("done")

if __name__ == "__main__":
    run()
