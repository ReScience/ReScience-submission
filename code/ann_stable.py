import matplotlib as mp
mp.use('TKAgg')
import matplotlib.pyplot as plt
import numpy as np
from ANNarchy import *
setup(dt=1.0)
# import of the model and learning rule, defined in net.py
from net_homeostatic import *

"""
Python script to reproduce the stable learning task
of the Clopath et al. 2010 publication.
See Fig.5 a in original publication.
Algorithm based on Matlab code of the Clopath et al. 2010 model.
Available on modelDB:
https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=144566
"""
#--------------- define the presynaptic neuron model --------------------------#
"""
Because of the learning rule, we need a additional layer, that contains the
necessary variables for the learning. This population is one to one connected
with the Poisson input layer and spiked for ever corresponding neuron in the
Poisson layer.
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
#-----------------------global variables---------------------------------------#
nb_pre = 500 # number of input neurons
nb_post= 1 # number of post synaptic neuron
duration = 100 #ms # number of time steps per epoch in ms
nb_epochs = 1000 # number of epochs per input pattern
#-----------------------population defintions----------------------------------#
pre_pop = Population(geometry=nb_pre,neuron=neuron)
post_pop= Population(geometry=nb_post,neuron=spkNeurV1)
#-----------------------projection definitions---------------------------------#
# Projection object to initialise the synapse with the learning rule
projInp_N = Projection(
    pre = pre_pop,
    post= post_pop,
    target='Exc',
    synapse = ffSyn
).connect_all_to_all(weights = Uniform(0.0,2.0))
#----------------------------define input--------------------------------------#
# input parameters as in Matlab source code
sigma = 10
in_max = 0.015
in_min = 0.0001
nb_pattern = 10
patterns = np.zeros((nb_epochs,duration))
for i in range(nb_epochs):
    patterns[i,:] = np.floor(np.random.rand()*nb_pattern)
patterns = np.reshape(patterns,nb_epochs*duration)

# initialise the gaussian input
ind=np.linspace(0,nb_pre-1,nb_pre)
gau= in_min + in_max*np.exp( - ( ind - nb_pre/2.)**2 / (2*sigma**2))
gau = np.append(gau,gau)
input_patterns = np.zeros((nb_pattern,nb_pre))
for i in range(nb_pattern):
    mup = 1+(i)*nb_pre/nb_pattern;
    input_patterns[i,:] = gau[mup:mup+nb_pre]


compile()# Compile the network
# set parameters analoug to the parameters in the Matlab source code
projInp_N.transmit = 4.0
 # learning speed ten times higher to speed up the learning
projInp_N.aLTP = 10*0.00008
projInp_N.aLTD = 10*0.00014
projInp_N.wMax =3.0

# monitor object to save the weight, after each epoch to save memory
monW = Monitor(projInp_N,'w',period=duration)

# start the simulation
for t in range(1,duration*nb_epochs):
    inp = ((np.random.rand(nb_pre))< input_patterns[int(patterns[t])])*1
    pre_pop.g_vm = -60+(inp*30)
    simulate(1)

# get the weights from monitor
w = monW.get('w')
# plot the weights during time
plt.figure()
plt.imshow(np.squeeze(w).T)
plt.xlabel('Number of epoch')
plt.ylabel('Synapse index')
plt.savefig('weights_stable.png')

print('Finish with simulation!')
