#----------------------imports and environment---------------------------------#
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
from ANNarchy import *
import numpy as np
#from net_fix import *
from net_homeostatic import *
from matplotlib import cm
from cmap import myCmap
"""
Python script for reproduce the rate code task from the Clopath et al. 2010
publication (Fig. 4 a). Network consists of ten, recurrent connected neurons.
Every neuron receives input from one neuron with a Poisson distributed
activity as stimulating input.
Poisson neurons firing with different frequencies (2Hz, 4Hz, 6Hz, ..)
as described in the original publication.
In the original publication, the resulting weights are averaged over 100s.
"""

###global parameter###
# duration of one presentation time
duration = 1000 #ms

#-----------------------population definitions----------------------------------#
"""
The 'PoissonPopulation' object is used to create a population with a Poisson distributed
activity with 10 neurons and a initial firing rate of 100 Hz.
The population of the ten AdEx neurons with the neuron model is created with the
'Population' object.
"""
poisPop = PoissonPopulation(geometry=10,rates=100.0)
pop_Ten = Population(geometry=10,neuron=AdExNeuron, name="pop_Ten")
#-----------------------projection definitions---------------------------------#
"""
Every neuron in the 'PoissonPopulation' is with one neuron in the neuron
population connected. Every spike of a Poisson neuron leads to a spike
of the corresponding AdEx neuron.
"""
projInp_Ten = Projection(
    pre = poisPop,
    post= pop_Ten,
    target='Exc'
).connect_one_to_one(weights = 30.0)

# create a Projection for the recurrent connections
projTen_Ten = Projection(
    pre= pop_Ten,
    post= pop_Ten,
    target= 'Exc',
    synapse= ffSyn
).connect_all_to_all(weights = 0.1,allow_self_connections=True)
# set network parameter
projTen_Ten.wMax= 0.25
projTen_Ten.uref= 60.0
#------------------------------main function-----------------------------------#
def run():
    # compile the network with ANNarchy
    compile()
    # number of experiment repetitions
    repeats = 20
    w = np.zeros((repeats,10,10)) #
    print('Start rate code experiment')
    for r in xrange(repeats):
        # repeat in 100 times to get the 100s as in Clopath et al. 2010
        for i in xrange(100):
            # set the different firint rates
            poisPop.rates = np.linspace(20,2,10)
            # simulate the network for 1000 ms
            simulate(duration)
            # save the changes in the weights
            w[r,:,:] += projTen_Ten.w
            # reset the network for the next 1000 ms
            reset()
        w[r,:,:] = w[r,:,:]/100


    img = np.ones((10,10))
    w = np.mean(w,axis=0)

    """
    Adapt the output like in Clopath et al. 2010.
    Depending on the connection, set another number to get another color.

    prepare a matrix of weights with different values for the different
    connections as mentioned in the Clopath et al., 2010 publication
    weak connections (< (2/3 of max. weight)) == 0
    strong unidirectional (> (2/3 of max. weight)) connections == 1.0
    strong bidirectional (> (2/3 of max. weight)) connections == 2.0
    """
    # Weak connection (under 2/3 of maximum weight) have the value = 0.0
    maxima = (np.nanmax(w)*2./3.)
    idx_b = np.where(w < maxima)
    img[idx_b] = 0.0

    # Strong biidirectional connections (> 2/3 of maximum weight) = 2.0
    idx_r = np.asarray(np.where(w >=maxima))
    for i in range(len(idx_r[0])):
        ix = (idx_r[0,i],idx_r[1,i])
        for j in range(len(idx_r[0])):
            ix2 = (idx_r[0,j],idx_r[1,j])
            if ix2 == (ix[1],ix[0]):
                img[ix[0],ix[1]] = 2.0
                img[ix[1],ix[0]] = 2.0

    # Strong unidirectional connections (> 2/3 of maximum weight) are every else

    # set the main diagonal to nan
    for i in xrange(10):
        w[i,i] = np.nan
        img[i,i] = np.nan

    ##--- start plotting ---##
    plt.figure()
    plt.imshow(img.T,interpolation='none',cmap= myCmap(),vmin=0,vmax=2)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel('Neuron Post',fontsize=20)
    plt.ylabel('Neuron Pre',fontsize=20)
    plt.savefig('rate_Code_Weights.png',bbox_inches='tight')

    print("finish")
#------------------------------------------------------------------------------#
if __name__ == "__main__":
    run()
