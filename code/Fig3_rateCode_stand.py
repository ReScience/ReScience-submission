"""
Python script for reproduce rate code task from the Clopath et al. (2010)
with the STDP rule from Song and Abbott (2001) as mentioned in Clopath et al. (2010)(Fig. 4 c).
Network consists of ten, recurrent connected neurons.
Every neuron receives input from one neuron with a Poisson distributed
activity as stimulating input.
Poisson neurons firing with different frequencies (2Hz, 4Hz, 6Hz, ..)
as described in the original publication.
In the original publication, the resulting weights are averaged over 100s.
We use the in ANNarchy implemented version of the Song and Abbott (2001) learning rule.
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from ANNarchy import *
setup(dt=1.0,seed=314)
from network import *
from cmap import myCmap

# Define the pair based STDP model after Song, S., and Abbott, L.F. (2001)
# Taken from the ANNarchy source code
STDP_Song = Synapse(
    parameters = """
        tau_plus = 7.0  : projection
        tau_minus = 10.0 : projection
        A_plus = 0.1    : projection
        A_minus = 0.1   : projection
        w_min = 0.0      : projection
        w_max = 0.75      : projection
    """,
    equations = """
        tau_plus  * dx/dt = -x : event-driven
        tau_minus * dy/dt = -y : event-driven
    """,
    pre_spike="""
        g_target += w
        x += A_plus * w_max
        w = clip(w + y, w_min , w_max)
    """,
    post_spike="""
        y -= A_minus * w_max
        w = clip(w + x, w_min , w_max)
    """
)

def run_Rate():

    # Global parameters
    duration = 1000 #ms

    # Populations
    """
    The 'PoissonPopulation' object is used to create a population with a Poisson distributed
    activity with 10 neurons and a initial firing rate of 100 Hz.
    The population of the ten AdEx neurons with the neuron model is created with the
    'Population' object.
    """
    poisPop = PoissonPopulation(geometry=10,rates=100.0)
    pop_Ten = Population(geometry=10,neuron=AdExNeuron, name="pop_Ten")

    # Projections
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

    # Create a Projection for the recurrent connections
    # use the Song and Abbott (2001) learning rule
    projTen_Ten = Projection(
        pre= pop_Ten,
        post= pop_Ten,
        target= 'Exc',
        synapse= STDP_Song
    ).connect_all_to_all(weights =0.1,allow_self_connections=True)


    # Compile the network with ANNarchy
    compile()
    # Number of repetitions
    repeats = 20
    w = np.zeros((repeats,10,10)) #
    print('Start rate code experiment')
    for r in range(repeats):
        # Repeat in 100 times to get the 100s as in Clopath et al. 2010
        for i in range(100):
            # Set the different firint rates
            poisPop.rates = np.linspace(20,2,10)
            # Simulate the network for 1000 ms
            simulate(duration)
            # Save the changes in the weights
            w[r,:,:] += projTen_Ten.w
            # Reset the network for the next 1000 ms
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

    # Strong bidirectional connections (> 2/3 of maximum weight) = 2.0
    idx_r = np.asarray(np.where(w >=maxima))
    for i in range(len(idx_r[0])):
        ix = (idx_r[0,i],idx_r[1,i])
        for j in range(len(idx_r[0])):
            ix2 = (idx_r[0,j],idx_r[1,j])
            if ix2 == (ix[1],ix[0]):
                img[ix[0],ix[1]] = 2.0
                img[ix[1],ix[0]] = 2.0

    # Strong unidirectional connections (> 2/3 of maximum weight) are every else

    # Set the main diagonal to nan
    for i in range(10):
        w[i,i] = np.nan
        img[i,i] = np.nan

    # Start plotting
    plt.figure()
    plt.imshow(img.T,interpolation='none',cmap= myCmap(),vmin=0,vmax=2)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel('Neuron Post',fontsize=20)
    plt.ylabel('Neuron Pre',fontsize=20)
    plt.savefig('Fig3_rateCode_standardSTDP.png',bbox_inches='tight')
    plt.show()

    print("done")


if __name__ == "__main__":
    run_Rate()
