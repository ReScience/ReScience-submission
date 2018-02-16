#----------------------imports and environment---------------------------------#
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
from ANNarchy import *
import numpy as np
from net import *
from matplotlib import cm
"""
Python script for reproduce the rate code task from the Clopath et al. 2010
publication. Netowrk consists of ten, recurrent connected neurons.
Every neuron receives input from one neuron with a poisson distributed
activity as stimulating input.
Poisson neurons firing with different frequencies (2Hz, 4Hz, 6Hz, ..)
as described in the original publication.
In the original publication, the resulting weights are averaged over 100s.
"""

###global parameter###
# duration of one presentation time
duration = 1000 #ms

#-----------------------population defintions----------------------------------#
poisPop = PoissonPopulation(geometry=10,rates=100.0)
pop_Ten = Population(geometry=10,neuron=spkNeurV1, name="pop_Ten")
#-----------------------projection definitions---------------------------------#
projInp_Ten = Projection(
    pre = poisPop,
    post= pop_Ten,
    target='Exc'
).connect_one_to_one(weights = 30.0)

projTen_Ten = Projection(
    pre= pop_Ten,
    post= pop_Ten,
    target= 'Exc',
    synapse= ffSyn               #Uniform(0.0,0.01)
).connect_all_to_all(weights = 0.1,allow_self_connections=True)
# set network parameter
projTen_Ten.wMax= 0.3
projTen_Ten.vmean= 80.0

#------------------------------main function-----------------------------------#
def run():
    # compile the network with ANNarchy
    compile()
    # number of experiment repetitions
    repeats = 1
    w = np.zeros((repeats,10,10)) #
    print('Start rate code experiment')
    for wi in xrange(repeats):
        # repeat in 100 times to get the 100s as in Clopath et al. 2010
        for i in xrange(100):
            # set the different firint rates
            poisPop.rates = np.linspace(2,20,10)
            # simulate the network for 1000 ms
            simulate(duration)
            # save the changes in the weights
            w[wi,:,:] += projTen_Ten.w
            # reset the network for the next 1000 ms
            reset(populations=True, projections=True, synapses=True)
        w[wi,:,:] = w[wi,:,:]/100


    img = np.ones((10,10))
    w = np.mean(w,axis=0)

    """
    Adapt the output like in Clopath et al. 2010.
    Depending on the connection, set another number to get another color.
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

    # plot the image like in Clopath 2010
    # note! here we use an other color map:
    # weak weights are yellow
    # unidirectional are ligth green
    # bidirectional are dark green
    plt.figure()
    plt.imshow(img.T,interpolation='none',cmap= mp.cm.get_cmap('summer_r'))
    plt.colorbar()
    plt.xlabel('Neuron Pre')
    plt.ylabel('Neuron Post')
    plt.savefig('rate_Code_Weights.png')

    # plot the original weight matrix
    plt.figure()
    plt.imshow(w,interpolation='none',cmap=plt.get_cmap('summer_r'))
    plt.colorbar()
    plt.xlabel('Neuron Post')
    plt.ylabel('Neuron Pre')
    plt.savefig('rate_Code_Matrix.png')
    print("finish")
#------------------------------------------------------------------------------#
if __name__ == "__main__":
    run()
