#----------------------imports and environment---------------------------------
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
from ANNarchy import *
import numpy as np
from net import *

###global parameter###
duration = 1000 #ms

#-----------------------population defintions-----------------------------------#
poisPop = PoissonPopulation(geometry=10,rates=100.0)
pop_Ten = Population(geometry=10,neuron=spkNeurV1, name="pop_Ten")
#-----------------------projection definitions----------------------------------
projInp_Ten = Projection(
    pre = poisPop,
    post=pop_Ten,
    target='Exc'
).connect_one_to_one(weights = 30.0)

projTen_Ten = Projection(
    pre=pop_Ten, 
    post=pop_Ten, 
    target='Exc',
    synapse=ffSyn               #Uniform(0.0,0.01)
).connect_all_to_all(weights = 0.1,allow_self_connections=True)
projTen_Ten.wMax=0.3
projTen_Ten.vmean=80.0


#------------------------------main function------------------------------------
def run():
    
    compile()
    repeats = 1
    w = np.zeros((repeats,10,10))
    for wi in xrange(repeats):
        for i in xrange(100):
            poisPop.rates = np.linspace(2,20,10)
            simulate(duration)
            w[wi,:,:] += projTen_Ten.w
            reset(populations=True, projections=True, synapses=True)
        w[wi,:,:] = w[wi,:,:]/100
        

    img = np.ones((10,10))


    w = np.mean(w,axis=0)

    maxima = (np.nanmax(w)*2./3.)
    idx_b = np.where(w < maxima)
    img[idx_b] = 0.0

    idx_r = np.asarray(np.where(w >=maxima))
    for i in range(len(idx_r[0])):
        ix = (idx_r[0,i],idx_r[1,i])
        for j in range(len(idx_r[0])):
            ix2 = (idx_r[0,j],idx_r[1,j])
            if ix2 == (ix[1],ix[0]):
                img[ix[0],ix[1]] = 2.0
                img[ix[1],ix[0]] = 2.0

    for i in xrange(10):
        w[i,i] = np.nan
        img[i,i] = np.nan

    plt.figure()
    plt.imshow(w.T,interpolation='none',cmap=plt.get_cmap('summer_r',3))
    plt.colorbar()
    plt.xlabel('Neuron Pre')
    plt.ylabel('Neuron Post')
    plt.savefig('rate_Code_Poisson.png')

    plt.figure()
    plt.imshow(img.T,interpolation='none',cmap=plt.get_cmap('summer_r',3))
    plt.colorbar()
    plt.xlabel('Neuron Post')
    plt.ylabel('Neuron Pre')
    plt.savefig('rate_Code_Poisson_Transpose.png')
    print("finish")
#------------------------------------------------------------------------------------
run()
