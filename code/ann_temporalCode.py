#----------------------imports and environment---------------------------------
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
from ANNarchy import *
import numpy as np
from net import *

###global parameter###
duration = 200 #ms
#----------------------defint time points of spikes-----------------------------#
spike_times =[[1],[2],[3],[4],[5],[6],[7],[8],[9],[10]]
#-----------------------population defintions-----------------------------------#
inpPop = SpikeSourceArray(spike_times=spike_times)
pop_Ten = Population(geometry=10,neuron=spkNeurV1, name="pop_Ten")
#-----------------------projection definitions----------------------------------
projInp_Ten = Projection(
    pre = inpPop,
    post=pop_Ten,
    target='Exc'
).connect_one_to_one(weights = 30.0)

projTen_Ten = Projection(
    pre=pop_Ten, 
    post=pop_Ten, 
    target='Exc',
    synapse=ffSyn
).connect_all_to_all(weights = 0.5,allow_self_connections=True)
projTen_Ten.wMax=0.3

#------------------------------main function------------------------------------
def run():

    compile()
    for i in xrange(500):
        spkT_N1 = [0+(i*duration)]
        spkT_N2 = [20+(i*duration)]
        spkT_N3 = [40+(i*duration)]
        spkT_N4 = [60+(i*duration)]
        spkT_N5 = [80+(i*duration)]
        spkT_N6 = [100+(i*duration)]
        spkT_N7 = [120+(i*duration)]
        spkT_N8 = [140+(i*duration)]
        spkT_N9 = [160+(i*duration)]
        spkT_N10 = [180+(i*duration)]

        inpPop.spike_times=[spkT_N1,spkT_N2,spkT_N3,spkT_N4,spkT_N5,spkT_N6,spkT_N7,spkT_N8,spkT_N9,spkT_N10]
        
        simulate(duration)
    w = projTen_Ten.w
    img = np.ones((10,10))
    maxima = (np.nanmax(w)*2./3.)
    idx = np.where(w < maxima)
    img[idx[0],idx[1]] = 0.0

    for i in xrange(10):
        w[i][i] = np.nan
        img[i,i]= np.nan
    print(np.shape(w))
    plt.figure()
    plt.imshow(w,interpolation='none',cmap=plt.get_cmap('summer_r',3))
    plt.colorbar()
    plt.savefig('temporal_Code_w.png')


    plt.figure()
    plt.imshow(img,interpolation='none',cmap=plt.get_cmap('summer_r',3))
    plt.colorbar()
    plt.savefig('temporal_Code.png')
    print("finish")
#------------------------------------------------------------------------------------
run()
