#----------------------imports and environment---------------------------------
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
from ANNarchy import *
setup(dt=1)
import numpy as np
from net_fix import *

"""
Python script to reproduce the STDP window protocol.
Record the change in the synaptic weight for different time intervalls
between pre- and postsynaptic spike.
"""

###global parameter###
duration = 40 # duration time of 35 ms
dt = np.linspace(-15,15,31) # 31 spiking pairs for a time difference dt between a pre and post synaptic spike from -15 to 15 ms
print(dt)
initW = 0.012
#----------------------defint time points of spikes-----------------------------#
# create a list of the different spike times for the pre synaptic neuron
spike_times1 =np.asarray(range(50,duration+50,50))

# create the list for the different spike times for the post synaptic neuron
spike_times2 =np.asarray(range(50,duration+50,50))
#-----------------------population defintions-----------------------------------#
# two SpikeSourceArray populations to determine the spike timings of the
# pre and post neuron
inpPop1 = SpikeSourceArray(spike_times=spike_times1.tolist())
inpPop2 = SpikeSourceArray(spike_times=spike_times2.tolist())
popN1 = Population(geometry=1,neuron=spkNeurV1, name="N1")
popN2 = Population(geometry=1,neuron=spkNeurV1, name="N2")
#-----------------------projection definitions----------------------------------
projST1_V1 = Projection(
    pre=inpPop1,
    post=popN1,
    target='Exc'
).connect_one_to_one(weights = 25.0)

projST2_V1 = Projection(
    pre=inpPop2,
    post=popN2,
    target='Exc'
).connect_one_to_one(weights = 25.0)

projV1_V1 = Projection(
    pre=popN1,
    post=popN2,
    target='Exc',
    synapse=ffSyn
).connect_one_to_one(weights = initW)

#---- parameter adjustments ----#
#projV1_V1.thetaLTD = -65.6
projV1_V1.vmean = 70.0
#projV1_V1.transmit = 1.0 # to activate the transmission over the synapse
#------------------------------main function------------------------------------
def run():

    compile()

    #------- neuron Monitors --------#
    # monitor to save changes in the synapse
    dendrite = projV1_V1.dendrite(0)
    m_d = Monitor(dendrite, ['w','deltaW','ltdTerm','ltpTerm'])
    n_pairs = 31

    w = np.zeros(n_pairs)
    dW = np.zeros(n_pairs)
    for i in range(n_pairs):
        reset()
        projV1_V1.w = initW       
        inpPop1.spike_times = [16] # presynaptic neuron always spikes at t=20 ms
        inpPop2.spike_times = [16+dt[i]] # add the time difference to estimate the postsynaptic spike time
        simulate(duration)
        d_w = m_d.get('w')
        delta_w = m_d.get('deltaW')
        w[i] = d_w[-1]#np.mean(d_w)
        dW[i] = np.sum(delta_w)

    #---get the recorded data-----#

    ltd_w = m_d.get('ltdTerm')
    ltp_w = m_d.get('ltpTerm')

    w = w/initW *100.

    fig,ax = plt.subplots(figsize=(13,9))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    plt.plot(np.linspace(0,13,14),w[0:14],color='tomato',lw=4.0)
    plt.plot(np.linspace(13,17,5),np.linspace(w[13],w[17],5),color='black')
    plt.plot(np.linspace(17,30,14),w[17:31],color='steelblue',lw=4.0)

    plt.axhline(y=100.0, color='k',linestyle='--')
    plt.axvline(x=15,color='k',linestyle='--')
    plt.ylim(ymin=55,ymax=145)
#    plt.yticks(np.linspace(50,140,5),np.linspace(50,140,5),fontsize=20)
    plt.ylabel('Normalized weight (%)',fontsize=25)
    #plt.xlim(xmin=-3,xmax=33)
    plt.xticks(np.linspace(5,25,3),np.linspace(-10,10,3),fontsize=20)
    plt.xlabel('T (ms)',fontsize=25)

    plt.savefig('windW.png',bbox_inches='tight', pad_inches = 0.1)
   


    print("finish")
#------------------------------------------------------------------------------------
run()
