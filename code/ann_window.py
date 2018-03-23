#----------------------imports and environment---------------------------------
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
from ANNarchy import *
import numpy as np
from net_fix import *

"""
Python script to reproduce the STDP window protocol.
Record the change in the synaptic weight for different time intervalls
between pre- and postsynaptic spike.
"""

###global parameter###
duration = 1700 #ms
d_t = 15
#----------------------defint time points of spikes-----------------------------#
# create a list of the different spike times for the pre synaptic neuron
spike_times1 =np.asarray(range(50,duration+50,50))

l = np.asarray(range(-d_t,d_t+1,1))
# create the list for the different spike times for the post synaptic neuron
spike_times2 =spike_times1[0:31]+l #[[40,50,60]]
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
).connect_one_to_one(weights = 30.0)

projST2_V1 = Projection(
    pre=inpPop2,
    post=popN2,
    target='Exc'
).connect_one_to_one(weights = 30.0)

projV1_V1 = Projection(
    pre=popN1,
    post=popN2,
    target='Exc',
    synapse=ffSyn
).connect_one_to_one(weights = 0.1)

#---- parameter adjustments ----#
projV1_V1.thetaLTD = -60.6
projV1_V1.vmean = 80.0
projV1_V1.transmit = 1.0 # to activate the transmission over the synapse
#------------------------------main function------------------------------------
def run():

    compile()

    #------- neuron Monitors --------#
    # monitor to save changes in the synapse
    dendrite = projV1_V1.dendrite(0)
    m_d = Monitor(dendrite, ['w','deltaW','ltdTerm','ltpTerm'])

    simulate(duration)

    #---get the recorded data-----#
    d_w = m_d.get('w')
    delta_w = m_d.get('deltaW')
    ltd_w = m_d.get('ltdTerm')
    ltp_w = m_d.get('ltpTerm')

    # sort the data dependent of the spiking order
    idx_LTD = np.asarray(np.where(delta_w[0:duration/2] < 0))
    ltd = delta_w[idx_LTD]
    idx_LTP = np.asarray(np.where(delta_w[duration/2:duration] > 0))
    ltp = delta_w[idx_LTP+duration/2]

    fig,ax = plt.subplots()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    plt.plot(np.linspace(0,14,len(ltd[0])),ltd[0],color='tomato',lw=6.0)
    plt.plot(np.linspace(16,31,len(ltp[0])),ltp[0],color='steelblue',lw=6.0)
    plt.plot((14,16), (ltd[0][-1],ltp[0][0]),color='k')
    plt.axhline(y=0.0, color='k',linestyle='--')
    plt.axvline(x=15,color='k',linestyle='--')
    plt.ylim(ymin=-0.005,ymax=0.005)
    plt.yticks(np.linspace(-0.004,0.004,5),np.linspace(60,140,5),fontsize=20)
    plt.ylabel('Normalized weight (%)',fontsize=25)
    plt.xlim(xmin=-3,xmax=33)
    plt.xticks(np.linspace(5,25,3),np.linspace(-10,10,3),fontsize=20)
    plt.xlabel('T (ms)',fontsize=25)
    plt.savefig('deltaW.png',bbox_inches='tight', pad_inches = 0.1)

    print("finish")
#------------------------------------------------------------------------------------
run()
