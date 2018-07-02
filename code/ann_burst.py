#----------------------imports and environment---------------------------------
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from ANNarchy import *
import numpy as np
from net_fix import *


initW1 = 0.007
initW2 = 0.0055
initW3 = 0.008

"""
Python script to reproduce the STDP window protocol.
Record the change in the synaptic weight for different time intervalls
between pre- and postsynaptic spike.
"""

###global parameter###
duration = 240 #ms
d_t = 10
t_1 = 110
#----------------------defint time points of spikes-----------------------------#
# create a list of the different spike times for the pre synaptic neuron
spike_times1 = [t_1]

# create the list for the different spike times for the post synaptic neuron
spike_times2 = [t_1]
#-----------------------population defintions-----------------------------------#
# two SpikeSourceArray populations to determine the spike timings of the
# pre and post neuron
inpPop1 = SpikeSourceArray(spike_times=spike_times1)
inpPop2 = SpikeSourceArray(spike_times=spike_times2)
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

projN1_N2 = Projection(
    pre=popN1,
    post=popN2,
    target='Exc',
    synapse=ffSyn
).connect_one_to_one(weights = 0.01)

#---- parameter adjustments ----#
#projN1_N2.thetaLTP = -50.3
projN1_N2.vmean = 80.0
#projN1_N2.transmit = 1.0 # to activate the transmission over the synapse
#------------------------------main function------------------------------------
def run():

    compile()

    #------- neuron Monitors --------#
    # monitor to save changes in the synapse
    dendrite = projN1_N2.dendrite(0)
    m_d = Monitor(dendrite, ['w','deltaW'])

    dWSpk_pos = np.zeros(3)
    deltaWSpk_pos = np.zeros(3)
    for i in range(3):
        projN1_N2.w = initW1
        inpPop2.spike_times = np.linspace(t_1+d_t,t_1+d_t+20*(i),i+1).tolist()
        simulate(duration)
        d_w = m_d.get('w')
        dWSpk_pos[i] = d_w[-1]#np.mean(d_w)
        delta_w = m_d.get('deltaW')        
        deltaWSpk_pos[i] = delta_w[-1]#np.mean(delta_w)
        reset(populations=True,projections=True,synapses=True)
    print(dWSpk_pos/initW1*100)
    dWSpk_neg = np.zeros(3)
    deltaWSpk_neg = np.zeros(3)
    for i in range(3):
        projN1_N2.w = initW1
        inpPop2.spike_times = np.linspace(t_1-d_t,t_1-(d_t+20*(i)),i+1).tolist()
        simulate(duration)
        d_w = m_d.get('w')
        dWSpk_neg[i] = np.mean(d_w)
        delta_w = m_d.get('deltaW')
        deltaWSpk_neg[i] = np.mean(delta_w)
        reset(populations=True,projections=True,synapses=True)
  

    ####################################################
    n_freq = 10
    rates = np.linspace(20,100,n_freq)
    dWBurst_pos = np.zeros(len(rates))
    deltaWBurst_pos = np.zeros(len(rates))
    for i in range(n_freq):
        d_t2 = (1000./rates[i])
        projN1_N2.w = initW2
        inpPop2.spike_times = np.linspace(t_1+d_t,t_1+d_t+(d_t2*2),3).tolist()
        simulate(duration)
        d_w = m_d.get('w')
        dWBurst_pos[i] = d_w[-1] #np.mean(d_w)
        delta_w = m_d.get('deltaW')        
        deltaWBurst_pos[i] = delta_w[-1]#np.mean(delta_w)
        reset(populations=True,projections=True,synapses=True)


    dWBurst_neg = np.zeros(len(rates))
    deltaWBurst_neg = np.zeros(len(rates))
    for i in range(n_freq):
        d_t2 = (1000./rates[i])
        projN1_N2.w = initW2
        inpPop2.spike_times = np.linspace(t_1-d_t,t_1-(d_t+d_t2*2),3).tolist()
        simulate(duration)
        d_w = m_d.get('w')
        dWBurst_neg[i] = np.mean(d_w)
        delta_w = m_d.get('deltaW')
        deltaWBurst_neg[i] = np.mean(delta_w)
        reset(populations=True,projections=True,synapses=True)
    
    #########################################################
    lags = np.linspace(-100,60,33)
    print(lags)
    dWLag_pos = np.zeros(len(lags))
    deltaWLag_pos = np.zeros(len(lags))
    for i in range(len(lags)):
        d_t2 = 20
        projN1_N2.w = initW3
        inpPop2.spike_times = np.linspace(t_1+lags[i],t_1+lags[i]+(d_t2*2),3).tolist()
        simulate(duration)
        d_w = m_d.get('w')
        dWLag_pos[i] = d_w[-1]#np.mean(d_w)
        delta_w = m_d.get('deltaW')        
        deltaWLag_pos[i] = delta_w[-1] #np.mean(delta_w)
        reset(populations=True,projections=True,synapses=True)

    dWSpk_pos = dWSpk_pos/initW1 * 100
    dWSpk_pos[dWSpk_pos>250] = 250
    dWSpk_neg = dWSpk_neg/initW1 * 100

    dWBurst_pos = dWBurst_pos/initW2 *100
    dWBurst_pos[dWBurst_pos>250] = 250
    dWBurst_neg = dWBurst_neg/initW2 *100

    dWLag_pos = dWLag_pos/initW3 *100

    #---plot data-----#
    fig = plt.figure(figsize=(12,10))
    gs=GridSpec(6,4)
    ax0= plt.subplot(gs[0:4,0:2])
    ax0.plot(dWSpk_pos,'x',color='black',lw=3,ms=15)
    ax0.plot(dWSpk_neg,'x',color='black',lw=3,ms=15)
    ax0.hlines(100,-0.2,3,colors='k')
    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)
    ax0.xaxis.set_ticks_position('bottom')
    ax0.yaxis.set_ticks_position('left')
    plt.xticks(np.linspace(0,2,3),np.linspace(1,3,3))
    plt.xlim(-0.2,2.2)
    plt.xlabel('Number of Spikes')
    plt.ylabel('Normalized weight (%)')
    plt.ylim(0.0,250)

    ax1= plt.subplot(gs[0:4,2:4])
    ax1.plot(dWBurst_pos,'--',color='black',lw=3,ms=15)
    ax1.plot(dWBurst_neg,'--',color='black',lw=3,ms=15)
    ax1.hlines(100,0,n_freq,colors='k')
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    plt.xticks(np.linspace(0,n_freq-1,5),np.linspace(20,100,5))
    plt.ylim(0.0,250)
    plt.xlabel('Frequency (Hz)')

    ax2 = plt.subplot(gs[5:6,:])
    ax2.plot(dWLag_pos,'-',color='black',lw=3)
    ax2.hlines(100,0,33,colors='k')
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    plt.ylim(0.0,300)
    fig.savefig('burst.png',bbox_inches='tight')
    print("finish")
#------------------------------------------------------------------------------------
run()
