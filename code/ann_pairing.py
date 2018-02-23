#----------------------imports and environment---------------------------------
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
from ANNarchy import *
import numpy as np
from net import *

"""
Python code to reproduce the pairing repition task in Clopath et al. 2010.
Pairs of pre-post and post-pre spikes for different pairing repition frequencies.
The original experiment is from Sjoestroem et al. 2001.
Between the spikes of each pair elapse 10 ms.
For the correct timing, every neuron recive input from a extra neuron,
which spikes to a certain time point.
"""

###global parameter###
duration = 1000 #ms
#----------------------defint time points of spikes-----------------------------#
spike_times1 =[[0]]
spike_times2 =[[10]]
#-----------------------population defintions-----------------------------------#
inpPop1 = SpikeSourceArray(spike_times=spike_times1)
inpPop2 = SpikeSourceArray(spike_times=spike_times2)#PoissonPopulation(geometry=1,rates=20)
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
projV1_V1.vmean = 120.0

#------------------------------main function------------------------------------
def run():

    compile()

    #------- neuron Monitors --------#
    # create a single dendrite object to record the weight of this dendrite
    dendrite = projV1_V1.dendrite(0)
    m_d = Monitor(dendrite, ['w','deltaW','ltdTerm','ltpTerm'])

    # set max repetition frequency
    max_freq = 50
    # time between a pre and a post spike (or post and pre spike)
    td = 10#ms

    defW = 0.125

    # save the weight change (dw) for pre post spike pairs
    dW_prePost =[]
    for f in np.arange(0.1,max_freq):
        #reset the network#
        reset()
        projV1_V1.w = defW
        spike_times1 = np.linspace(0,duration,f+1)
        spike_times2 = np.linspace(0+td,duration+td,f+1)
        # set the spike times with the actual repetition frequency f
        inpPop1.spike_times = spike_times1.tolist()
        inpPop2.spike_times = spike_times2.tolist()
        simulate(duration)
        #save records#
        dW_prePost.append(m_d.get('w'))

    # save the weight change (dw) for post pre spike pairs
    dW_postPre =[]
    for f in np.arange(0.1,max_freq):
        #reset the network#
        reset()
        projV1_V1.w = defW
        spike_times1 = np.linspace(20,duration-30,f+1)
        spike_times2 = np.linspace(20-td,duration-td-30,f+1)
        # set the spike times with the actual repetition frequency f
        inpPop1.spike_times = spike_times1.tolist()
        inpPop2.spike_times = spike_times2.tolist()
        simulate(duration)
        #save records#
        dW_postPre.append(m_d.get('w'))

    # estimate the total change per repition frequency
    sumdW_prePost = np.zeros(max_freq)
    sumdW_postPre = np.zeros(max_freq)
    for f in range(len(dW_postPre)):
        sumdW_prePost[f] = np.mean(dW_prePost[f])
        sumdW_postPre[f] = np.mean(dW_postPre[f])

    # plot the figure 2b in Clopath et al. 2010
    fig,ax = plt.subplots()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    plt.plot(sumdW_prePost,color='steelblue',lw=3)
    plt.plot(sumdW_postPre,color='tomato',lw=3)

    upB = defW/100.* 150.
    loB = defW/100. * 50.0


    plt.xlabel(r'$\rho$ [Hz]',fontsize=12)
    plt.ylabel('Normalized weight (%)',fontsize=12)
    plt.xlim(0.0,50.0)
    plt.yticks(np.linspace(loB,upB,3),range(50,200,50) )
    plt.savefig('pairing.png',bbox_inches='tight', pad_inches = 0.1)

    print("finish")
#------------------------------------------------------------------------------------
run()
