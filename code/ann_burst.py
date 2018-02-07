#----------------------imports and environment---------------------------------
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
from ANNarchy import *
import numpy as np
from net import *

###global parameter###
duration = 130 #ms
#----------------------defint time points of spikes-----------------------------#
spike_times1 =[[10]]
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
).connect_one_to_one(weights = 1.0)
projV1_V1.vmean=100.0
#------------------------------main function------------------------------------
def run():

    compile()

    ##set parameters##
#    ffSyn.aLTD = 0.0014
#    ffSyn.aLTP = 0.0008

    #------- neuron Monitors --------#
    m_N1 = Monitor(popN1,['spike'])
    m_N2 = Monitor(popN2,['spike'])    
    m_ST1 = Monitor(inpPop1,['spike'])
    m_ST2 = Monitor(inpPop2,['spike'])
    dendrite = projV1_V1.dendrite(0)
    m_d = Monitor(dendrite, ['w','deltaW','ltdTerm','ltpTerm'])

    spks_N1 = []
    spks_N2 = []
    dW_Burst =[]
    #--set the spike times for first test--#
    for i in xrange(3):
        spike_times1 = [10]
        spike_times2 = np.linspace(20,20+(i*20),i+1)
        inpPop1.spike_times = spike_times1
        inpPop2.spike_times = spike_times2.tolist()
        simulate(duration)
        #save records#
        spk_N1 = m_N1.get('spike')
        t_n1,n_1 = m_N1.raster_plot(spk_N1)      
        spks_N1.append(t_n1)
        spk_N2 = m_N2.get('spike')
        t_n2,n_2 = m_N2.raster_plot(spk_N2)      
        spks_N2.append(t_n2)  
        dW_Burst.append(m_d.get('deltaW'))
        #reset the network#
        reset()
        projV1_V1.w=1.0


    dW_Burst_2 =[]
    #--set the spike times for first test--#
    for i in xrange(3):
        spike_times1 = [60]
        spike_times2 = np.linspace(50-(i*20),50,i+1)
        inpPop1.spike_times = spike_times1
        inpPop2.spike_times = spike_times2.tolist()
        simulate(duration)
        #save records#  
        dW_Burst_2.append(m_d.get('deltaW'))
        #reset the network#
        reset()
        projV1_V1.w=1.0

    meandW_Burst = np.zeros(3)
    meandW_Burst_2 = np.zeros(3)
    for i in xrange(3):
        meandW_Burst[i] = np.sum(dW_Burst[i])
        meandW_Burst_2[i] = np.sum(dW_Burst_2[i])

    plt.figure()
    plt.plot(meandW_Burst,'o',color='steelblue')
    plt.plot(meandW_Burst_2,'o',color='tomato')
    plt.savefig('burst_1.png')

####------- task2 ------####
    dW_Burst =[]
    dW_Burst_1 =[]
    #--set the spike times for first test--#
    ##--20 Hz--##
    spike_times1 = [10]
    spike_times2 = np.linspace(20,120,3)
    inpPop1.spike_times = spike_times1
    inpPop2.spike_times = spike_times2.tolist()
    simulate(duration)
    #save records#
    dW_Burst.append(m_d.get('deltaW'))
    #reset the network#
    reset()
    projV1_V1.w=1.0
    ##--50 Hz--##
    spike_times1 = [10]
    spike_times2 = np.linspace(20,60,3)
    inpPop1.spike_times = spike_times1
    inpPop2.spike_times = spike_times2.tolist()
    simulate(duration)
    #save records#
    dW_Burst.append(m_d.get('deltaW'))
    #reset the network#
    reset()
    projV1_V1.w=1.0
    ##--10 Hz--##
    spike_times1 = [10]
    spike_times2 = np.linspace(20,40,3)
    inpPop1.spike_times = spike_times1
    inpPop2.spike_times = spike_times2.tolist()
    simulate(duration)
    #save records#
    dW_Burst.append(m_d.get('deltaW'))
    #reset the network#
    reset()
    projV1_V1.w=1.0

    #--set the spike times for first test--#
    ##--20 Hz--##
    spike_times1 = [120]
    spike_times2 = np.linspace(10,110,3)
    inpPop1.spike_times = spike_times1
    inpPop2.spike_times = spike_times2.tolist()
    simulate(duration)
    #save records#
    dW_Burst_1.append(m_d.get('deltaW'))
    #reset the network#
    reset()
    projV1_V1.w=1.0
    ##--50 Hz--##
    spike_times1 = [120]
    spike_times2 = np.linspace(70,110,3)
    inpPop1.spike_times = spike_times1
    inpPop2.spike_times = spike_times2.tolist()
    simulate(duration)
    #save records#
    dW_Burst_1.append(m_d.get('deltaW'))
    #reset the network#
    reset()
    projV1_V1.w=1.0
    ##--10 Hz--##
    spike_times1 = [120]
    spike_times2 = np.linspace(90,110,3)
    inpPop1.spike_times = spike_times1
    inpPop2.spike_times = spike_times2.tolist()
    simulate(duration)
    #save records#
    dW_Burst_1.append(m_d.get('deltaW'))
    #reset the network#
    reset()
    projV1_V1.w=1.0

    meandW_Burst = np.zeros(3)
    meandW_Burst_1 = np.zeros(3)
    for i in xrange(3):
        meandW_Burst[i] = np.sum(dW_Burst[i])
        meandW_Burst_1[i] = np.sum(dW_Burst_1[i])

    plt.figure()
    plt.plot(meandW_Burst,'o',color='steelblue')
    plt.plot(meandW_Burst_1,'o',color='tomato')
    plt.savefig('burst_2.png')
    print("finish")
#------------------------------------------------------------------------------------
run()
