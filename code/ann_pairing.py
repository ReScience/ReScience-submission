#----------------------imports and environment---------------------------------
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
from ANNarchy import *
import numpy as np
from net import *

###global parameter###
duration = 1040 #ms
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
).connect_one_to_one(weights = 0.001)
projV1_V1.vmean = 80.0

#------------------------------main function------------------------------------
def run():

    compile()

    

    #------- neuron Monitors --------#
    m_N1 = Monitor(popN1,['spike'])
    m_N2 = Monitor(popN2,['spike'])    
    m_ST1 = Monitor(inpPop1,['spike'])
    m_ST2 = Monitor(inpPop2,['spike'])
    dendrite = projV1_V1.dendrite(0)
    m_d = Monitor(dendrite, ['w','deltaW','ltdTerm','ltpTerm'])

    #--do not send any excitaory potential during the two V1 neurons--#
    max_freq = 51#50
    td = 10#ms
    t_spk1 = []
    t_spk2 = []
    dW_prePost =[]
    for f in range(max_freq):
        #reset the network#
        reset()
        projV1_V1.w=1.0
        spike_times1 = np.linspace(0,duration-30,f+1)
        spike_times2 = np.linspace(0+td,duration+td-30,f+1)    
        inpPop1.spike_times = spike_times1.tolist()
        inpPop2.spike_times = spike_times2.tolist()
        simulate(duration)
        #save records#
        spk_N1 = m_N1.get('spike')
        t_n1,n_1 = m_N1.raster_plot(spk_N1)
        t_spk1.append(t_n1)
        spk_N2 = m_N2.get('spike')
        t_n2,n_2 = m_N2.raster_plot(spk_N2)
        t_spk2.append(t_n2)
        dW_prePost.append(m_d.get('deltaW'))

       
    dW_postPre =[]
    for f in range(max_freq):
        #reset the network#
        reset()
        projV1_V1.w=1.0
        spike_times1 = np.linspace(20,duration-30,f+1)
        spike_times2 = np.linspace(20-td,duration-td-30,f+1)
        inpPop1.spike_times = spike_times1.tolist()
        inpPop2.spike_times = spike_times2.tolist()
        simulate(duration)
        #save records#
        spk_N1 = m_N1.get('spike')
        t_n1,n_1 = m_N1.raster_plot(spk_N1)
        t_spk1.append(t_n1)
        spk_N2 = m_N2.get('spike')
        t_n2,n_2 = m_N2.raster_plot(spk_N2)
        t_spk2.append(t_n2)
        dW_postPre.append(m_d.get('deltaW'))


    meandW_prePost = np.zeros(max_freq)
    meandW_postPre = np.zeros(max_freq)
    for f in range(max_freq):
        meandW_prePost[f] = np.mean(dW_prePost[f])
        meandW_postPre[f] = np.mean(dW_postPre[f])

    fig,ax = plt.subplots()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    plt.plot(meandW_prePost[0:-1:3],color='steelblue',lw=3)
    plt.plot(meandW_postPre[0:-1:3],color='tomato',lw=3)
    plt.savefig('pairing.png',bbox_inches='tight', pad_inches = 0.1)

    #--- sort the different spike times of the pre and post--#

    print("finish")
#------------------------------------------------------------------------------------
run()
