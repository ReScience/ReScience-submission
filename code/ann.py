#----------------------imports and environment---------------------------------
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
from ANNarchy import *
import numpy as np
from net import *

#----------------------defint time points of spikes-----------------------------#
spike_times1 =[[49]] #[[10,20,30]]
spike_times2 =[[47]] #[[40,50,60]]
#-----------------------population defintions-----------------------------------#
inpPop1 = PoissonPopulation(geometry=1,rates=20)#SpikeSourceArray(spike_times=spike_times1)
inpPop2 = PoissonPopulation(geometry=1,rates=20)#SpikeSourceArray(spike_times=spike_times2)
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

#------------------------------main function------------------------------------
def run():

    compile()
    duration = 4000 #ms
    #------- neuron Monitors --------#
    m_N1 = Monitor(popN1,['vm','vmean','umeanLTD','umeanLTP','spike','xtrace','Spike'])
    m_N2 = Monitor(popN2,['vm','vmean','umeanLTD','umeanLTP','spike','xtrace','Spike'])
    m_ST1 = Monitor(inpPop1,['spike'])
    m_ST2 = Monitor(inpPop2,['spike'])
    dendrite = projV1_V1.dendrite(0)
    m_d = Monitor(dendrite, ['w','deltaW','ltdTerm','ltpTerm'])

    #--do not send any excitaory potential during the two V1 neurons--#
    #projV1_V1.transmission = False
    simulate(duration)


    #---get the recorded data-----#
    spk_N1 = m_N1.get('spike')
    t_n1,n_1 = m_N1.raster_plot(spk_N1)
    vm_N1 = m_N1.get('vm')
    vmean_N1=m_N1.get('vmean')
    vmeanLTD_N1=m_N1.get('umeanLTD')
    xtrace_N1 = m_N1.get('xtrace')
    spike_N1 = m_N1.get('Spike')

    spk_N2 = m_N2.get('spike')
    t_n2,n_2 = m_N2.raster_plot(spk_N2)
    vm_N2 = m_N2.get('vm')
    vmean_N2=m_N2.get('vmean')
    vmeanLTD_N2=m_N2.get('umeanLTD')
    xtrace_N2 = m_N2.get('xtrace')
    spike_N2 = m_N2.get('Spike')

    spk_ST1 = m_ST1.get('spike')
    t_st1,n_st1 = m_ST1.raster_plot(spk_ST1)

    spk_ST2 = m_ST2.get('spike')
    t_st2,n_st2 = m_ST2.raster_plot(spk_ST2)    


    d_w = m_d.get('w')
    delta_w = m_d.get('deltaW')
    ltd_w = m_d.get('ltdTerm')
    ltp_w = m_d.get('ltpTerm')

    #--- sort the different spike times of the pre and post--#
    t_wind = 15 #ms
    for pre_spkt in t_st1:
        idx = np.where((t_st2 >= (pre_spkt-t_wind))&(t_st2 <= (pre_spkt+t_wind)))
        if len(idx[0]) >0:
            print(pre_spkt)
            print(t_st2[idx[0]])
            print('-------------')
    #----------plot simulation results---------------------#
    plt.figure(figsize=(10,10))
    plt.subplot(221)
    plt.plot(t_n1,n_1,'o')
    plt.title('Neuron 1')
    plt.xlabel('Time (ms)')
    plt.ylabel('# neuron')
    plt.xlim(0.0,duration)
    plt.subplot(222)
    plt.plot(t_n2,n_2,'o')
    plt.title('Neuron 2')
    plt.xlim(0.0,duration)
    plt.subplot(223)
    plt.plot(t_st1,n_st1,'o')
    plt.title('Input 1')
    plt.xlim(0.0,duration)
    plt.subplot(224)
    plt.plot(t_st2,n_st2,'o')
    plt.title('Input 2')
    plt.xlim(0.0,duration)
    plt.savefig('spike_test.png')

    plt.figure(figsize=(10,20))
    plt.subplot(321)
    plt.plot(vm_N1)
    plt.title('Neuron 1')
    plt.ylabel('vm')
    plt.subplot(322)
    plt.plot(vm_N2)
    plt.title('Neuron 2')
    plt.subplot(323)
    plt.plot(vmeanLTD_N1)
    plt.ylabel('vmmean')
    plt.subplot(324)
    plt.plot(vmeanLTD_N2)
    plt.subplot(325)
    plt.plot(spike_N1)
    plt.ylabel('x trace')
    plt.xlabel('Time (ms)')
    plt.subplot(326)
    plt.plot(spike_N2)
    plt.xlabel('Time (ms)')
    plt.savefig('vm_test.png')


    plt.figure(figsize=(12,12))
    plt.subplot(221)
    plt.plot(d_w)
    plt.xlabel('time [ms]')
    plt.ylabel('weight')
    plt.subplot(222)
    plt.plot(delta_w)
    plt.xlabel('time [ms]')
    plt.ylabel('delta weight')
    plt.subplot(223)
    plt.plot(ltd_w)
    plt.ylabel('ltd')
    plt.subplot(224)
    plt.plot(ltp_w)
    plt.ylabel('ltp')
    plt.savefig('weight.png')
    print("finish")
#------------------------------------------------------------------------------------
run()
