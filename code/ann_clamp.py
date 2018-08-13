#----------------------imports and environment---------------------------------
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
from ANNarchy import *
import numpy as np
from net_fix import *

"""
Python script to reproduce the voltage clamp experiment.
The membrane potential is fixed to values between -80 mV and 0 mV
and the presynaptic neuron is active with a firing rate of 25 Hz for 50 s.
As in the original publication, the experiment is done to fit
data in the visual cortex and the hippocampus.
See Figure 1. h in Clopath et al. (2010).
"""

###global parameter###
duration = 50*1000 #ms
initW = 0.0001
#----------------------defint time points of spikes-----------------------------#
# create a list of spiking time points. 25 Hz for 50 seconds :
spike_times1 =np.asarray(range(0,duration,1000/4))
print(len(spike_times1))
#-----------------------population defintions-----------------------------------#
"""
To control the spike timings of the AdEx neurons, an additional input population
is used. The spike timing of the SpikeSourceArray can be determined with a
list of time points. """
inpPop1 = SpikeSourceArray(spike_times=spike_times1.tolist())
popN1 = Population(geometry=1,neuron=AdExNeuron, name="N1")
popN2 = Population(geometry=1,neuron=AdExNeuron, name="N2")
#------------------------------------------------------------------------------#
"""
Create a own equation for the learning as mentioned by Clopath et. al. (2010).
For more information about define the learning rule and synapses see
'net_fix.py' or 'net_homeostatic.py'.
"""
equatSTDP_clamp = """
    ltdTerm = if w>wMin : (aLTD*pre.Spike * pos(u_clamp - thetaLTD)) else : 0.0
    ltpTerm = if w<wMax : (aLTP * pos(u_clamp - thetaLTP) *(pre.xtrace)* pos(u_clamp - thetaLTD)) else : 0.0
      deltaW = ( -ltdTerm + ltpTerm)
        dw/dt = deltaW*0.0 :min=0.0,explicite"""

parameter_clamp="""
thetaLTD = -70.6
thetaLTP = -45.3
aLTD = 0.00014
aLTP = 0.00008
wMin = 0.0
wMax =3.0
u_clamp= -80.0
"""

Syn_clamp = Synapse( parameters = parameter_clamp,
    equations= equatSTDP_clamp,
    pre_spike='''g_target += w''')
#-----------------------projection definitions----------------------------------
projST1_C1 = Projection(
    pre=inpPop1,
    post=popN1,
    target='Exc'
).connect_one_to_one(weights = 30.0)

projC1_C2 = Projection(
    pre=popN1,
    post=popN2,
    target='Exc',
    synapse=Syn_clamp
).connect_one_to_one(weights = initW)

#---- parameter adjustments ----#
projC1_C2.transmit = 1.0 # to activate the transmission over the synapse
#------------------------------main function------------------------------------
def run():

    compile()
    # create a list of 100 values from -80 mV to 0 mV for the postsynaptic
    # membrane potential
    post_memb = np.linspace(-80,0,100)

    #------- neuron Monitors --------#
    # monitor to save changes in the synapse
    dendrite = projC1_C2.dendrite(0)
    m_d = Monitor(dendrite, ['deltaW'])#,period=duration)

    rec_dW_norm = np.zeros(len(post_memb))
    rec_W_norm = np.zeros(len(post_memb))

    print('Start voltage clamp experiment.')
    """
    Voltage clamp experiment with the parameter set for the visual cortex
    (standard parameter set).
    """
    for i in range(len(post_memb)):
        projC1_C2.u_clamp=post_memb[i]
        projC1_C2.w = initW
        simulate(duration)
        delta_w = m_d.get('deltaW')
        rec_dW_norm[i] = np.mean(delta_w)
        rec_W_norm[i] = np.mean(projC1_C2.w)+rec_dW_norm[i]
        reset()


    """
    Voltage clamp experiment with the parameter set for the hippocampal
    """
    projC1_C2.thetaLTD = -41.0
    projC1_C2.thetaLTP =-38.0
    projC1_C2.aLTD = 0.00038
    projC1_C2.aLTP = 0.00002

    rec_dW_hippo = np.zeros(len(post_memb))
    rec_W_hippo = np.zeros(len(post_memb))
    spike_times1 =np.asarray(range(0,duration,1000/10))
    inpPop1.spike_times = spike_times1.tolist()

    for i in range(len(post_memb)):
        projC1_C2.u_clamp=post_memb[i]
        projC1_C2.w = initW
        simulate(duration)
        delta_w = m_d.get('deltaW')
        rec_dW_hippo[i] = np.mean(delta_w)#
        rec_W_hippo[i] = np.mean(projC1_C2.w)+rec_dW_hippo[i]
        reset()

    rec_W_norm[rec_W_norm>np.max(rec_W_hippo)] = np.max(rec_W_hippo)

    rec_W_norm = rec_W_norm/initW*100.
    rec_W_hippo = rec_W_hippo/initW*100.

    ##--- start plotting ---##

    # get the index of the theta values to draw the lines in the final plot
    theta_m = -70.6
    theta_pl= -45.3
    ix1 = np.where(post_memb<theta_m)
    ixM = np.argmax(post_memb[ix1])
    ix1 = np.where(post_memb<theta_pl)
    ixPL = np.argmax(post_memb[ix1])

    fig,ax = plt.subplots()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    plt.plot(rec_W_norm,'--',color='steelblue',lw=3.0)
    plt.plot(rec_W_hippo,color='tomato',lw=3.0)
    plt.xlabel('Voltage (mv)')
    #plt.ylim(ymin=0.00005,ymax=0.0003)
    #plt.yticks(np.linspace(0.00005,0.0003,6),np.linspace(50,300,6))
    #plt.xlim(0,20)
    plt.ylabel('Normalized weight (%)')
    plt.xticks(np.linspace(0,len(post_memb)-1,5),np.linspace(-80,0,5))
    plt.axvline(ixM,color='k', linestyle='--')
    plt.text(ixM,np.max(rec_W_hippo),r'$\theta_{-}$',fontsize=20)
    plt.axvline(ixPL,color='k', linestyle='--')
    plt.text(ixPL,np.max(rec_W_hippo),r'$\theta_{+}$',fontsize=20)
    plt.savefig('W_hippo',bbox_inches='tight')
    print("finish")
#------------------------------------------------------------------------------------
run()
