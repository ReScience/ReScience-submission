"""
Python code to reproduce the pairing repetition task in Clopath et al. 2010 (Fig. 2 b).
Pairs of pre-post and post-pre spikes for different pairing repetition frequencies.
The original experiment is from Sjoestroem et al. 2001.
Between the spikes of each pair elapse 10 ms.
For the correct timing, every neuron receive input from a extra neuron,
which spikes to a certain time point.
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from ANNarchy import *
setup(dt=1)
from network import *

# Global parameters
duration = 1000 #ms == 1 s

# Spike times
spike_times1 =[[0]]
spike_times2 =[[10]]

# Populations
"""
To control the spike timings of the AdEx neurons, two additional input populations
are used. The spike timing of the SpikeSourceArray can be determined with a
list of time points.
"""
inpPop1 = SpikeSourceArray(spike_times=spike_times1)
inpPop2 = SpikeSourceArray(spike_times=spike_times2)
popN1 = Population(geometry=1,neuron=AdExNeuron, name="N1")
popN2 = Population(geometry=1,neuron=AdExNeuron, name="N2")

# Projections
"""
Define simple projection from the input SpikeSourceArray populations
to the neuron populations.
If the neuron in the input population spikes,
1 ms later a spike in the connected AdEx neuron population is triggered.
"""

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

"""
Create the projection between the two AdEx neurons
"""
projV1_V1 = Projection(
    pre=popN1,
    post=popN2,
    target='Exc',
    synapse=ffSyn
).connect_one_to_one(weights = 0.1)



projV1_V1.vmean_fix = 100.0
projV1_V1.set_fix = 1.0 # use a fix apmlitude for the LTD term

def run():
    print('Start experiment to reproduce the pairing frequency data.')
    # Compile command to create the ANNarchy network
    compile()

    # Create a single dendrite object to record the weight of this dendrite
    dendrite = projV1_V1.dendrite(0)
    m_d = Monitor(dendrite, ['w','deltaW','ltdTerm','ltpTerm'])

    # Set max repetition frequency
    max_freq = 50

    # Time between a pre and a post spike (or post and pre spike)
    td = 10#ms

    # Inital weight value
    initW = 0.125

    # Save the weight change (dw) for pre post spike pairs
    dW_prePost =[]
    for f in np.arange(0.1,max_freq):
        # Reset the network
        reset()
        projV1_V1.w = initW
        spike_times1 = np.linspace(0, int(duration), int(f+1))
        spike_times2 = np.linspace(int(0+td), int(duration+td), int(f+1))
        # Set the spike times with the actual repetition frequency f
        inpPop1.spike_times = spike_times1.tolist()
        inpPop2.spike_times = spike_times2.tolist()
        simulate(duration)
        # Save records
        dW_prePost.append(m_d.get('w'))

    # Save the weight change (dw) for post pre spike pairs
    dW_postPre =[]
    for f in np.arange(0.1,max_freq):
        # Reset the network
        reset()
        projV1_V1.w = initW
        spike_times1 = np.linspace(20, int(duration-30), int(f+1))
        spike_times2 = np.linspace(int(20-td), int(duration-td-30), int(f+1))
        # Set the spike times with the actual repetition frequency f
        inpPop1.spike_times = spike_times1.tolist()
        inpPop2.spike_times = spike_times2.tolist()
        simulate(duration)
        # Save records
        dW_postPre.append(m_d.get('w'))

    # Estimate the total change per repetition frequency
    sumdW_prePost = np.zeros(max_freq)
    sumdW_postPre = np.zeros(max_freq)
    for f in range(len(dW_postPre)):
        sumdW_prePost[f] = np.mean(dW_prePost[f])
        sumdW_postPre[f] = np.mean(dW_postPre[f])

    # Start plotting
    fig,ax = plt.subplots()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    plt.plot(sumdW_prePost,color='steelblue',lw=6)
    plt.plot(sumdW_postPre,color='tomato',lw=6)

    upB = initW/100.* 150.
    loB = initW/100. * 50.0


    plt.xlabel(r'$\rho$ [Hz]',fontsize=25)
    plt.ylabel('Normalized weight (%)',fontsize=25)
    plt.xlim(0.0,50.0)
    plt.yticks(np.linspace(loB,upB,3),range(50,200,50),fontsize=20 )
    plt.xticks(fontsize=20)
    plt.savefig('Fig1_pairing.png',bbox_inches='tight', pad_inches = 0.1)
    plt.show()
    print("Done with the experiment.")

if __name__ == "__main__":
    run()
