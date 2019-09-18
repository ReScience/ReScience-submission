"""
Python script to reproduce the STDP window protocol.
Record the change in the synaptic weight for different time intervals
between pre- and postsynaptic spike. See Fig. 2 a in Clopath et al. (2010).
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from ANNarchy import *
setup(dt=1)
from network import *

# Global parameters
duration = 40 # duration time of 40 ms
initW = 0.012 # initial weight

# Spike times
spike_times1 =[[0]]
spike_times2 =[[10]]

# Population definitions
"""
To control the spike timings of the AdEx neurons, two additional input populations
are used. The spike timing of the SpikeSourceArray can be determined with a
list of time points. """
inpPop1 = SpikeSourceArray(spike_times=spike_times1)
inpPop2 = SpikeSourceArray(spike_times=spike_times2)
popN1 = Population(geometry=1,neuron=AdExNeuron, name="N1")
popN2 = Population(geometry=1,neuron=AdExNeuron, name="N2")

# Projection definitions
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
).connect_one_to_one(weights = initW)

# Parameter adjustments
projV1_V1.vmean_fix = 70.0
projV1_V1.set_fix = 1.0 # use a fix apmlitude for the LTD term

def run():
    "Runs the STDP window protocol"
    print('Start the experiment to reproduce the STDP window')
    # 31 spiking pairs for a time difference dt between a pre and post
    # synaptic spike from -15 to 15 ms
    dt = np.linspace(-15,15,31)
    compile()

    # Monitor to save changes in the synapse
    dendrite = projV1_V1.dendrite(0)
    m_d = Monitor(dendrite, ['w','deltaW','ltdTerm','ltpTerm'])
    n_pairs = 31

    w = np.zeros(n_pairs)
    dW = np.zeros(n_pairs)
    for i in range(n_pairs):
        #reset the network#
        reset()
        projV1_V1.w = initW # set weight back to initial weight value

        inpPop1.spike_times = [16] # presynaptic neuron always spikes at t=16 ms

        # add the time difference to estimate the postsynaptic spike time
        inpPop2.spike_times = [16+dt[i]]
        simulate(duration)
        d_w = m_d.get('w')
        delta_w = m_d.get('deltaW')
        w[i] = d_w[-1]
        dW[i] = np.sum(delta_w)

    # Get the recorded data
    ltd_w = m_d.get('ltdTerm')
    ltp_w = m_d.get('ltpTerm')

    w = (initW+dW)/initW*100.

    # Start plotting
    fig,ax = plt.subplots(figsize=(13,9))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    plt.plot(np.linspace(0,13,14),w[0:14],color='tomato',lw=10.0)
    plt.plot(np.linspace(13,17,5),np.linspace(w[13],w[17],5),color='black',lw=2)
    plt.plot(np.linspace(17,30,14),w[17:31],color='steelblue',lw=10.0)

    plt.axhline(y=100.0, color='k',linestyle='--')
    plt.axvline(x=15,color='k',linestyle='--')
    plt.ylim(ymin=55,ymax=145)
    plt.ylabel('Normalized weight (%)',fontsize=30)
    plt.xticks(np.linspace(5,25,3),np.linspace(-10,10,3),fontsize=25)
    plt.yticks(fontsize=25)
    plt.xlabel('T (ms)',fontsize=30)

    plt.savefig('Fig1_window.png',bbox_inches='tight', pad_inches = 0.1)
    plt.show()
    print("Done with the experiment.")

if __name__ == "__main__":
    run()
