# Script allowing to reproduce Fig. 3 of:
#
#   Laje, R. and Buonomano, D.V. (2013). Robust timing and motor patterns by taming chaos in recurrent neural networks. Nat Neurosci.
#
# Author: Julien Vitay (julien.vitay@informatik.tu-chemnitz.de)
# Licence: MIT
from __future__ import print_function
import numpy as np
import scipy.stats
import time

# Import the definition of the network
from RecurrentNetwork import RecurrentNetwork

###############
# Parameters
###############
nb_learning_trials_rec = 20 # Number of learning trials for the recurrent weights
nb_learning_trials_readout = 10 # Number of learning trials for the readout weights
nb_networks = 10 # Number of different networks

stimulus_amplitude = 5.0 # Amplitude of the input pulse
t_offset = 200 # Time to wait before the stimulation
d_stim= 50 # Duration of the stimulation
t_relax = 150 # Duration to relax after the trajectory

target_baseline = 0.2 # Baseline of the target function
target_amplitude = 1. # Maximal value of the target function
target_width = 30. # Width of the Gaussian

###################
# Main procedure
###################

# Vary the timing interval
delays = [250, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]

# Store the Pearson correlation coefficients
pearsons = []

# Iterate over the delays
for target_time in delays:
    print('*'*60)
    print('Learning a delay of', target_time)
    print('*'*60)
    d_trajectory = target_time + 150 # Duration of the desired trajectory
    trial_duration = t_offset + d_stim + d_trajectory + t_relax # Total duration of a trial
    pearsons.append([])

    for n in range(nb_networks): # 10 networks per delay
        print('*'*60)
        print('Network', n+1)
        print('*'*60)

        # Create a new network (note: networks are reused in the original article)
        net= RecurrentNetwork(
            Ni = 2, # Number of inputs
            N = 800, # Number of recurrent neurons
            No = 1, # Number of read-out neurons
            tau = 10.0, # Time constant of the neurons
            g = 1.5, # Synaptic strength scaling
            pc = 0.1, # Connection probability
            Io = 0.001, # Noise variance
            delta = 1.0, # Initial diagonal value of the P matrix
            P_plastic = 0.6, # Percentage of neurons receiving plastic synapses
        )

        # Impulse input after 200 ms
        impulse = np.zeros((net.Ni, 1, trial_duration))
        impulse[0, 0, t_offset:t_offset+d_stim] = stimulus_amplitude

        # Target output for learning the readout weights
        target = np.zeros((net.No, 1, trial_duration))
        time_axis = np.linspace(0, trial_duration, trial_duration)
        target[0, 0, : ] = target_baseline + (target_amplitude - target_baseline) * np.exp(-(t_offset + d_stim + target_time - time_axis)**2/target_width**2)

        # Initial trial to determine the innate trajectory
        print('Initial trial to determine a trajectory (without noise)')
        trajectory, initial_output = net.simulate(stimulus=impulse, noise=False)

        # 20 trials of learning for the recurrent weights
        for i in range(nb_learning_trials_rec):
            print('Learning trial recurrent', i+1)
            _, _ = net.simulate(stimulus=impulse, trajectory=trajectory,
                                learn_start=t_offset+d_stim, learn_stop=t_offset+d_stim+d_trajectory)

        # 10 trials of learning for the readout weights
        for i in range(nb_learning_trials_readout):
            print('Learning trial readout', i+1)
            _, _ = net.simulate(stimulus=impulse, trajectory=target,
                                learn_start=t_offset+d_stim, learn_stop=t_offset+d_stim+d_trajectory,
                                learn_readout=True)

        # Test trial
        print('Test trial')
        reproduction, final_output = net.simulate(stimulus=impulse)

        # Pearson correlation coefficient
        pred = final_output[t_offset+d_stim:t_offset+d_stim+d_trajectory, 0, 0]
        desired = target[0, 0, t_offset+d_stim:t_offset+d_stim+d_trajectory ]
        r, p = scipy.stats.pearsonr(desired, pred)
        pearsons[-1].append(r)

# Save the results
np.savez('../data/timingcapacity.npz', r=np.array(pearsons))

##################
# Visualization
##################
correlation_mean = np.mean(pearsons**2, axis=1)
correlation_std = np.std(pearsons**2, axis=1)
plt.errorbar(np.array(delays)/1000., correlation_mean, correlation_std/np.sqrt(10), linestyle='-', marker='^')
plt.xlim((0., 8.5))
plt.ylim((-0.1, 1.1))
plt.xlabel('Interval (s)')
plt.ylabel('Performance ($R^2$)')
plt.show()
