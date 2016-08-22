# Script allowing to reproduce Fig. 2 of:
#
#   Laje, R. and Buonomano, D.V. (2013). Robust timing and motor patterns by taming chaos in recurrent neural networks. Nat Neurosci.
#
# Author: Julien Vitay (julien.vitay@informatik.tu-chemnitz.de)
# Licence: MIT
from __future__ import print_function
import numpy as np
import scipy.io as sio
import time

# Import the definition of the network
from RecurrentNetwork import RecurrentNetwork

###############
# Parameters
###############
nb_learning_trials_rec = 30 # Number of learning trials for the recurrent weights
nb_learning_trials_readout = 10 # Number of learning trials for the readout weights
nb_test_trials = 5 # Number of test trials
nb_perturbation_trials = 5 # Number of perturbation trials

stimulus_amplitude = 2.0 # Amplitude of the input pulse
t_offset = 200 # Time to wait before the stimulation
d_stim= 50 # Duration of the stimulation
t_relax = 150 # Duration to relax after the trajectory

perturbation_amplitude = 0.2 # Amplitude of the perturbation pulse
t_perturbation = 300 # Offset for the perturbation
d_perturbation = 10 # Duration of the perturbation


####################
# Create the network
####################
net = RecurrentNetwork(
    Ni = 4, # Number of inputs
    N = 800, # Number of recurrent neurons
    No = 2, # Number of read-out neurons
    tau = 10.0, # Time constant of the neurons
    g = 1.5, # Synaptic strength scaling
    pc = 0.1, # Connection probability
    Io = 0.001, # Noise variance
    delta = 1.0, # Initial value of the P matrix
    P_plastic = 0.6, # Percentage of neurons receiving plastic synapses
)

###################
# Input definitions
###################
# Retrieve the targets and reformat them
try:
    targets = sio.loadmat('../data/DAC_handwriting_output_targets.mat')
except:
    print("You have to download the handwriting data first.")
    print("Go to the data/ folder and run the script get_handwriting.sh:")
    print("$ bash get_handwriting.sh")
    exit(-1)

chaos = targets['chaos']
neuron = targets['neuron']

# Durations
_, d_chaos = chaos.shape
_, d_neuron = neuron.shape

# Impulses
impulse_chaos = np.zeros((t_offset + d_stim + d_chaos + t_relax, net.Ni, 1))
impulse_chaos[t_offset:t_offset+d_stim, 0, 0] = stimulus_amplitude
impulse_neuron = np.zeros((t_offset + d_stim + d_neuron + t_relax, net.Ni, 1))
impulse_neuron[t_offset:t_offset+d_stim, 2, 0] = stimulus_amplitude

# Perturbation
perturbation_chaos = np.zeros((t_offset + d_stim + d_chaos + t_relax, net.Ni, 1))
perturbation_chaos[t_offset:t_offset+d_stim, 0, 0] = stimulus_amplitude
perturbation_chaos[t_offset + t_perturbation: t_offset + t_perturbation + d_perturbation, 1, 0] = perturbation_amplitude
perturbation_neuron = np.zeros((t_offset + d_stim + d_neuron + t_relax, net.Ni, 1))
perturbation_neuron[t_offset:t_offset+d_stim, 2, 0] = stimulus_amplitude
perturbation_neuron[t_offset + t_perturbation: t_offset + t_perturbation + d_perturbation, 3, 0] = perturbation_amplitude

# Targets
target_chaos = np.zeros((t_offset + d_stim + d_chaos + t_relax, net.No, 1))
target_chaos[t_offset + d_stim: t_offset + d_stim + d_chaos, :, 0] = chaos.T
target_neuron = np.zeros((t_offset + d_stim + d_neuron + t_relax, net.No, 1))
target_neuron[t_offset + d_stim: t_offset + d_stim + d_neuron, :, 0] = neuron.T

###################
# Main procedure
###################
tstart = time.time()

# Initial trial to determine the innate trajectory
print('Initial chaos trial')
trajectory_chaos, initial_chaos_output = net.simulate(stimulus=impulse_chaos, noise=False)
print('Initial neuron trial')
trajectory_neuron, initial_neuron_output = net.simulate(stimulus=impulse_neuron, noise=False)

# 30 trials of learning for the recurrent weights
for i in range(nb_learning_trials_rec):
    print('Learning trial recurrent', i+1)
    _, _ = net.simulate(stimulus=impulse_chaos, trajectory=trajectory_chaos,
                        learn_start=t_offset+d_stim, learn_stop=t_offset+d_stim+d_chaos)
    _, _ = net.simulate(stimulus=impulse_neuron, trajectory=trajectory_neuron,
                        learn_start=t_offset+d_stim, learn_stop=t_offset+d_stim+d_neuron)

# 10 trials of learning for the readout weights
for i in range(nb_learning_trials_readout):
    print('Learning trial readout', i+1)
    _, _ = net.simulate(stimulus=impulse_chaos, trajectory=target_chaos,
                        learn_start=t_offset+d_stim, learn_stop=t_offset+d_stim+d_chaos,
                        learn_readout=True)
    _, _ = net.simulate(stimulus=impulse_neuron, trajectory=target_neuron,
                        learn_start=t_offset+d_stim, learn_stop=t_offset+d_stim+d_neuron,
                        learn_readout=True)

# Save the network at the end of learning
net.save('network-chaosneuron.npz')

# Test trials
final_output_chaos = []
final_output_neuron = []
for _ in range(nb_test_trials):
    print('Test chaos trial')
    _, o = net.simulate(stimulus=impulse_chaos)
    final_output_chaos.append(o)
    print('Test neuron trial')
    _, o = net.simulate(stimulus=impulse_neuron)
    final_output_neuron.append(o)

# Perturbation trials
perturbation_output_chaos = []
perturbation_output_neuron = []
for _ in range(nb_perturbation_trials):
    print('Perturbation chaos trial')
    _, o = net.simulate(stimulus=perturbation_chaos)
    perturbation_output_chaos.append(o)
    print('Perturbation neuron trial')
    _, o = net.simulate(stimulus=perturbation_neuron)
    perturbation_output_neuron.append(o)

print('Simulation done in', time.time() - tstart, 'seconds.')

##################
# Visualization
##################
import matplotlib.pyplot as plt

subsampling_chaos = (t_offset + d_stim + np.linspace(0, d_chaos, 20)).astype(np.int32)
subsampling_neuron = (t_offset + d_stim + np.linspace(0, d_neuron, 20)).astype(np.int32)

ax = plt.subplot2grid((2,2),(0, 0))
ax.plot(chaos[0, :], chaos[1, :], linewidth=2.)
for i in range(nb_perturbation_trials):
    ax.plot(final_output_chaos[i][t_offset + d_stim: t_offset + d_stim + d_chaos, 0, 0],
            final_output_chaos[i][t_offset + d_stim: t_offset + d_stim + d_chaos, 1, 0])
    ax.plot(final_output_chaos[i][subsampling_chaos, 0, 0],
            final_output_chaos[i][subsampling_chaos, 1, 0], 'bo')
ax.set_xlabel('x')
ax.set_ylabel('y')

ax = plt.subplot2grid((2,2),(0, 1))
ax.plot(neuron[0, :], neuron[1, :], linewidth=2.)
for i in range(nb_perturbation_trials):
    ax.plot(final_output_neuron[i][t_offset + d_stim: t_offset + d_stim + d_neuron, 0, 0],
            final_output_neuron[i][t_offset + d_stim: t_offset + d_stim + d_neuron, 1, 0])
    ax.plot(final_output_neuron[i][subsampling_neuron, 0, 0],
            final_output_neuron[i][subsampling_neuron, 1, 0], 'bo')
ax.set_xlabel('x')
ax.set_ylabel('y')

ax = plt.subplot2grid((2,2),(1, 0))
ax.plot(chaos[0, :], chaos[1, :], linewidth=2.)
for i in range(nb_perturbation_trials):
    ax.plot(perturbation_output_chaos[i][t_offset + d_stim: t_offset + d_stim + d_chaos, 0, 0],
            perturbation_output_chaos[i][t_offset + d_stim: t_offset + d_stim + d_chaos, 1, 0])
    ax.plot(perturbation_output_chaos[i][subsampling_chaos, 0, 0],
            perturbation_output_chaos[i][subsampling_chaos, 1, 0], 'bo')
ax.set_xlabel('x')
ax.set_ylabel('y')

ax = plt.subplot2grid((2,2),(1, 1))
ax.plot(neuron[0, :], neuron[1, :], linewidth=2.)
for i in range(nb_perturbation_trials):
    ax.plot(perturbation_output_neuron[i][t_offset + d_stim: t_offset + d_stim + d_neuron, 0, 0],
        perturbation_output_neuron[i][t_offset + d_stim: t_offset + d_stim + d_neuron, 1, 0])
    ax.plot(perturbation_output_neuron[i][subsampling_neuron, 0, 0],
            perturbation_output_neuron[i][subsampling_neuron, 1, 0], 'bo')
ax.set_xlabel('x')
ax.set_ylabel('y')

plt.show()
