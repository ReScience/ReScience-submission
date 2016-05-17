from __future__ import print_function
import numpy as np
import time

from RecurrentNetwork import RecurrentNetwork

###############
# Parameters
###############
nb_learning_trials_rec = 20 # Number of learning trials for the recurrent weights
nb_learning_trials_readout = 10 # Number of learning trials for the readout weights
nb_perturbation_trials = 5 # Number of perturbation trials

stimulus_amplitude = 5.0 # Amplitude of the input pulse
t_offset = 200 # Time to wait before the stimulation
d_stim= 50 # Duration of the stimulation
d_trajectory = 2150 # Duration of the desired trajectory
t_relax = 550 # Duration to relax after the trajectory
trial_duration = t_offset + d_stim + d_trajectory + t_relax # Total duration of a trial

perturbation_amplitude = 0.5 # Amplitude of the perturbation pulse
t_perturbation = 500. # Offset for the perturbation
d_perturbation = 10 # Duration of the perturbation

target_baseline = 0.2 # Baseline of the target function
target_amplitude = 1. # Maximal value of the target function
target_width = 30. # Width of the Gaussian
target_time = d_trajectory - 150 # Peak time within the learning interval

####################
# Create the network
####################
net = RecurrentNetwork(
    Ni = 2, # Number of inputs
    N = 800, # Number of recurrent neurons
    No = 1, # Number of read-out neurons
    tau = 10.0, # Time constant of the neurons
    g = 1.8, # Synaptic strength scaling
    pc = 0.1, # Connection probability
    Io = 0.001, # Noise variance
    delta = 1.0, # Initial diagonal value of the P matrix
    P_plastic = 0.6, # Percentage of neurons receiving plastic synapses
)

###################
# Input definitions
###################
# Impulse after 200 ms
impulse = np.zeros((net.Ni, 1, trial_duration))
impulse[0, 0, t_offset:t_offset+d_stim] = stimulus_amplitude

# Perturbation during the trial
perturbation = np.zeros((net.Ni, 1, trial_duration))
perturbation[0, 0, t_offset : t_offset + d_stim] = stimulus_amplitude
perturbation[1, 0, t_offset + t_perturbation: t_offset + t_perturbation + d_perturbation] = perturbation_amplitude

# Target output for learning the readout weights
target = np.zeros((net.No, 1, trial_duration))
time_axis = np.linspace(0, trial_duration, trial_duration)
target[0, 0, : ] = target_baseline + (target_amplitude - target_baseline) * np.exp(-(t_offset + d_stim + target_time - time_axis)**2/target_width**2)

###################
# Main procedure
###################
tstart = time.time()

# Initial trial to determine the innate trajectory
print('Initial trial to determine a trajectory (without noise)')
trajectory, initial_output = net.simulate(stimulus=impulse, noise=False)

# Pre-training test trial
print('Pre-training test trial')
pretraining, pretraining_output = net.simulate(stimulus=impulse)

# Perturbation trial
print(nb_perturbation_trials, 'perturbation trials')
perturbation_initial = []
for i in range(nb_perturbation_trials):
    _, perturbation_output = net.simulate(stimulus=perturbation)
    perturbation_initial.append(perturbation_output)

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

# Save the network at the end of learning
net.save('network-simple.npz')

# Test trial
print('Test trial')
reproduction, final_output = net.simulate(stimulus=impulse)

# Perturbation trial
print(nb_perturbation_trials, 'perturbation trials')
perturbation_final = []
for i in range(nb_perturbation_trials):
    _, o = net.simulate(stimulus=perturbation)
    perturbation_final.append(o)

print('Simulation done in', time.time() - tstart, 'seconds.')

##################
# Visualization
##################
import matplotlib.pyplot as plt
import matplotlib.patches as patches

ax = plt.subplot2grid((4,2),(0, 0), colspan=2)
ax.imshow(trajectory[:, :100, 0].T, aspect='auto', origin='lower')
ax.set_xticks([0, 500, 1000, 1500, 2000, 2500])
ax.set_xticklabels([0, 0.5, 1, 1.5, 2, 2.5])
ymin, ymax = ax.get_ylim()
ax.add_patch(patches.Rectangle((t_offset, ymin), d_stim, ymax-ymin, alpha=0.1))
ax.set_title('Trajectory')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Recurrent units')

ax = plt.subplot2grid((4,2),(1, 0))
ax.plot(trajectory[:, 0, 0] + 1, 'b')
ax.plot(trajectory[:, 1, 0] + 3, 'b')
ax.plot(trajectory[:, 2, 0] + 5, 'b')
ax.plot(pretraining[:, 0, 0] + 1, 'r')
ax.plot(pretraining[:, 1, 0] + 3, 'r')
ax.plot(pretraining[:, 2, 0] + 5, 'r')
ax.set_xticks([0, 500, 1000, 1500, 2000, 2500])
ax.set_xticklabels([0, 0.5, 1, 1.5, 2, 2.5])
ax.set_yticks([1, 3, 5])
ax.set_yticklabels([0, 0, 0])
ymin, ymax = ax.get_ylim()
ax.add_patch(patches.Rectangle((t_offset, ymin), d_stim, ymax-ymin, alpha=0.7))
ax.add_patch(patches.Rectangle((t_offset+d_stim, ymin), d_trajectory, ymax-ymin, alpha=0.1))
ax.set_title('Pre-training')
ax.set_ylabel('Firing rate r')

ax = plt.subplot2grid((4,2),(1, 1))
ax.plot(trajectory[:, 0, 0] + 1, 'b')
ax.plot(trajectory[:, 1, 0] + 3, 'b')
ax.plot(trajectory[:, 2, 0] + 5, 'b')
ax.plot(reproduction[:, 0, 0] + 1, 'r')
ax.plot(reproduction[:, 1, 0] + 3, 'r')
ax.plot(reproduction[:, 2, 0] + 5, 'r')
ax.set_xticks([0, 500, 1000, 1500, 2000, 2500])
ax.set_xticklabels([0, 0.5, 1, 1.5, 2, 2.5])
ax.set_yticks([1, 3, 5])
ax.set_yticklabels([0, 0, 0])
ymin, ymax = ax.get_ylim()
ax.add_patch(patches.Rectangle((t_offset, ymin), d_stim, ymax-ymin, alpha=0.7))
ax.add_patch(patches.Rectangle((t_offset+d_stim, ymin), d_trajectory, ymax-ymin, alpha=0.1))
ax.set_title('Post-training')

ax = plt.subplot2grid((4,2),(2, 0))
ax.plot(initial_output[:, 0, 0], 'b')
ax.plot(pretraining_output[:, 0, 0], 'r')
ax.plot(target[0, 0, :], 'k')
ax.set_xticks([0, 500, 1000, 1500, 2000, 2500])
ax.set_xticklabels([0, 0.5, 1, 1.5, 2, 2.5])
ax.set_yticks([-2, -1, 0, 1, 2])
ax.set_ylim((-2, 2))
ymin, ymax = ax.get_ylim()
ax.add_patch(patches.Rectangle((t_offset, ymin), d_stim, ymax-ymin, alpha=0.7))
ax.add_patch(patches.Rectangle((t_offset+d_stim, ymin), d_trajectory, ymax-ymin, alpha=0.1))
ax.set_ylabel('Output (test)')

ax = plt.subplot2grid((4,2),(2, 1))
ax.plot(initial_output[:, 0, 0], 'b')
ax.plot(final_output[:, 0, 0], 'r')
ax.plot(target[0, 0, :], 'k')
ax.set_xticks([0, 500, 1000, 1500, 2000, 2500])
ax.set_xticklabels([0, 0.5, 1, 1.5, 2, 2.5])
ax.set_yticks([-2, -1, 0, 1, 2])
ax.set_ylim((-2, 2))
ymin, ymax = ax.get_ylim()
ax.add_patch(patches.Rectangle((t_offset, ymin), d_stim, ymax-ymin, alpha=0.7))
ax.add_patch(patches.Rectangle((t_offset+d_stim, ymin), d_trajectory, ymax-ymin, alpha=0.1))

ax = plt.subplot2grid((4,2),(3, 0))
for i in range(nb_perturbation_trials):
    ax.plot(perturbation_initial[i][:, 0, 0])
ax.plot(target[0, 0, :], 'k')
ax.set_xticks([0, 500, 1000, 1500, 2000, 2500])
ax.set_xticklabels([0, 0.5, 1, 1.5, 2, 2.5])
ax.set_yticks([-2, -1, 0, 1, 2])
ax.set_ylim((-2, 2))
ymin, ymax = ax.get_ylim()
ax.add_patch(patches.Rectangle((t_offset, ymin), d_stim, ymax-ymin, alpha=0.7))
ax.add_patch(patches.Rectangle((t_offset + t_perturbation, ymin), d_perturbation, ymax-ymin, alpha=0.7))
ax.add_patch(patches.Rectangle((t_offset+d_stim, ymin), d_trajectory, ymax-ymin, alpha=0.1))
ax.set_xlabel('Time (s)')
ax.set_ylabel('Output (perturbed)')

ax = plt.subplot2grid((4,2),(3, 1))
for i in range(nb_perturbation_trials):
    ax.plot(perturbation_final[i][:, 0, 0])
ax.plot(target[0, 0, :], 'k')
ax.set_xticks([0, 500, 1000, 1500, 2000, 2500])
ax.set_xticklabels([0, 0.5, 1, 1.5, 2, 2.5])
ax.set_yticks([-2, -1, 0, 1, 2])
ax.set_ylim((-2, 2))
ymin, ymax = ax.get_ylim()
ax.add_patch(patches.Rectangle((t_offset, ymin), d_stim, ymax-ymin, alpha=0.7))
ax.add_patch(patches.Rectangle((t_offset + t_perturbation, ymin), d_perturbation, ymax-ymin, alpha=0.7))
ax.add_patch(patches.Rectangle((t_offset+d_stim, ymin), d_trajectory, ymax-ymin, alpha=0.1))
ax.set_xlabel('Time (s)')

plt.show()
