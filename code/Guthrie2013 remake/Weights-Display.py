# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2015, Meropi Topalidou
# Distributed under the (new) BSD License.
#
# Contributors: Meropi Topalidou (Meropi.Topalidou@inria.fr)
# -----------------------------------------------------------------------------

# Display of the evolution of weights during trial in one
# simulation (mean out of 250 simulation).
# Run first 250_simulations-Guthrie.py or 250_simulations-Piron.py
# -----------------------------------------------------------------------------

import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


path = '../../data/Guthrie results/'
import sys
sys.path.append(path)
from parameters import *

suptitle = 'Protocol: Guthrie'
inverse = input('\nDo you want to have an inverse of Probabilities during the simulation or stop rewards?\nChoose: \n0 for No change\n1 for Reverse of Probabilities\n2 for Stopping Reward\n')
folder = path + 'Results'
if inverse == 1:
	inverse_trial = input('\nAfter how many trials will be the reverse?\n')
	inverse_all = input('\nDo you want to reverse all probabilities or just the middle ones?\nChoose:\n0 for middle ones\n1 for all\n')
	folder += '-inverse-after' + str(inverse_trial)
	folder += 'all' if inverse_all else 'middle cues'#-NoCortLearn-HalfParam
	title = 'Inverse probabilities of\n'
	title += 'all ' if inverse_all else 'middle cues '
	title += 'cues after %s trials' %str(inverse_trial)
elif inverse == 2:
	inverse_trial = input('\nAfter how many trials, the reward will stop?\n')
	inverse_all = 0
	folder += '-StopRewards-after' + str(inverse_trial)
	title = 'Stop rewards after %s' %str(inverse_trial)
else:
	folder += ''
	title = ''



file = folder + '/Weights_Str.npy'
load = np.load(file)
MeanWeights = (load.mean(axis = 0) - Wmin) / (Wmax - Wmin)
t = np.arange(MeanWeights[:,0].shape[0])
StdWeights = load.std(axis = 0)
file = folder+ '/MeanWeights-Str.npy'
np.save(file, MeanWeights)

fig, ax = plt.subplots(1)
ax.plot(MeanWeights[:,0], c = 'b', label = '1st cue')
ax.fill_between(t, MeanWeights[:,0]-StdWeights[:,0], MeanWeights[:,0]+StdWeights[:,0], facecolor='grey')
plt.plot(MeanWeights[:,1], c = 'r', label = '2nd cue')
ax.fill_between(t, MeanWeights[:,1]-StdWeights[:,1], MeanWeights[:,1]+StdWeights[:,1], facecolor='grey')
plt.plot(MeanWeights[:,2], c = 'g', label = '3rd cue')
ax.fill_between(t, MeanWeights[:,2]-StdWeights[:,2], MeanWeights[:,2]+StdWeights[:,2], facecolor='grey')
plt.plot(MeanWeights[:,3], c = 'm', label = '4th cue')
ax.fill_between(t, MeanWeights[:,3]-StdWeights[:,3], MeanWeights[:,3]+StdWeights[:,3], facecolor='grey')
plt.legend(loc='upper left')

plt.ylabel("Weights")
plt.xlabel("Trial number")
temp_title = 'Weights from\nCognitive Cortex to \nCognitive Striatum'
plt.title(temp_title, fontsize=12)
plt.title(suptitle, loc='left', fontsize=12)
plt.title(title, loc='right', fontsize=12)
#plt.ylim(0.4,0.80)
file = folder + "/Weights-Str.png"
fig.savefig(file)

plt.show()
