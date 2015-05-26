# -----------------------------------------------------------------------------
# Copyright (c) 2015, Meropi Topalidou
# Distributed under the (new) BSD License.
#
# Contributors: Meropi Topalidou (Meropi.Topalidou@inria.fr)
#				Nicolas Rougier  (Nicolas.Rougier@inria.fr)
# -----------------------------------------------------------------------------

# Simulate number of experiments that is given in parameters.py of the different
# models. Each simulation is a number of trials under Guthrie protocol.
# -----------------------------------------------------------------------------

if __name__ == "__main__":
	import numpy as np
	import os
	from model import *
	from learning import *
	from testing import *
	from parameters import *
	reverse = input('\nDo you want to have an reverse of Probabilities during the simulation or stop rewards?\nChoose: \n0 for No change\n1 for Reverse of Probabilities\n2 for Stopping Reward\n')
	path = '../../data/Guthrie results/Results'
	if reverse == 1:
		reverse_trial = input('\nAfter how many trials will be the reverse?\n')
		reverse_all = input('\nDo you want to reverse all probabilities or just the middle ones?\nChoose:\n0 for middle ones\n1 for all\n')
		path += '-reverse-after' + str(reverse_trial)
		path += 'all' if reverse_all else 'middle cues'#-NoCortLearn-HalfParam
	elif reverse == 2:
		reverse_trial = input('\nAfter how many trials, the reward will stop?\n')
		reverse_all = 0
		path += '-StopRewards-after' + str(reverse_trial)
	else:
		reverse_trial = 0
		reverse_all = 0
	#path += '-9000'
	if not os.path.exists(path):
		os.makedirs(path)
	debugging = path + '/Debugging.txt'#.txt'
	#f = open(debugging, 'wb')


	CVtotal = np.zeros((simulations, n))
	WtotalSTR = np.zeros((simulations, n_trials, n))

	P = np.zeros((simulations, n_trials))
	RT = np.zeros((simulations, n_trials))
	D = np.zeros(simulations)
	RP = np.zeros((simulations, n))
	AP = np.zeros((simulations, n))
	mBc = np.zeros(simulations)
	ABC = np.zeros(simulations)
	NoMove = np.zeros(simulations)


	for i in range(simulations):
		print 'Experiment: ', i + 1
		reset(protocol = 'Guthrie')

		P[i,:], RT[i,:], D[i], RP[i,:], AP[i,:], mBc[i], ABC[i], NoMove[i], WtotalSTR[i] = learning_trials(inversable = reverse, reverse_all = reverse_all, reverse_trial = reverse_trial, protocol = 'Guthrie', trials = n_trials, debugging = False, save = True)

		CVtotal[i, :] = CUE["value"]
		print
		print


	debug_total(P, D, ABC, NoMove, AP, RP, CVtotal, WtotalSTR)
	file = path + '/Weights_Str.npy'
	np.save(file,WtotalSTR)

	#f.close()

	file = path + '/MeanCuesValues.npy'
	np.save(file,CVtotal)
	file = path + '/NoMove.npy'
	np.save(file,NoMove)
	file = path + '/RT.npy'
	np.save(file,RT)
	file = path + '/Performance.npy'
	np.save(file,P)
	file = path + '/DifferentChoices.npy'
	np.save(file,D)
	file = path + '/MotBefCog.npy'
	np.save(file,mBc)
	file = path + '/ActBefCues.npy'
	np.save(file,ABC)
