# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2015, Meropi Topalidou
# Distributed under the (new) BSD License.
#
# Contributors: Meropi Topalidou (Meropi.Topalidou@inria.fr)
#				Nicolas Rougier  (Nicolas.Rougier@inria.fr)
# -----------------------------------------------------------------------------
# References:
#
# * Interaction between cognitive and motor cortico-basal ganglia loops during
#   decision making: a computational study. M. Guthrie, A. Leblois, A. Garenne,
#   and T. Boraud. Journal of Neurophysiology, 109:3025–3040, 2013.
# -----------------------------------------------------------------------------
import numpy as np
import random
from trial import *
from parameters import *

def learning(reverse = False, reverse_all = True, f = None, trial_n = 0, debugging = True, protocol = 'Piron', familiar = True, learn = True, hist = False, P = [], D = [], mBc = [], ABC = [], NoMove = [], RT = [], RP = None, AP = None):
	if hist:
		histor, time = trial(reverse = reverse, reverse_all = reverse_all, f = f, trial_n = trial_n, learn = learn, protocol = protocol, hist = hist, familiar = familiar, debugging = debugging, P = P, D = D, RP = RP, AP = AP, mBc = mBc, ABC = ABC, NoMove = NoMove)
	else:
		time = trial(reverse = reverse, reverse_all = reverse_all, f = f, trial_n = trial_n, learn = learn, protocol = protocol, hist = hist, familiar = familiar, debugging = debugging, P = P, D = D, RP = RP, AP = AP, mBc = mBc, ABC = ABC, NoMove = NoMove)
	RT.append(time)
	if not len(P) == trial_n+1:
		P.append(0)
	if hist:
		return histor, RT, P, D, RP, AP, mBc, ABC, NoMove
	else:
		return RT, P, D, RP, AP, mBc, ABC, NoMove

def learning_trials(inversable = 0, reverse_all = True, reverse_trial = 50, less_trained_trials = 20, f = None, hist = False, trials = n_trials, debugging = True, save = False, protocol = 'Piron', familiar = True, type = 'learning', W_COG = None, W_MOT = None, W_STR = None, trained = False, j = -1):

	P, D, mBc, ABC, NoMove, RT = [], [], [], [], [], []
	RP = np.zeros(n)
	AP = np.zeros(n)
	wStr = np.zeros((trials,n))
	while not trained:
		P, D, mBc, ABC, NoMove, RT = [], [], [], [], [], []
		RP = np.zeros(n)
		AP = np.zeros(n)
		wStr = np.zeros((trials,n))
		reset(protocol = protocol, W_STR = W_STR)
		reverse = False
		for j in range(less_trained_trials):

			if debugging:
				print 'Trial: ', j + 1
			if hist:
				histor, RT, P, D, RP, AP, mBc, ABC, NoMove = learning(f = f, trial_n = j, debugging = debugging, protocol = protocol, familiar = familiar, hist = hist, P = P, D = D, mBc = mBc, ABC = ABC, NoMove = NoMove, RT = RT, RP = RP, AP = AP)
			else:
				RT, P, D, RP, AP, mBc, ABC, NoMove = learning(f = f, trial_n = j, debugging = debugging, protocol = protocol, familiar = familiar, hist = hist, P = P, D = D, mBc = mBc, ABC = ABC, NoMove = NoMove, RT = RT, RP = RP, AP = AP)

			wStr[j,:] = connections["CTX.cog -> STR.cog"].weights

		if type == 'learning':
			if np.mean(P) > 0.70:
				trained = True
		else:
			trained = True
	else:
		reverse = False
		for i in range(j+1, trials):

			if inversable == 1:
				if i == reverse_trial:
					if reverse_all:
						if protocol == 'Piron':
							CUE["reward"] = rewards_Piron_reverse
						elif protocol == 'Guthrie':
							CUE["reward"] = rewards_Guthrie_reverse_all
					else:
						CUE["reward"] = rewards_Guthrie_reverse_middle
					reverse = True
			elif inversable == 2:
				if i == reverse_trial:
					CUE["reward"] = 0., 0., 0., 0.
				reverse = False
			else:
				reverse = False
			if debugging:
				print 'Trial: ', i + 1
			if hist:
				histor, RT, P, D, RP, AP, mBc, ABC, NoMove = learning(reverse = reverse, reverse_all = reverse_all, f = f, trial_n = i, debugging = debugging, protocol = protocol, familiar = familiar, hist = hist, P = P, D = D, mBc = mBc, ABC = ABC, NoMove = NoMove, RT = RT, RP = RP, AP = AP)
			else:
				RT, P, D, RP, AP, mBc, ABC, NoMove = learning(reverse = reverse, reverse_all = reverse_all, f = f, trial_n = i, debugging = debugging, protocol = protocol, familiar = familiar, hist = hist, P = P, D = D, mBc = mBc, ABC = ABC, NoMove = NoMove, RT = RT, RP = RP, AP = AP)
			wStr[i,:] = connections["CTX.cog -> STR.cog"].weights
	debug(f = f, RT = RT, P = P, D = D, RP = RP, AP = AP, mBc = mBc, ABC = ABC, NoMove = NoMove)
	debug_learning(wStr[-1,:], cues_value = CUE["value"], f = f)
	if save:
		return P, RT, np.array(D).mean(), RP, AP, np.array(mBc).mean(), np.array(ABC).mean(), len(NoMove)/float(n_trials),wStr
	if hist:
		return histor, P
	else:
		return P
