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
from model import *
from display import *
from trial import *
from parameters import *

def sigmoid(V,Vmin=Vmin,Vmax=Vmax,Vh=Vh,Vc=Vc):
    return  Vmin + (Vmax-Vmin)/(1.0+np.exp((Vh-V)/Vc))

def testing_trials(f = None, trials = n_trials, debugging = True, save = False, protocol = 'Piron', familiar = True, hist = False):
	P, D, mBc, ABC, NoMove, RT = [], [], [], [], [], []
	AP = np.zeros(n)
	for j in range(trials):
		if debugging:
			print 'Trial: ', j + 1
		if hist:
			histor, time = trial(f = f, trial_n = j, learn = False, protocol = protocol, hist = hist, familiar = familiar, debugging = debugging, P = P, D = D, AP = AP, mBc = mBc, ABC = ABC, NoMove = NoMove)
		else:
			time = trial(f = f, trial_n = j, learn = False, protocol = protocol, hist = hist, familiar = familiar, debugging = debugging, P = P, D = D, AP = AP, mBc = mBc, ABC = ABC, NoMove = NoMove)
		RT.append(time)
		if not len(P) == j+1:
			P.append(0)
	debug(f = f, P = P, D = D, AP = AP, mBc = mBc, ABC = ABC, NoMove = NoMove)
	print
	if save:
		return P, RT, np.array(D).mean(), AP, np.array(mBc).mean(), np.array(ABC).mean(), len(NoMove)/float(trials)
	if hist:
		return histor
