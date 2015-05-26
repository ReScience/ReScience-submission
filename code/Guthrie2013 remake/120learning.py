# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2015, Meropi Topalidou
# Distributed under the (new) BSD License.
#
# Contributors: Meropi Topalidou (Meropi.Topalidou@inria.fr)
#				Nicolas Rougier  (Nicolas.Rougier@inria.fr)
# -----------------------------------------------------------------------------

# Testing learning for each model under Guthrie protocol
# -----------------------------------------------------------------------------
if __name__ == "__main__":
	from model import *
	from display import *
	from learning import *
	reverse = input('\nDo you want to have an reverse of Probabilities during the simulation or stop rewards?\nChoose: \n0 for No change\n1 for Reverse of Probabilities\n2 for Stopping Reward\n')
	if reverse == 1:
		reverse_trial = input('\nAfter how many trials will be the reverse?\n')
		reverse_all = input('\nDo you want to reverse all probabilities or just the middle ones?\nChoose:\n0 for middle ones\n1 for all\n')
	elif reverse == 2:
		reverse_trial = input('\nAfter how many trials, the reward will stop?\n')
		reverse_all = 0
	else:
		reverse_trial = 0
		reverse_all = 0

	reset(protocol = 'Guthrie')
	hist, P = learning_trials(inversable = reverse, reverse_all = reverse_all, reverse_trial = reverse_trial, hist = True, protocol = 'Guthrie')
	if 0: display_all(hist, 3.0)#, "single-trial-all.pdf")
	if 0: display_ctx(hist, 3.0)
	if 0: display_all(hist, 3.0)#, "single-trial-all.pdf")
