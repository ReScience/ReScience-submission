#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------------
# Copyright (c) 2015, Guillaume Viejo
# Distributed under the (new) BSD License.
#
# Contributors : Guillaume Viejo (guillaume.viejo@isir.upmc.fr)
# -------------------------------------------------------------------------
# References:
#
# Speed/accuracy trade-off between the habitual and the goal-directed processes.
# M. Keramati, A. Dezfouli, P. Piray. Plos Comp Bio, 7(5), 2011
# -------------------------------------------------------------------------
from __future__ import print_function
try:
	import xrange as range # Overrides range in Python2, does nothing in Python 3
except Exception:
	pass
import numpy as np
from fonctions import *
from models import *
from matplotlib import *
from pylab import *
import sys


# -------------------------------------------------------------------------
# PARAMETERS + INITIALIZATION
# -------------------------------------------------------------------------
eta 		= 0.0001    # variance of evolution noise v / Same as the original article
var_obs 	= 0.05   	# variance of observation noise n / Same as the original article
beta 		= 1.0       # rate of exploration / Same as the original article
gamma 		= 0.95     	# discount factor / Same as the original article
sigma 		= 0.02     	# updating rate of the average reward / Same as the original article
rau 		= 0.1       # update rate of the reward function / Same as the original article
tau 		= 0.08      # time step for graph exploration / Same as the original article
depth 		= 3        	# depth of search when computing the goal value / Same as the original article	
phi 		= 0.5       # update rate of the transition function / Not mentionned in the original article
init_cov 	= 1.0   	# initialisation of covariance matrice / Not mentionned in the original article
kappa 		= 0.1      	# unscentered transform parameters / Not mentionned in the original article

nb_iter_test 	= 500   # number of iterations in test
nb_iter_mod 	= 100	# number of iterations in moderate training
deval_mod_time 	= 40	# timing of devaluation in moderate training
nb_iter_ext 	= 300	# number of iterations in extensive training
deval_ext_time 	= 240	# timing of devaluation in extensive training
nb_blocs 		= 25 	# number of repetitions of the experiment for averages

states 		= ['s0', 's1']
actions 	= ['pl', 'em']
rewards 	= createQValuesDict(states, actions) 
kalman 		= KalmanQLearning("", states, actions, gamma, beta, eta, var_obs, init_cov, kappa) 	# Instance of kalman filter
selection 	= Keramati(kalman, depth, phi, rau, sigma, tau) 									# Instance of Keramati model of selection

# -------------------------------------------------------------------------
# Variables to plot
# -------------------------------------------------------------------------
data = {}
for i,l in zip(['mod', 'ext'],[nb_iter_mod,nb_iter_ext]):
	data[i] = {	'vpi'	:np.zeros((nb_blocs,l,2)),
				'r'		:np.zeros((nb_blocs,l)),
				'p'		:np.zeros((nb_blocs,l,4)),
				'q'		:np.zeros((nb_blocs,l))
				}				

# -------------------------------------------------------------------------
# Training + devaluation
# -------------------------------------------------------------------------
for i in xrange(nb_blocs):
	for exp, nb_trials, deval_time in zip(['mod','ext'], [nb_iter_mod, nb_iter_ext], [deval_mod_time, deval_ext_time]):
		kalman.initialize()
		selection.initialize()
		state = 's0'
		rewards[0][rewards[('s1','em')]] = 1.0
		print(exp, nb_trials, deval_time)
		for j in xrange(nb_trials):
			#Setting Reward
			if j == deval_time:
				rewards[0][rewards[('s1','em')]] = 0.0
				selection.rfunction[0][selection.rfunction[('s1', 'em')]] = -1.0
			#Learning
			while True:
				action = selection.chooseAction(state)        
				next_state = transitionRules(state, action)
				selection.updateValues(rewards[0][rewards[(state, action)]], next_state)				
				if state == 's1' and action == 'em':
					#Retrieving data					
					data[exp]['vpi'][i,j] = computeVPIValues(kalman.values[0][kalman.values['s0']],kalman.covariance['cov'].diagonal()[kalman.values['s0']])
					data[exp]['r'][i,j] = selection.rrate * selection.tau			
					data[exp]['p'][i,j] = testQValues(states, selection.values, kalman.beta, 0, nb_iter_test)
					data[exp]['q'][i,j] = kalman.values[0][kalman.values[('s0','pl')]]-kalman.values[0][kalman.values[('s0','em')]]
					state = next_state
					break
				else:
					state = next_state

meandata = {i:{j:np.mean(data[i][j], 0) for j in ['vpi', 'r', 'p']} for i in ['mod','ext']}

# -----------------------------------
# Plot
# -----------------------------------
dashes = ['--', '-.', '-']
colors = ['black', 'grey']


fig1 = figure(figsize = (15,9))

subplot(221)
for s in ['s0']:
	for a,i in zip(actions, range(len(actions))):
		plot(meandata['mod']['vpi'][:,selection.values[(s,a)]], linestyle = '-', color = colors[i], label = "VPI("+s+","+a+")", linewidth = 2)
plot(meandata['mod']['r'], color = 'black', label = "R*tau", linestyle = '--', linewidth = 2)
axvline(deval_mod_time, color='black', linewidth = 2)
legend()
grid()
ylim(0,0.1)
xlabel("Trial")

subplot(223)
for s in ['s0']:
	for a,i in zip(actions, range(len(actions))):
		plot(meandata['mod']['p'][:,selection.values[(s,a)]], linestyle = '-', color = colors[i], label = "P("+s+","+a+")", linewidth = 1.5)
axvline(deval_mod_time, color='black', linewidth = 2)
ylim(0.1,0.9)
xlabel("Trial")
ylabel("P(s,a)")
yticks(np.arange(0.2, 0.9, 0.1))
grid()
legend()

subplot(222)
for s in ['s0']:
	for a,i in zip(actions, range(len(actions))):
		plot(meandata['ext']['vpi'][:,selection.values[(s,a)]], linestyle = '-', color = colors[i], label = "VPI("+s+","+a+")", linewidth = 2)
plot(meandata['ext']['r'], linestyle = '--', color = 'black', label = "R*tau", linewidth = 2)
axvline(deval_ext_time, color='black', linewidth = 2)
legend()
xlabel("Trial")
grid()
ylim(0,0.1)


subplot(224)
for s in ['s0']:
	for a,i in zip(actions, range(len(actions))):
		plot(meandata['ext']['p'][:,selection.values[(s,a)]], linestyle = '-', color = colors[i], label = "P("+s+","+a+")", linewidth = 1.5)
axvline(deval_ext_time, color='black', linewidth = 2)
ylim(0.1,0.9)
grid()
xlabel("Trial")
ylabel("P(s,a)")
yticks(np.arange(0.2, 0.9, 0.1))
legend()

subplots_adjust(left = 0.08, wspace = 0.3, right = 0.86, hspace = 0.35)
figtext(0.12, 0.93, "Moderate pre-devaluation training", fontsize = 18)
figtext(0.55, 0.93, "Extensive pre-devaluation training", fontsize = 18)
figtext(0.06, 0.92, 'A', fontsize = 20)
figtext(0.06, 0.45, 'C', fontsize = 20)
figtext(0.49, 0.92, 'B', fontsize = 20)
figtext(0.49, 0.45, 'D', fontsize = 20)

fig1.savefig('fig.pdf', bbox_inches='tight')

show()
# -----------------------------------








