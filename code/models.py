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

import sys
import os
import numpy as np
from fonctions import *
from scipy.stats import norm

class KalmanQLearning():
    """ Class that implement a KalmanQLearning : 
    Kalman Temporal Differences : The deterministic case, Geist & al, 2009
    """

    def __init__(self, name, states, actions, gamma, beta, eta, var_obs, init_cov, kappa):
        self.name=name
        self.gamma=gamma;self.beta=beta;self.eta=eta;self.var_obs=var_obs;self.init_cov=init_cov;self.kappa=kappa
        self.values = createQValuesDict(states, actions)
        self.covariance = createCovarianceDict(len(states)*len(actions), self.init_cov, self.eta)
        self.states = states
        self.actions = actions
        self.state = None
        self.action = None
        self.reward = None

    def getAllParameters(self):        
        return dict({#'gamma':[0.1,self.gamma,0.9],
                     #'beta':[1,self.beta,5],
                     'eta':[1.0e-6,self.eta,0.001],
                     'var_obs':[0.01,self.var_obs,0.07]})
                     #'init_cov':[5,self.init_cov,15]})

    def setAllParameters(self, dict_p):
        for i in dict_p.iterkeys():
            self.setParameter(i,dict_p[i][1])

    def getParameter(self, name):
        if name == 'gamma':
            return self.gamma
        elif name == 'beta':
            return self.beta
        elif name == 'eta':
            return self.eta
        else:
            print "Parameters not found"
            sys.exit(0)     

    def setParameter(self, name, value):
        if name == 'gamma':
            self.gamma = value
        elif name == 'beta':
            self.beta = value
        elif name == 'eta':
            self.eta = value
        elif name == 'var_obs':
            self.var_obs = value
        elif name == 'init_cov':
            self.init_cov = value
        else:
            print "Parameters not found"
            sys.exit(0)

    def initialize(self):
        self.values = createQValuesDict(self.states, self.actions)
        self.covariance = createCovarianceDict(len(self.states)*len(self.actions), self.init_cov, self.eta)
        self.state = None
        self.action = None
        self.reward = None

    def chooseAction(self, state):
        self.state = state        
        self.covariance['noise'] = self.covariance['cov']*self.eta
        self.covariance['cov'][:,:] = self.covariance['cov'][:,:] + self.covariance['noise']
        self.action = getBestActionSoftMax(state, self.values, self.beta, 0)
        return action

    def updateValue(self, reward):
        r = (reward==1)*1        
        s = self.state
        a = self.action
        sigma_points, weights = computeSigmaPoints(self.values[0], self.covariance['cov'], self.kappa)
        rewards_predicted = (sigma_points[:,self.values[(s,a)]]-self.gamma*np.max(sigma_points[:,self.values[s]], 1)).reshape(len(sigma_points), 1)
        reward_predicted = np.dot(rewards_predicted.flatten(), weights.flatten())
        cov_values_rewards = np.sum(weights*(sigma_points-self.values[0])*(rewards_predicted-reward_predicted), 0)
        cov_rewards = np.sum(weights*(rewards_predicted-reward_predicted)**2) + self.var_obs
        kalman_gain = cov_values_rewards/cov_rewards
        self.values[0] = self.values[0] + kalman_gain*(reward-reward_predicted)
        self.covariance['cov'][:,:] = self.covariance['cov'][:,:] - (kalman_gain.reshape(len(kalman_gain), 1)*cov_rewards)*kalman_gain

    def predictionStep(self):
        self.covariance['noise'] = self.covariance['cov']*self.eta
        self.covariance['cov'][:,:] = self.covariance['cov'][:,:] + self.covariance['noise']

    def updatePartialValue(self, s, a, n_s, reward):
        #FOr keramati model, action selection is made outside the class
        self.state = s
        self.action = a        
        r = (reward==1)*1
        sigma_points, weights = computeSigmaPoints(self.values[0], self.covariance['cov'], self.kappa)
        rewards_predicted = (sigma_points[:,self.values[(s,a)]]-self.gamma*np.max(sigma_points[:,self.values[n_s]], 1)).reshape(len(sigma_points), 1)
        reward_predicted = np.dot(rewards_predicted.flatten(), weights.flatten())
        cov_values_rewards = np.sum(weights*(sigma_points-self.values[0])*(rewards_predicted-reward_predicted), 0)
        cov_rewards = np.sum(weights*(rewards_predicted-reward_predicted)**2) + self.var_obs
        kalman_gain = cov_values_rewards/cov_rewards
        self.values[0] = self.values[0] + kalman_gain*(reward-reward_predicted)
        self.covariance['cov'][:,:] = self.covariance['cov'][:,:] - (kalman_gain.reshape(len(kalman_gain), 1)*cov_rewards)*kalman_gain
        

class Keramati():
    """Class that implement Keramati models for action selection
    """
    
    def __init__(self, kalman,depth,phi, rau, sigma, tau):
        self.kalman = kalman
        self.depth = depth; self.phi = phi; self.rau = rau;self.sigma = sigma; self.tau = tau
        self.actions = kalman.actions; self.states = kalman.states
        self.values = createQValuesDict(kalman.states, kalman.actions)
        self.rfunction = createQValuesDict(kalman.states, kalman.actions)
        self.vpi = dict.fromkeys(self.states,list())
        self.rrate = 0.0
        self.state = None
        self.action = None
        self.transition = createTransitionDict(['s0','s0','s1','s1'],
                                               ['pl','em','pl','em'],
                                               ['s1','s0','s0','s0'], 's0') #<====VERY BAD==============    NEXT_STATE = TRANSITION[(STATE, ACTION)]        
    def initialize(self):
        self.values = createQValuesDict(self.states, self.actions)
        self.rfunction = createQValuesDict(self.states, self.actions)
        self.vpi = dict.fromkeys(self.states,list())
        self.rrate = 0.0
        self.state = None
        self.action = None
        self.transition = createTransitionDict(['s0','s0','s1','s1'],
                                               ['pl','em','pl','em'],
                                               ['s1','s0','s0','s0'], 's0')
                
    def chooseAction(self, state):
        self.state = state
        self.kalman.predictionStep()
        vpi = computeVPIValues(self.kalman.values[0][self.kalman.values[self.state]], self.kalman.covariance['cov'].diagonal()[self.kalman.values[self.state]])
        for i in range(len(vpi)):
            if vpi[i] >= self.rrate*self.tau:
                depth = self.depth
                self.values[0][self.values[(self.state, self.actions[i])]] = self.computeGoalValue(self.state, self.actions[i], depth)
            else:
                self.values[0][self.values[(self.state, self.actions[i])]] = self.kalman.values[0][self.kalman.values[(self.state,self.actions[i])]]

        self.action = getBestActionSoftMax(state, self.values, self.kalman.beta)
        return self.action

    def updateValues(self, reward, next_state):
        self.updateRewardRate(reward, delay = 0.0)
        self.kalman.updatePartialValue(self.state, self.action, next_state, reward)
        self.updateRewardFunction(self.state, self.action, reward)
        self.updateTransitionFunction(self.state, self.action)

    def updateRewardRate(self, reward, delay = 0.0):        
        self.rrate = ((1-self.sigma)**(1+delay))*self.rrate+self.sigma*reward

    def updateRewardFunction(self, state, action, reward):
        self.rfunction[0][self.rfunction[(state, action)]] = (1-self.rau)*self.rfunction[0][self.rfunction[(state, action)]]+self.rau*reward

    def updateTransitionFunction(self, state, action):
        #This is cheating since the transition is known inside the class
        #Plus assuming the transition are deterministic
        nextstate = self.transition[(state, action)]
        for i in [nextstate]:
            if i == nextstate:
                self.transition[(state, action, nextstate)] = (1-self.phi)*self.transition[(state, action, nextstate)]+self.phi
            else:
                self.transition[(state, action, i)] = (1-self.phi)*self.transition[(state, action, i)]
        
    def computeGoalValue(self, state, action, depth):
        next_state = self.transition[(state, action)]
        tmp = np.max([self.computeGoalValueRecursive(next_state, a, depth-1) for a in self.values[next_state]])
        value =  self.rfunction[0][self.rfunction[(state, action)]] + self.kalman.gamma*self.transition[(state, action, next_state)]*tmp
        return value

    def computeGoalValueRecursive(self, state, a, depth):
        action = self.values[(state, self.values[state].index(a))]
        next_state = self.transition[(state, action)]
        if depth:
            tmp = np.max([self.computeGoalValueRecursive(next_state, a, depth-1) for a in self.values[next_state]])
            return self.rfunction[0][self.rfunction[(state, action)]] + self.kalman.gamma*self.transition[(state, action, next_state)]*tmp
        else:
            return self.rfunction[0][self.rfunction[(state, action)]] + self.kalman.gamma*self.transition[(state, action, next_state)]*np.max(self.kalman.values[0][self.kalman.values[(state, action)]])        
        
