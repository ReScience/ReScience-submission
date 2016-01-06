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

import numpy as np
from copy import deepcopy
from scipy.stats import norm

# -----------------------------------
# FONCTIONS
# -----------------------------------
def createQValuesDict(states, actions):
    #return a big dict
    values = {0:np.zeros(len(states)*len(actions))}
    tmp = 0
    for s in states:
        values[s] = []
        n = 0
        for a in actions:
            values[(s,a)] = tmp
            values[s].append(tmp)
            values[(s,n)] = a
            tmp = tmp + 1
            n = n + 1
    return values

def transitionRules(state, action):
    if state == 's0' and action == 'pl':
        return 's1'
    else:
        return 's0'



def createRewardRateDict():
    return dict({'rate':0,
                 'reward':0,
                 'tau':0,
                 'step':0,
                 'tau_list':[0],
                 'r_list':[0],
                 ('s0', 'pl'):8,
                 ('s0', 'em'):12,
                 ('s1', 'pl'):12,
                 ('s1', 'em'):0,
                 'r*tau':[]})

def createCovarianceDict(n, init_cov, eta):
    #cov = np.eye(n)*init_cov+(1-np.eye(n,n))*((1-init_cov)/(n-1))    
    cov = np.eye(n)
    return dict({'cov':cov,
                 'noise':np.eye(n,n)*init_cov*eta})


def getBestAction(state, values, ind = 0):
    return values[(state, np.where(values[ind][values[state]] == np.max(values[ind][values[state]]))[0][0])]

def getBestActionSoftMax(state, values, beta, ind = 0):
    tmp = np.exp(values[ind][values[state]]*float(beta))
    tmp = tmp/float(np.sum(tmp))
    tmp = [np.sum(tmp[0:i]) for i in range(len(tmp))]
    return values[(state, np.sum(np.array(tmp) < np.random.rand())-1)]

def SoftMax(values, beta):
    tmp = np.exp(values*float(beta))
    tmp = tmp/float(np.sum(tmp))
    tmp = [np.sum(tmp[0:i]) for i in range(len(tmp))]
    return np.sum(np.array(tmp) < np.random.rand())-1

def computeEntropy(values, beta):
    tmp = np.exp(values*float(beta))
    tmp = tmp/float(np.sum(tmp))
    return -np.sum(tmp*np.log2(tmp))

def computeVPIValues(mean, variance):
    #WARNING input and output very specific
    # mean = array(current state), variance = array(current state)
    #vpi = array(current state)
    vpi = np.zeros((len(mean)))
    ind = np.argsort(mean)
    vpi[ind[-1]] = (mean[ind[-2]]-mean[ind[-1]])*norm.cdf(mean[ind[-2]], mean[ind[-1]], np.sqrt(variance[ind[-1]])) + (np.sqrt(variance[ind[-1]])/np.sqrt(2*np.pi))*np.exp(-(mean[ind[-2]]-mean[ind[-1]])**2/(2*variance[ind[-1]]))
    for i in range(len(mean)-2, -1, -1):
        vpi[ind[i]] = (mean[ind[i]]-mean[ind[-1]])*(1-norm.cdf(mean[ind[-1]], mean[ind[i]], np.sqrt(variance[ind[i]]))) + (np.sqrt(variance[ind[i]])/np.sqrt(2*np.pi))*np.exp(-(mean[ind[-1]]-mean[ind[i]])**2/(2*variance[ind[i]]))        
    return vpi
       
def updateRewardRate(reward_rate, sigma, delay = 0.0):
    return ((1-sigma)**(1+delay))*reward_rate['rate']+sigma*reward_rate['reward']

def updateQValuesHabitual(values, delta, alpha):
    return values+delta*alpha

def computeSigmaPoints(values, covariance, kappa=0.5):
    n = len(values)
    point = np.zeros((2*n+1,n))
    point[0] = values
    c = np.linalg.cholesky((n+kappa)*covariance)
    point[range(1,n+1)] = values+np.transpose(c)
    point[range(n+1, 2*n+1)] = values-np.transpose(c)
    weights = np.zeros((2*n+1,1))
    weights[0] = kappa/(n+kappa)
    weights[1:2*n+1] = 1/(2*n+kappa)
    return point, weights

def computeExplorationCost(reward_rate, tau, transition):
    return tau*reward_rate[transition]
    
def createTransitionDict(state1, action, state2, stopstate):
    transition = dict()
    n = float(len(set(state1)))
    transition[None] = stopstate
    for i,j,k in zip(state1, action, state2):
        transition[(i,j,k)] = 1/n
        transition[(i,j)] = k
    return transition

def updateTransition(transition_values, transition, phi):
    transition_values[transition] = (1-phi)*transition_values[transition] + phi
    for i in transition_values.iterkeys():
        if i != transition and i != None and len(i) == 3 and i[0] == transition[0]:
            transition_values[i] = (1-phi)*transition_values[i]
    return transition_values

def updateRewardsFunction(rewards_function, state, action, rau):
    rewards_function[1] = (1-rau)*rewards_function[1]+rau*rewards_function[0][rewards_function[(state, action)]]
    return rewards_function

def computeGoalValue(values, state, action, rewards, gamma, depth, phi, rau):
    rewards_function = deepcopy(rewards)
    rewards_function[1] = rewards_function[0].copy()
    rewards_function = updateRewardsFunction(rewards_function, state, action, rau)
    transition = createTransitionDict(['s0','s0','s1','s1'],['pl','em','pl','em'],['s1','s0','s0',None], 's0') #<====VERY BAD==============    NEXT_STATE = TRANSITION[(STATE, ACTION)]
    next_state = transition[(state, action)]
    if next_state == None:
        return rewards_function[1][rewards_function[(state, action)]] + gamma*transition[(state, action, next_state)]*np.max(values[0][values[transition[None]]])        
    else:
        transition = updateTransition(transition, (state, action, next_state), phi)
        tmp = np.max([computeGoalValueRecursive(values, next_state, a, rewards_function.copy(), transition.copy(), gamma, depth-1, phi, rau) for a in values[next_state]])
        value = rewards_function[1][rewards_function[(state, action)]] + gamma*transition[(state, action, next_state)]*tmp
        return value

def computeGoalValueRecursive(values, state, a, rewards_function, transition, gamma, depth, phi, rau):
    action = values[(state, values[state].index(a))]
    next_state = transition[(state, action)]
    rewards_function = updateRewardsFunction(rewards_function, state, action, rau)
    transition = updateTransition(transition, (state, action, next_state), phi)
    if next_state == None:
        return rewards_function[1][rewards_function[(state, action)]] + gamma*transition[(state, action, next_state)]*np.max(values[0][values[transition[None]]])        
    elif depth == 0:
        return rewards_function[1][rewards_function[(state, action)]] + gamma*transition[(state, action, next_state)]*np.max(values[0][values[next_state]])
    else:
        tmp = np.max([computeGoalValueRecursive(values, next_state, a, rewards_function.copy(), transition.copy(), gamma, depth-1, phi, rau) for a in values[next_state]])
        return rewards_function[1][rewards_function[(state, action)]] + gamma*transition[(state, action, next_state)]*tmp


def testQValues(states, values, beta, ind, niter):
    tmp = np.zeros((4)) #VBAAAAAAD
    for i in range(niter):
        for s in states:
            a = getBestActionSoftMax(s, values, beta, ind)
            tmp[values[(s,a)]] = tmp[values[(s, a)]] + 1 
    tmp = tmp/float(niter)
    return tmp
        
