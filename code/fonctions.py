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
from scipy.stats import norm

# -----------------------------------
# FONCTIONS
# -----------------------------------
def createQValuesDict(states, actions):
    """
    Returns a dictionary containing the Q-values of all state-action pairs.
    The keys are either :
    *states return action index
    *(state, action index) return action name
    *(state, action) return action index
    * 0 return values
    """
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
    """
    Set before runs
    Return next state
    """
    if state == 's0' and action == 'pl':
        return 's1'
    else:
        return 's0'

def createCovarianceDict(n, init_cov, eta):
    """
    Initialise covariance matrix 
    """
    cov = np.eye(n)
    return dict({'cov':cov,
                 'noise':np.eye(n,n)*init_cov*eta})

def getBestActionSoftMax(state, values, beta, ind = 0):
    """
    Compute soft-max function and return action
    """
    tmp = np.exp(values[ind][values[state]]*float(beta))
    tmp = tmp/float(np.sum(tmp))
    tmp = [np.sum(tmp[0:i]) for i in range(len(tmp))]
    return values[(state, np.sum(np.array(tmp) < np.random.rand())-1)]

def computeVPIValues(mean, variance):
    """
    Return VPI for each action from Kalman Q-Learning
    Input and output are very specific
    mean = array(current state), variance = array(current state)
    vpi = array(current state)    
    """
    vpi = np.zeros((len(mean)))
    ind = np.argsort(mean)
    tmp1 = (mean[ind[-2]]-mean[ind[-1]])*norm.cdf(mean[ind[-2]], mean[ind[-1]], np.sqrt(variance[ind[-1]]))
    tmp2 = (np.sqrt(variance[ind[-1]])/np.sqrt(2*np.pi))*np.exp(-(mean[ind[-2]]-mean[ind[-1]])**2/(2*variance[ind[-1]]))
    vpi[ind[-1]] = tmp1 + tmp2
    for i in range(len(mean)-2, -1, -1):
        tmp1 = (mean[ind[i]]-mean[ind[-1]])*(1-norm.cdf(mean[ind[-1]], mean[ind[i]], np.sqrt(variance[ind[i]])))
        tmp2 = (np.sqrt(variance[ind[i]])/np.sqrt(2*np.pi))*np.exp(-(mean[ind[-1]]-mean[ind[i]])**2/(2*variance[ind[i]]))
        vpi[ind[i]] = tmp1 + tmp2
    return vpi
       
def updateRewardRate(reward_rate, sigma, delay = 0.0):
    return ((1-sigma)**(1+delay))*reward_rate['rate']+sigma*reward_rate['reward']

def computeSigmaPoints(values, covariance, kappa=0.5):
    """
    Determinitics sampling of sigma-points around the current Q-VAlues of Kalman Q-Learning
    return array of sigma-points and weight
    See Kalman Temporal Differences: The deterministic case, Geist et al, 2009 for explanation of the method
    """
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
    
def createTransitionDict(state1, action, state2, stopstate):
    """
    Transition dictionnary
    Key (current state, action, next state) return transition probability
    Key (current state, action) return next state
    """
    transition = dict()
    n = float(len(set(state1)))
    transition[None] = stopstate
    for i,j,k in zip(state1, action, state2):
        transition[(i,j,k)] = 1/n
        transition[(i,j)] = k
    return transition


def testQValues(states, values, beta, ind, niter):
    """
    Called at the end of a trial to evaluate the probability of actions
    """
    tmp = np.zeros((4)) #VBAAAAAAD
    for i in range(niter):
        for s in states:
            a = getBestActionSoftMax(s, values, beta, ind)
            tmp[values[(s,a)]] = tmp[values[(s, a)]] + 1 
    tmp = tmp/float(niter)
    return tmp
        
