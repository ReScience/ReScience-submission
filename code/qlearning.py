#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Alexandre Coninx
    ISIR - CNRS / Sorbonne Université
    27/04/2017
    Modified by Loic Cressot (02/2018)
""" 

import numpy as np
from numpy.linalg import norm
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression, Ridge


import matplotlib.pyplot as plt
import time


###################################################################################################################
### Fitted Q-Iteration ###
###################################################################################################################


########################################################
class Fitted_QIteration:
	""" Fitted Q-Iteration """
	def __init__(self,n_rbf=100, n_actions=25):
		self.n_rbf = n_rbf
		self.n_actions = n_actions
		self.fspacesize = n_rbf*n_actions
		self.beta = np.random.normal(scale=(1./(self.fspacesize**2)),size=self.fspacesize) # Glorot-style initialization
		self.beta2 = np.random.normal(scale=(1./(self.fspacesize**2)),size=self.fspacesize) # Glorot-style initialization
		self.mapping = None
		#Learning stuff
		self.state_rbf_out = None
		self.input_features = None
		# all_features if a 3D matrix of shape (n_states, n_actions, self.fspacesize)
		# containing features for all (s,a) indices. Used to speed up computations
		self.all_features = None


	### Training code ###

	def compute_istate_action_feature(self,i_state,action):
		""" RBF activations associated to one S-A pair """
		# Build feature S-A feature space
		out = np.zeros(self.fspacesize)
		out[action*self.n_rbf:(action+1)*self.n_rbf] = self.state_rbf_out[i_state]		
		return out


	def map_states(self,states,actions):
		""" Compute S-A mappings for all states """
		self.state_rbf_out = np.zeros((states.shape[0], self.n_rbf))
		self.input_features = np.zeros((states.shape[0], self.fspacesize))
		for i in range(states.shape[0]):
			self.state_rbf_out[i] = rbf_function(states[i],self.mapping[0],self.mapping[1])
			self.input_features[i] = self.compute_istate_action_feature(i,actions[i])
		self.input_features = self.input_features[:-1] # Remove last one to make fit work


	def fit_sk(self, states, actions, rewards, gamma=0.9, max_epochs=10, recompute_mapping=False, alpha_l2=1.0, conv_threshold=1.0):
		"""
		Computes K-means and fitted Q iteration
		Return learnt policy

		Params :
		conv_threshold indicates a beta convergence treshold (>=0) for stopping learning 

		Operations :
		- Computes self.state_rbf_out
		- Fits the parameter self.beta
		"""
		assert(0.0<gamma and gamma<1.0)
		assert(0.0<=alpha_l2 and alpha_l2<=1.0)
		assert(conv_threshold>=0.0)
		
		print("\nStarting fitted Q iterations :\n")

		if type(states)==type(tuple()): # weird, states seems to come in a tuple
			states=states[0]		

		flatstates = flatten_batch(states)
		print("Generating features mapping...")
		if recompute_mapping or not self.mapping:			
			self.mapping = gen_rbf_mapping(flatstates, self.n_rbf)

		# computes self.state_rbf_out
		self.map_states(flatstates, actions)

		# create regression
		if alpha_l2==0:
			lr = LinearRegression(fit_intercept=False)
		else:
			lr = Ridge(fit_intercept=False, alpha=alpha_l2)
		lr2 = Ridge(fit_intercept=False, alpha=alpha_l2)
		
		# iteration params
		ndata = len(actions)
		self.q_estimates = np.zeros(ndata-1)			
		divide = (ndata-1)/10 # for printing iterations
		# self.losses = list()
		# self.betas = list()
				
		print('Beginning fitting iterations with max_epochs :' + str(max_epochs))

		# repeat state_rbf_out n_actions times
		repeat_rbf_out = np.tile(self.state_rbf_out, self.n_actions)

		# align data for fast computation
		repeat_rbf_out = repeat_rbf_out[1:]
		rewards = rewards[:-1]


		for e in range(max_epochs):
			t1=time.time()			

			# keep value of beta for beta convergence computation :
			old_beta = self.beta
										
			# fast matricial computation of q_estimates							
			mult_beta = np.multiply( repeat_rbf_out, self.beta ) # element-wise multiplication of matrix with vector
			array_list = np.split( mult_beta , self.n_actions, 1) # split in n_actions sub matrices
			sums = np.column_stack([np.sum(e,1) for e in array_list]) # stack column vectors into matrix
			max_sums = np.max(sums,1)

			self.q_estimates = rewards + gamma*max_sums

			lr.fit( self.input_features, self.q_estimates )

			# compute score on all dataset
			score = lr.score( self.input_features, self.q_estimates)

			self.beta = lr.coef_

			# compute beta convergence
			beta_conv = norm(self.beta-old_beta)/norm(old_beta)

			t2=time.time()
			print("Epoch : %3d   |   Fitting score : %0.3f / 1.00   |   Beta convergence : %0.3f   |    time : %0.3f seconds" % (e, score, beta_conv, t2-t1) )
			#print("Fitting score 2 "+str(e)+" : "+str(score2) + " in " + str(t2-t1) + " s\n")
			if beta_conv <= conv_threshold:
				print('Beta conv threhsold ' + str(conv_threshold) + ' reached, stopping iterations')
				break

		print('\nFitting done.\n')
		# return learnt policy
		return self.compute_Qmax
			
	# ### Functions for use in policy ###
	def compute_Q(self,s,a):
		""" Q(s,a) """			
		return np.dot(activations, self.beta)
	
	def compute_Qmax(self,s):
		""" Qmax(s) """
		# compute f_rbf activations
		activations = f_rbf(s,self.mapping)
		# element-wise multiplication of self.beta and activations repeated n_actions times
		mult_beta = np.multiply( np.tile(activations, self.n_actions), self.beta )
		# split array into n_actions sub arrays, compute their sums and return argmax
		return np.argmax( [ np.sum(e) for e in np.split( mult_beta , self.n_actions) ] )


	def compute_error(self,states, actions, rewards, gamma=0.9):
		""" RMS error of the fit """
		flatstates = flatten_batch(states)
		ndata = len(actions)
		squerror = 0.
		for i in range(ndata-1):
			q_estimate = rewards[i] + gamma*self.compute_Qmax(flatstates[i+1])
			squerror += (q_estimate - self.compute_Q(flatstates[i],actions[i]))**2
		return np.sqrt(squerror/(ndata-1))

#### end class Fitted_QIteration ##################################################


def stateKMeansClusterization(states,n_clusters=100):
	""" Clusterize state space points """
	#TODO see why n_jobs=-2 does not work on OSX 
	#km = KMeans(n_clusters=100,n_jobs=-2)
	km = KMeans(n_clusters=n_clusters, n_jobs=1)
	km.fit(states)
	centers = km.cluster_centers_
	return centers

def compute_mean_intercenter_distance(centers):
	""" Compute mean distance between clusters centroids """
	n = centers.shape[0]
	squared_distance_matrix = np.ones((n,n))*np.inf
	for i in range(n):
		for j in range(i):
			squared_distance_matrix[i,j] = np.sum((centers[i] - centers[j])**2)
	mindists = np.sqrt(squared_distance_matrix.min(axis=0)[:-1])
	return np.mean(mindists)


def rbf_function(x,mu,sigma):
	""" RBF function """
	squared_dist = np.sum((x-mu)**2, axis=1)
	return np.exp(-squared_dist/(2*(sigma**2)))


def gen_rbf_mapping(states,n_rbf=100):
	""" Get the centers and sigma of a clusterization of the state space """
	centers = stateKMeansClusterization(states, n_clusters=n_rbf)
	sigma = compute_mean_intercenter_distance(centers)
	return centers, sigma

def f_rbf(s, mapping):
	""" RBF activations associated to one S-A pair """
	activations = rbf_function(s.flatten(),mapping[0],mapping[1]) # Get activation
	return activations / activations.sum() # Normalize to 1t

def flatten_batch(states):
	""" Flatten an array except the first dim """
	collapseddims = np.product(states.shape[1:])
	return states.reshape((states.shape[0],collapseddims))



