#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Loic Cressot
    ISIR - CNRS / Sorbonne Université
    10/2018
""" 

"""
    jonschkowski_priors2015 : module implementing the method from Jonschkowski and Brock
    in their 2015 paper 'Learning State Representations with Robotic Priors'
"""

import os, sys

import numpy as np 

import keras
from keras.models import Model, Sequential
from keras.layers import Dense, Input, Reshape, Flatten, Lambda, Concatenate, Add, Subtract, Multiply, GaussianNoise
from keras.optimizers import Adam, SGD


from keras import backend as K
from keras.optimizers import Adam
from keras import regularizers


class Priors_model():

    def __init__(self, obs_shape, state_dim, seed=1, learning_rate=0.001,
                l1_reg=0.001, loss_weights=(1.0,5.0,1.0,5.0),
                noise_stddev=1e-4, tanh_units=None
                ):

        self.state_dim = state_dim
        self.loss_weights = [K.constant(w) for w in loss_weights] # list of constant scalar tensors

        self.mean_obs = None
        self.std_obs = None

        # seed random number generator
        np.random.RandomState(seed)

        ## INPUTS
        # define the four observation inputs : observations at t1 and t2 and their respective next observations
        obs_input = Input(obs_shape, name='obs_input', dtype='float32')        
        next_obs_input = Input(obs_shape, name='next_obs_input', dtype='float32') 
        # and the two same_action and same_reward inputs
        same_action_input = Input((None,), name='same_action_input', dtype='bool')
        same_a_diff_r_input = Input((None,), name='same_a_diff_r_input', dtype='bool')


        ## MODEL
        # build a linear model, with L1 regularization and a noise layer on top of it (noise is only added at training time)
        #   initiliaze with uniform drawn in [-0.05; 0.05] and don't use bias
        self.linear_model = Sequential()
        self.linear_model.add( Flatten(input_shape=obs_shape) )
        # if specified, add tanh hidden unit to learn from more complex observations
        if tanh_units :
            self.linear_model.add( Dense(tanh_units, activation='tanh', kernel_regularizer=regularizers.l1(l1_reg), use_bias=False,
                                         kernel_initializer = keras.initializers.RandomUniform(minval=-0.05, maxval=0.05, seed=None)) ) 
        self.linear_model.add( Dense(state_dim, activation='linear', kernel_regularizer=regularizers.l1(l1_reg), use_bias=False,
                                     kernel_initializer = keras.initializers.RandomUniform(minval=-0.05, maxval=0.05, seed=None)) )
        self.linear_model.add( GaussianNoise(noise_stddev) ) # noise layer, only active at training

        ## LOSS
        ## As this is unsupervised learning, we can write the loss function with a Lambda layer, directly inside the model
        
        # compute state output for each observation input
        states = self.linear_model(obs_input)
        next_states = self.linear_model(next_obs_input)

        # compute difference of states :
        delta_states = Subtract()([next_states, states])

        # the output loss layer
        temporalcoherence_loss = Lambda( self.temporalcoherence_loss, output_shape=(1,), name='temporalcoherence_loss')(delta_states)
        proportionality_loss = Lambda( self.proportionality_loss, output_shape=(1,), name='proportionality_loss')([delta_states, same_action_input])
        causality_loss = Lambda( self.causality_loss, output_shape=(1,), name='causality_loss')([states, same_a_diff_r_input])
        repeatability_loss = Lambda( self.repeatability_loss, output_shape=(1,), name='repeatability_loss')([states, delta_states, same_action_input])

        # create a Lamdba layer to perform scalar mutliplication
        scalar_multiplication = lambda scalar, name : Lambda(lambda x : x*scalar, output_shape=(1,) , name=name)

        # create a list of the four weigthed loss which will be optimized jointly
        weighted_losses = [ scalar_multiplication(self.loss_weights[0],'weighted_temporalcoherence')([temporalcoherence_loss]),
                            scalar_multiplication(self.loss_weights[1],'weighted_proportionality')([proportionality_loss]),
                            scalar_multiplication(self.loss_weights[2],'weighted_causality')([causality_loss]),
                            scalar_multiplication(self.loss_weights[3],'weighted_repeatability')([repeatability_loss ]) ]

        ## COMPILE MODEL
        # build a model with the six inputs and the list of outputs loss layers, and compile it
        #self.model = Model(inputs=[obs_input, same_action_input, same_a_diff_r_input], outputs=weighted_losses)
        self.model = Model(inputs=[obs_input, next_obs_input, same_action_input, same_a_diff_r_input], outputs=weighted_losses)

        #   Keras needs a final loss function with (y_true, y_pred) parameters
        #   As the method requires an expectation over the batch lossess, we return the mean of the loss in y_pred
        #   y_true is not used as labels here but only for regularizing the conditional expectations in proportionality, causality and repeatability losses
        def mean_loss(y_true, y_pred):
            return K.mean(y_pred*y_true, axis=0)#/K.sum(K.switch())
        
        self.model.compile(loss=mean_loss, optimizer=Adam(learning_rate))


    def temporalcoherence_loss(self, delta_states):
        # temporalcoherence_loss        
        return K.sum(delta_states**2, axis=1)


    def proportionality_loss(self, args):
        delta_states_row, same_action_input = args

        delta_states_sqrt_row = K.sqrt(K.epsilon() + K.sum(delta_states_row**2, axis=1))
        delta_states_sqrt_row = K.tf.reshape(delta_states_sqrt_row, [-1,1])
        delta_states_sqrt_col = K.tf.transpose(delta_states_sqrt_row)

        pairwise_diff_sq = K.switch(same_action_input, (delta_states_sqrt_row - delta_states_sqrt_col)**2, K.zeros_like(same_action_input,dtype='float32') )

        return K.sum(pairwise_diff_sq, axis=1)

    def causality_loss(self, args):
        states_row, same_a_diff_r_input = args
        states_row1, states_row2 = K.tf.split(states_row, 2, axis=1)
        states_col1 = K.tf.transpose(states_row1)
        states_col2 = K.tf.transpose(states_row2)

        pairwise_exp = K.switch( same_a_diff_r_input, K.exp( -((states_row1-states_col1)**2 + (states_row2-states_col2)**2) ), K.zeros_like(same_a_diff_r_input,dtype='float32') )
        
        return K.sum(pairwise_exp,axis=1)

    def repeatability_loss(self, args):
        states_row, delta_states_row, same_action_input = args
        
        states_row1, states_row2 = K.tf.split(states_row, 2, axis=1)
        states_col1 = K.tf.transpose(states_row1)
        states_col2 = K.tf.transpose(states_row2)

        delta_states_row1, delta_states_row2 = K.tf.split(delta_states_row, 2, axis=1)
        delta_states_col1 = K.tf.transpose(delta_states_row1)
        delta_states_col2 = K.tf.transpose(delta_states_row2)

        pairwise_delta_diff_sq = K.switch(same_action_input, (delta_states_row1 - delta_states_col1)**2 + (delta_states_row2 - delta_states_col2)**2, K.zeros_like(same_action_input,dtype='float32'))
        pairwise_diff_sq = K.switch(same_action_input, K.exp(-((states_row1 - states_col1)**2 + (states_row2 - states_col2)**2)), K.zeros_like(same_action_input,dtype='float32') )

        return pairwise_delta_diff_sq*pairwise_diff_sq



    def training_generator(observations, actions, rewards, episode_starts, batch_size):
        """
        data generator for training the model batch by batch, each batch being indepently sampled
            this means that with this generator, the same data can be sampled several times within the same epoch
        yields : (X, Y) tuple for training, with null Y yielded only for keras to work
        """
        n_obs = observations.shape[0]

        if batch_size > n_obs-1 or batch_size < 1:
            raise ValueError('Wrong batch size')

        # if actions are one dimensional, add another dimension to avoid a bug later when performing product over axis 1
        if len(actions.shape)==1:
            actions = actions[:,np.newaxis];
        elif len(actions.shape)>2:
            raise ValueError('actions should be either a one are a two dimensional array')

        # indices for all time steps where the episode continues (to be able to count ∆s)
        all_indices = np.array([i for i in range(n_obs-1) if not episode_starts[i + 1]], dtype='int32')

        # finds all indices in a batch which actions are the same than a given action and with indices > index 
        find_same_actions = lambda index, batch: np.prod(actions[batch] == actions[batch[index]], axis=1) *\
                                                 (np.arange(batch.shape[0])>index)

        # finds all indices in a batch which actions are the same than a given action but rewards+1 are different and with indices > index
        find_same_actions_diff_rewards = lambda index, batch: np.prod(actions[batch] == actions[batch[index]], axis=1) *\
                                                                     (rewards[batch + 1] != rewards[batch[index] + 1]) *\
                                                                     (np.arange(batch.shape[0])>index)

        # take a 2D np array and add the line index to every value (gives 3D array with (line index, value) pairs )
        #add_line_index = lambda array2D : np.array([[ [i,array2D[i,j]] for j in range(array2D.shape[1])  ] for i in range(array2D.shape[0]) ]  )

        while True:
            # randomly sample a batch of indices :
            rand_indices = np.random.choice(all_indices, batch_size, replace=False)

            # build half matrix composed on lines of indices with the same action (respectively same action and different reward) than the action of the matrix line index
            same_actions = np.array([ find_same_actions(i, rand_indices) for i in range(batch_size)], dtype='bool')
            same_actions_diff_rewards = np.array([ find_same_actions_diff_rewards(i, rand_indices) for i in range(batch_size)], dtype='bool')
            
            ### TODO : matrix of smae_actions and diif

            X = [ observations[rand_indices], observations[rand_indices+1], same_actions, same_actions_diff_rewards ]

            # Y are not labels here because this is unsupervised learning.
            # We use Y to store a constant ratio corresponding to the samples ratio for computing conditional expectations
            # Indeed in keras we have to pass vectors of size batchsize, but when there are zeros we have to normalize the mean to have a right conditional mean
            # This is usefull for proportionality, causality and repeatability losses which have conditional expectations
            ratio_same_action = np.mean(same_actions.flatten()>0)*batch_size
            ratio_same_action_diff_reward = np.mean(same_actions_diff_rewards.flatten()>0)*batch_size
            Y = [np.ones(batch_size),
                 np.ones(batch_size)*ratio_same_action,
                 np.ones(batch_size)*ratio_same_action_diff_reward,
                 np.ones(batch_size)*ratio_same_action ]
            #Y = np.ones(batch_size)
            yield X,Y


    def learn(self, observations, actions, rewards, episode_starts, batch_size= 256, num_epochs=100, verbose=True):
        """
        Learn model on data
        """
        verbose = 2 if verbose else 0

        # center and scale observations
        self.mean_obs = np.mean(observations, axis=0, keepdims=True)
        self.std_obs = np.std(observations, ddof=1)
        observations = (observations - self.mean_obs) / self.std_obs
        
        num_samples = observations.shape[0]
        steps_per_epoch = int(num_samples/batch_size)

        train_gen = Priors_model.training_generator(observations, actions, rewards, episode_starts, batch_size)

        # fit model using batches generated by the training generator. No need to shuffle since it is already done by the generator
        # Note : it seems to be way faster with use_multiprocessing=False
        # returns the training History object
        return self.model.fit_generator(generator=train_gen, steps_per_epoch=steps_per_epoch, epochs=num_epochs,
                                 use_multiprocessing=False, shuffle=False, verbose=verbose )

        ## DEBUG


    def phi(self, observations, batch_size=256):
        """
        Compute transformation for one or more observations
        """    
        # center and scale observation    
        observations = (observations - self.mean_obs) / self.std_obs
        # compute and return associated state(s)
        return self.linear_model.predict(observations, batch_size=batch_size) 


