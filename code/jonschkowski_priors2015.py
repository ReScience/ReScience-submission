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
        obs_input_t1 = Input(obs_shape, name='obs_input_t1')
        obs_input_t1_next = Input(obs_shape, name='obs_input_t1_next')
        obs_input_t2 = Input(obs_shape, name='obs_input_t2')
        obs_input_t2_next = Input(obs_shape, name='obs_input_t2_next')
        # and the two same_action and same_reward inputs
        same_action_input = Input((1,), name='same_action_input', dtype='float32')
        same_a_diff_r_input = Input((1,), name='same_a_diff_r_input', dtype='float32')

        all_inputs = [obs_input_t1, obs_input_t1_next, obs_input_t2, obs_input_t2_next, same_action_input, same_a_diff_r_input]

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
        s_t1= self.linear_model(obs_input_t1)
        s_t1_next= self.linear_model(obs_input_t1_next)
        s_t2= self.linear_model(obs_input_t2)
        s_t2_next = self.linear_model(obs_input_t2_next)        

        delta_st1 = Subtract()([ s_t1_next, s_t1 ])
        delta_st2 = Subtract()([ s_t2_next, s_t2 ])

        self.norm2 = lambda x : K.sum(x**2, axis=1)

        # don't forget to squeeze these vectors to make them 1-dimensional
        #same_action_input = K.squeeze(same_action_input, axis=1)
        #same_a_diff_r_input = K.squeeze(same_a_diff_r_input, axis=1)

        # the output loss layer
        #loss_out = Lambda(self.composite_loss, output_shape=(4,), name='loss_out')([s_t1, s_t1_next, s_t2, s_t2_next, same_action_input, same_a_diff_r_input ])
        temporalcoherence_loss = Lambda( self.temporalcoherence_loss, output_shape=(1,), name='temporalcoherence_loss')([delta_st1, delta_st2])
        proportionality_loss = Lambda( self.proportionality_loss, output_shape=(1,), name='proportionality_loss')([delta_st1, delta_st2, same_action_input])
        causality_loss = Lambda( self.causality_loss, output_shape=(1,), name='causality_loss')([s_t1, s_t2, same_a_diff_r_input])
        repeatability_loss = Lambda( self.repeatability_loss, output_shape=(1,), name='repeatability_loss')([s_t1, s_t2, delta_st1, delta_st2, same_action_input])

        # scalar_multiplication = lambda scalar : Lambda(lambda x : x*scalar, output_shape=(1,) )

        # weighted_losses = [ scalar_multiplication(self.loss_weights[0])([temporalcoherence_loss]),
        #                     scalar_multiplication(self.loss_weights[1])([proportionality_loss]),
        #                     scalar_multiplication(self.loss_weights[2])([causality_loss]),
        #                     scalar_multiplication(self.loss_weights[3])([repeatability_loss ]) ]

        weighted_losses = [ temporalcoherence_loss,
                            proportionality_loss,
                            causality_loss,
                            repeatability_loss ]

        print(weighted_losses)

        ## COMPILE MODEL
        # build a model with the six inputs and the output loss layer, and compile it
        #self.model = Model(inputs=all_inputs, outputs=loss_out)
        self.model = Model(inputs=all_inputs, outputs=weighted_losses)
        #self.model = Model(inputs=all_inputs, outputs=proportionality_loss)

        #   Keras needs a final loss function with (y_true, y_pred) parameters
        #   As the method requires an expectation over the batch lossess, we return the mean of the loss in y_pred
        #   y_true is unused here as this is unsupervised learning
        def mean_loss(y_true, y_pred):
            return K.mean(y_pred, axis=0)
        
        self.model.compile(loss=mean_loss, optimizer=Adam(learning_rate))



    def temporalcoherence_loss(self, args):
        delta_st1, delta_st2 = args
        # temporalcoherence_loss
        #return self.loss_weights[0]* (self.norm2(delta_st1) + self.norm2(delta_st2))/2.0
        return self.loss_weights[0]* (self.norm2(delta_st1) + self.norm2(delta_st2))/2.0


    def proportionality_loss(self, args):
        delta_st1, delta_st2, same_action_input = args
        # don't forget to squeeze these vectors to make them 1-dimensional (need to make that inside loss function and not outside for keras history issues)
        same_action = K.squeeze(same_action_input, axis=1)
        # proportionality_loss
        #   Note the K.mean(same_action, axis=0) normalization, because this is a conditional expectation on same_action
        #return self.loss_weights[1]*  Multiply()( [ ( K.sqrt(K.epsilon()+self.norm2(delta_st1)) - K.sqrt(K.epsilon()+self.norm2(delta_st2)) )**2,
        return self.loss_weights[1]*  K.switch( same_action,
                                                ( K.sqrt(K.epsilon()+self.norm2(delta_st1)) - K.sqrt(K.epsilon()+self.norm2(delta_st2)) )**2,
                                                K.zeros_like(same_action)
                                               )/K.mean(same_action, axis=0)
        #return self.loss_weights[1]*( K.sqrt(K.epsilon()+self.norm2(delta_st1)) - K.sqrt(K.epsilon()+self.norm2(delta_st2)) )**2
    
    def causality_loss(self, args):
        s_t1, s_t2, same_a_diff_r_input = args
        # don't forget to squeeze these vectors to make them 1-dimensional (need to make that inside loss function and not outside for keras history issues)
        same_a_diff_r = K.squeeze(same_a_diff_r_input, axis=1)
        # causality_loss
        #   Note the K.mean(same_a_diff_r, axis=0) normalization, because this is a conditional expectation on same_a_diff_r
        #return self.loss_weights[2]* Multiply()( [ K.exp( -self.norm2(s_t2-s_t1) ), same_a_diff_r ])/K.mean(same_a_diff_r, axis=0)
        return self.loss_weights[2]* K.switch( same_a_diff_r, K.exp( -self.norm2(s_t2-s_t1) ), K.zeros_like(same_a_diff_r) )/K.mean(same_a_diff_r, axis=0)

    def repeatability_loss(self, args):
        s_t1, s_t2, delta_st1, delta_st2, same_action_input = args
        # don't forget to squeeze these vectors to make them 1-dimensional (need to make that inside loss function and not outside for keras history issues)
        same_action = K.squeeze(same_action_input, axis=1)
        # repeatability_loss
        #   Note the K.mean(same_action, axis=0) normalization, because this is a conditional expectation on same_action
        #return self.loss_weights[3]* Multiply()([ K.exp(-self.norm2( s_t2-s_t1)), self.norm2(delta_st1-delta_st2), same_action ])/K.mean(same_action, axis=0)
        return self.loss_weights[3]* K.switch( same_action, Multiply()([K.exp(-self.norm2( s_t2-s_t1)), self.norm2(delta_st1-delta_st2)]), K.zeros_like(same_action)
                                             )/K.mean(same_action, axis=0)

   

    # def composite_loss(self, args):

    #     s_t1, s_t1_next, s_t2, s_t2_next, same_action_input, same_a_diff_r_input  = args

    #     delta_st1 = Subtract()([ s_t1_next, s_t1 ])
    #     delta_st2 = Subtract()([ s_t2_next, s_t2 ])

    #     norm2 = lambda x : K.sum(x**2, axis=1)

    #     # don't forget to squeeze these vectors to make them 1-dimensional
    #     same_action_input = K.squeeze(same_action_input, axis=1)
    #     same_a_diff_r_input = K.squeeze(same_a_diff_r_input, axis=1)
        
    #     # temporalcoherence_loss
    #     temporalcoherence_loss = (norm2(delta_st1) + norm2(delta_st2))/2.0

    #     # proportionality_loss
    #     #   Note the K.mean(same_action_input, axis=0) normalization, because this is a conditional expectation on same_action_input
    #     proportionality_loss = Multiply()( [ ( K.sqrt(K.epsilon()+norm2(delta_st1)) - K.sqrt(K.epsilon()+norm2(delta_st2)) )**2,
    #                                            same_action_input ])/K.mean(same_action_input, axis=0)

    #     # causality_loss
    #     #   Note the K.mean(same_a_diff_r_input, axis=0) normalization, because this is a conditional expectation on same_a_diff_r_input
    #     causality_loss = Multiply()( [ K.exp( -norm2(s_t2-s_t1) ), same_a_diff_r_input ])/K.mean(same_a_diff_r_input, axis=0)

    #     # repeatability_loss
    #     #   Note the K.mean(same_action_input, axis=0) normalization, because this is a conditional expectation on same_action_input
    #     repeatability_loss = Multiply()([ K.exp(-norm2( s_t2-s_t1)),
    #                                       norm2(delta_st1-delta_st2),
    #                                       same_action_input
    #                                     ])/K.mean(same_action_input, axis=0)

    #     # return a vector (needed by keras) which will then be sumed up 
    #     composite_loss =  self.loss_weights[0]*temporalcoherence_loss +\
    #                       self.loss_weights[1]*proportionality_loss +\
    #                       self.loss_weights[2]*causality_loss +\
    #                       self.loss_weights[3]*repeatability_loss

    #     return composite_loss
        

    def training_generator(observations, actions, rewards, episode_starts, batch_size):
        """
        data generator for training the model batch by batch, each batch being indepently sampled
            this means that with this generator, the same data can be sampled several times within the same epoch
        yields : (X, Y) tuple for training, with null Y yielded only for keras to work
        """
        n_obs = observations.shape[0] - 1

        if batch_size > n_obs-1 or batch_size < 1:
            raise ValueError('Wrong batch size')

        # indices for all time steps where the episode continues (to be able to count ∆s)
        all_indices = np.array([i for i in range(n_obs) if not episode_starts[i + 1]], dtype='int32')

        while True:
            count=0
            while True:
                # randomly sample a batch of indices :
                indices_1 = np.random.choice(all_indices, batch_size, replace=False)
                indices_2 = np.random.choice(all_indices, batch_size, replace=False)                
                
                same_action = np.prod( actions[indices_1]==actions[indices_2], axis=1 ) # deal with non scalar actions, for instance a=[0,1]
                diff_reward = 1*(rewards[indices_1+1]!=rewards[indices_2+1])
                same_action_diff_reward = same_action*diff_reward

                # check if there is at least one same_action and different reward, i.e (a1 = a2) & (r1+1 ≠ r2+1)
                if np.any(same_action_diff_reward):
                    break # break loop if there is, and yield data            
                #if not, repeat and count
                count=count+1
                if count>10:
                    # raise exception if counting too much, typically if the data is weird or the batch_size too small
                    raise Exception('Waiting too much to find a convenient batch with at least one (a1 = a2) & (r1+1 ≠ r2+1),\n\
                                     please consider increasing the batch size or changing the dataset')

            X = [ observations[indices_1], observations[indices_1+1], observations[indices_2], observations[indices_2+1],
                  same_action, same_action_diff_reward ]
            Y = [np.zeros([batch_size,1]) for i in range(4)] # unused here because this is unsupervised learning, yet it is needed
            yield X,Y

    def training_generator2(observations, actions, rewards, episode_starts, batch_size, steps_per_epoch):
        """
        data generator for training the model batch by batch, for each epochs batches are a split of a permutation of all indices
         this means that with this generator all data are used once at each epoch
        yields : (X, Y) tuple for training
        """
        n_obs = observations.shape[0] - 1

        if batch_size > n_obs-1 or batch_size < 1:
            raise ValueError('Wrong batch size')

        # indices for all time steps where the episode continues (to be able to count ∆s)
        all_indices_1 = np.array([i for i in range(n_obs) if not episode_starts[i + 1]], dtype='int32')
        all_indices_2 = np.copy(all_indices_1)

        # define a funciton to create vectors of indices of exactly the right size
        def compute_batches(indices):
            # first shuffle indices in place
            np.random.shuffle(indices)
            if batch_size*steps_per_epoch <= len(indices):
                batches = indices[0:batch_size*steps_per_epoch]
            else :
                batches = indices[0:batch_size*(steps_per_epoch-1)]
                last_batch = np.random.choice(indices[batch_size*(steps_per_epoch-1):], batch_size)
                batches = np.concatenate([batches, last_batch])
            return batches

        while True:           
            # create vectors of indices of exactly the right size
            batches_1 = compute_batches(all_indices_1)
            batches_2 = compute_batches(all_indices_2)

            # compute same_action and same_a_diff_r boolean vectors
            same_action = np.prod( actions[batches_1]==actions[batches_2], axis=1 ) # deal with non scalar actions, for instance a=[0,1]
            same_reward = 1*(rewards[batches_1+1]==rewards[batches_2+1])
            same_a_diff_r = same_action*(1-same_reward)

            ind_range = np.arange(len(same_a_diff_r)) # usefull when a batch hasn't any same_a_diff_r, see in for loop

            for step in range(steps_per_epoch):

                inds = np.arange(step*batch_size,(step+1)*batch_size)

                # sometimes, same_a_diff_r_batch happens to be full of zeros
                # in that case, we take a random other pair of index where it not zero and add it in this batch
                if not any(same_a_diff_r[inds]):
                    # take a pair where we know same_a_diff_r=1
                    random_pair_index = np.random.choice(ind_range[same_a_diff_r==1],1)[0]
                    # insert it at 0 index in this batch
                    batches_1[inds[0]] = batches_1[random_pair_index]
                    batches_2[inds[0]] = batches_2[random_pair_index]
                    same_action[inds[0]] = same_a_diff_r[inds[0]] = 1
                    
                # finally, yield the batch
                X = [ observations[batches_1[inds]], observations[batches_1[inds]+1], observations[batches_2[inds]], observations[batches_2[inds]+1],
                      same_action[inds], same_a_diff_r[inds] ]
                Y = [np.zeros([batch_size,1]) for i in range(4)] # unused here because this is unsupervised learning, yet it is needed
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

        train_gen = Priors_model.training_generator2(observations, actions, rewards, episode_starts, batch_size, steps_per_epoch)
        #train_gen = Priors_model.training_generator(observations, actions, rewards, episode_starts, batch_size)

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


