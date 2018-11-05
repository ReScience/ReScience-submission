#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Cressot Loic & Merckling Astrid
    ISIR - CNRS / Sorbonne Université
    02/2018
"""

# math / ML
import numpy as np
from sklearn.decomposition import PCA

# Plotting
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def convert_actions_to_integers(data_actions, actions_mapping):
    """
        Takes actions and actions_mapping loaded from .npz data and returns the actions
        as a list of their integer representation. (used when actions are tuples for instance)
        Also returns the number of actions
    """
    # now return a list of actions represented as integers in range(0, num_actions)
    return [actions_mapping[tuple(a)] for a in data_actions], len(actions_mapping)


def plot_representation(states, rewards, colors=None, name='Learned State Representation', add_colorbar=True, state_dim=2, plot_dims=None,recordto="./models/",plot_name="",display=False, verbose=True):
    """
        Plot representations, perform PCA and projection if lower dimension.
        If specified, after PCA, project on given plot_dims dimensions. ex : plot_dims=(3,4) instead of default (1,2)
        Can also plot 2D (TODO 3D) vector field
        
        Parameters :
        ------------
        - states (np.ndarray) : states to plot
        - rewards (np.ndarray) : rewards corresponding to states
        - colors (np.array(ndata,3)) : array with RGB color as rows to color the points, else None
        - name (str) : plot name
        - add_colorbar (Bool) : add a color bar to plot
        - state_dim (int) : 2 or 3 for 2D / 3D plotting
        - plot_dims (tuple(int)) : dimensions to plot if the states dimensions is higher than the plotting dimension
        - recordto (str) : path of saving folder
        - plot_name (str)
        - display (bool) : display pot
        - verbose (bool)
    """
    assert(display or recordto) # return exception if this function is called without save nor display
    # perform PCA
    pca = PCA(n_components=states.shape[1])
    states_PCA = pca.fit_transform(states)
    # take PCA first components if state_dim < states.shape[1]
    if state_dim < states.shape[1]:
        states = states_PCA
    elif state_dim > states.shape[1]:
        raise ValueError('Error, cannot plot states on axis with higher dimensions than states')
    if verbose:
        print('\n pca variance ratios : {} \n'.format(str(pca.explained_variance_ratio_)))

    # check plot_dims parameter. If None or empty, fill it
    if not plot_dims:
        plot_dims = [i for i in range(1,state_dim+1) ]
    elif len(plot_dims) != state_dim:
        raise ValueError("Error, plot_dims tuple :" + str(plot_dims) + " must be of len state_dim :"+ str(state_dim))
    else:
        plot_dims = list(plot_dims)
    
    add_colorbar = add_colorbar and (colors is None)

    if add_colorbar:
        null_rewards = rewards==0
        neg_rewards = rewards<0
        pos_rewards = rewards>0

        states = np.concatenate( (states[null_rewards,:], states[neg_rewards,:], states[pos_rewards,:]) , 0 )
        rewards = np.concatenate( ( rewards[null_rewards], rewards[neg_rewards], rewards[pos_rewards]) , 0 )
        
        c=np.clip(rewards, -1, 1) # colors of the points

    elif colors is not None:
        neg_rewards = rewards<0
        other_rewards = rewards>=0

        states = np.concatenate( (states[other_rewards,:], states[neg_rewards,:]) , 0 )
        rewards = np.concatenate( (rewards[other_rewards], rewards[neg_rewards]) , 0 )        
        colors = np.concatenate( (colors[other_rewards,:], colors[neg_rewards,:]) , 0 )        
        c=colors
    else:
        c=None
    
    plt.ion() #to turn in interactive mode
    fig=plt.figure(name)
    plt.clf()
    #plt.hold(True)
    plt.axis('equal')   
    if state_dim==2:        
        plt.scatter(states[:, plot_dims[0]-1], states[:, plot_dims[1]-1], s=7, c=c, cmap='coolwarm', linewidths=0.1)        
        plt.xlabel('State dimension ' + str(plot_dims[0]) )        
        plt.ylabel('State dimension ' + str(plot_dims[1]))
        if add_colorbar:
           plt.colorbar(label='Reward')
    elif state_dim >=3:
        fg = fig.add_subplot(111, projection='3d')
        fg.scatter(states[:, plot_dims[0]-1], states[:, plot_dims[1]-1], states[:, plot_dims[2]-1], c=c, cmap='coolwarm', linewidths=0.1)
        fg.set_xlabel('State dimension ' + str(plot_dims[0]))
        fg.set_ylabel('State dimension ' + str(plot_dims[1]))
        fg.set_zlabel('State dimension ' + str(plot_dims[2]))

    else:
        raise Exception('Dimension error for plotting : '+str(state_dim))
    if recordto:
        fig.savefig(recordto+'representation'+plot_name+'.png')# save the figure to file   
    if display:
        plt.pause(0.0001)



def plot_observations(observations, name = 'Observation Samples',recordto="./models/",plot_name="",display=False, rows=8, cols=10, offset=0, show_all=False):
    if not (display or recordto): # return exception if this function is called without save nor display
        return ValueError('plot_observations must be called with either display or save parameter')
    if observations.dtype != np.uint8:
        raise TypeError('observations must be np.uint8')

    plt.ion() #to turn in interactive mode
    while offset < len(observations) - rows*cols:
        plt.figure(name)
        current_observations = observations[offset:offset+rows*cols]
        # reshape observations if images are encoded as lines
        if len(current_observations.shape)==2:
            imwidth = int(np.sqrt(current_observations.shape[1]/3))
            current_observations = np.reshape( current_observations, [current_observations.shape[0], imwidth, imwidth, 3])
        for i in range(rows*cols):
            plt.subplot(rows, cols, i+1)
            plt.imshow(current_observations[i], interpolation='nearest')
            plt.gca().invert_yaxis()
            plt.xticks([])
            plt.yticks([])            
        if not show_all:
            break
        ans = input('Enter any key to continue. q to quit : ')
        if ans=='Q' or ans=='q':
            plt.close()
            break
        offset += rows*cols
        plt.close()
    ## end while ##
    if recordto:
        plt.savefig(recordto+'observations'+plot_name+'.png')# save the figure to file 
    if display:
        plt.pause(0.0001)



def assert_npz(file):
    if not '.npz' in file:
        raise ValueError("File must be a .npz !")
    else:
        return True


def load_data(datafile):
    """
        Loads a .npz file containing squarred image observations in 'observations' scope
        Returns loaded data
    """    
    assert_npz(datafile)
    #loaded_data_npz = np.load('../data/'+datafile+'.npz')#simple_navigation_task_train
    loaded_data_npz = np.load(datafile)#simple_navigation_task_train
    # do this to convert loaded data to their initial type
    loaded_data = {k:loaded_data_npz[k] for k in loaded_data_npz.keys()}

    return loaded_data


def learn_states(   training_data,
                    model,                                      
                    num_epochs,
                    batch_size,
                    recordto='',
                    plot_dims=None,
                    display=True,
                    validation_ratio=0.1,
                    ):
    """
        - Loads data if training_data is a file path
        - Learns states for observations with model object
        - Displays results if asked
        - Returns learned states and loaded training_data
    """

    # Loads data if training_data is a file path
    if isinstance(training_data, str):
        training_data = load_data(training_data)

    # Plot / Record training data
    if recordto or display:
        plot_observations(  training_data['observations'], 
                            name="Observation Samples (Subset of Training Data) -- Simple Navigation Task",
                            recordto=recordto,
                            display=display,
                            )

    print('\nLearning a state representation ... \n')

    history = model.learn(observations=training_data['observations'], actions=training_data['actions'],
                rewards=training_data['rewards'], episode_starts=training_data['episode_starts'],
                batch_size=batch_size, num_epochs=num_epochs, verbose=True, validation_ratio=validation_ratio)

    # predict learned states
    training_states = model.phi(training_data['observations'])
    
    # Plot / Record learnt states
    if recordto or display:      
        # Plot / Record state representation
        plot_representation(training_states[1:],
                            # offset reward of 1 to match observations, and don't show previous rewards for episode_start steps
                            training_data['rewards'][:-1]*(training_data['episode_starts'][1:]==False),
                            name='Observation-State-Mapping Applied to Training Data -- Simple Navigation Task',
                            add_colorbar=True, 
                            state_dim=min(model.state_dim,3),
                            plot_dims=plot_dims,
                            plot_name='_train',
                            recordto=recordto,
                            display=display,
                            )

    # return learnt states along with traing history
    return training_states, history



def compute_states( test_data,
                    model,
                    recordto='', 
                    display=False,
                    plot_dims=None,
                   ):
    """
        Computes states from trained model
        - Loads data if test_data is a file path
        - Computes states with model object
        - Displays results if asked
        - Returns computed states
    """

    # Loads data if training_data is a file path
    if isinstance(test_data, str):
        test_data = load_data(test_data)
    else:
        test_data = test_data

    if display:
        print('\nDisplaying testing data ...\n')
        plot_observations(test_data['observations'],
                        name='Observation Samples (Subset of Test Data) -- Simple Navigation Task',
                        display=display,
                        )

    print('\nComputing states from learned representations ...\n')
    
    test_states = model.phi(test_data['observations'])


    plot_representation(test_states[1:],
                        # offset reward of 1 to match observations, and don't show previous rewards for episode_start steps
                        test_data['rewards'][:-1]*(test_data['episode_starts'][1:]==False),
                        name='Observation-State-Mapping Applied to Test Data -- Simple Navigation Task',
                        add_colorbar=True,
                        state_dim=min(model.state_dim,3),
                        plot_dims=plot_dims,
                        plot_name='_test',
                        display=display,
                        )

    if recordto:
        assert_npz(recordto)
        np.savez(recordto, data=test_data, states=test_states)

    return test_states



