#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Cressot Loic & Merckling Astrid
    ISIR - CNRS / Sorbonne Université
    02/2018
""" 

import numpy as np
import time
import scipy.misc
import random
import os, sys
from tqdm import tqdm


path = os.path.abspath(__file__)
dir_path = os.path.dirname(path)

sys.path.insert(0, dir_path+'/../')
import tools 


class Agent(object):
    """
    Classic Agent
    """
    def __init__(self, policy, transformation=None):
        """
        paramters :
        -----------
        - policy : policy  to apply
        - transformation : transformation to be applied on raw observations. Default is None
        """
        self._policy=policy
        self._transformation = transformation # observation mapping to features

    def act(self, observation, reward, done, previous_action=None):
        action=self._policy(observation)

    def message(self, env, reward, done, episode, step):
        """
        Send message to environnement if implemented, default is not implemented, see subclasses
        """
        return

    def run_in_env(self, env, n_ep, keep_screen=False, new_winsize=None, seed=None, n_seq=1):
        """
        Run `n_ep` episodes of maximum `env._max_episode_steps` length in environment `env`
        Keep trace of everything and return it

        Note : `action[t]` is the action performed AFTER observing `observation[t]` and `reward[t]`
        
        Parameters
        ----------
        - env (gym.env object) : AI gym env
        - n_ep (int) : number of episodes
        - keep_screen (bool): keep screen shots of environment
        - new_winsize ( tuple(int,int) ): specify different size for kept screen shots
        - seed (int or float): seed for env.seed
        - n_seq (int) : number of previous (observations,action,rewards) to feed to transformation

        Returns
        -------
        - stats (dict) : dictionnary containing recorded array-like variables :
            `reward_ep`, `rewards`, `observations`, `actions`, `num_actions`, `episode_starts`, `screens`, `ground_truth`,
            `discrete_actions`, `actions_int`, `centered_reduced_actions
        """

        print('\nRunning agent in environment for ' + str(n_ep) + ' epochs of max steps ' + str(env._max_episode_steps) + '\n')
        
        # init all used variables :
        reward = 0.0
        done = False
        max_list_length = n_ep*env._max_episode_steps
        # initialize lists with their maximum length to gain computation speed and remove unused cells afterwards
        screens=[0,]*max_list_length
        observations=[0,]*max_list_length
        rewards=[0,]*max_list_length
        actions=[0,]*max_list_length
        episode_starts=[0,]*max_list_length
        reward_ep=np.zeros((n_ep,1)) # rewards per episode
        ground_truth = [0,]*max_list_length 
        
        #### Prepare main loop
        raw_ob = env.reset()
        if keep_screen:            
            winsize=env.render(mode='rgb_array').shape

        # define some functions and shorcuts
        updatePolicy = self._policy.update

        #### Start main loop

        env.seed(seed)
        t=0
        new_episode=True
        t1= time.time()  

        for i in tqdm(range(n_ep), desc='agent runs in env'): # tqdm prints loop progress
            raw_ob = env.reset()
            if n_seq==1:
                transf_ob = self._transformation(raw_ob) if self._transformation else raw_ob
            else:    
                obs_a_r = [raw_ob, self._policy.sample(), 0.0]
                transf_ob = self._transformation(obs_a_r) if self._transformation else obs_a_r
            j=1; done=False; old_action=None
            t_ep=t

            while not done:
                # record raw observation
                observations[t] = raw_ob

                # record ground truth
                ground_truth[t] = env.env.ground_truth

                ### RL #################### -->

                # choose action from policy with current (previous) observation (transformed as states)
                action = self._policy(transf_ob)

                # perform a step in the environnment to observe the results of this action
                raw_ob, reward, done, _ = env.step(action)
                if n_seq==1:
                    new_transf_ob = self._transformation(raw_ob) if self._transformation else raw_ob
                else:  # pass to transformation the n_seq last obs, actions and rewards
                    t_ini = max(t_ep,t-n_seq)
                    obs_a_r = [ observations[t_ini:t+1]+[raw_ob], actions[t_ini:t+1]+[action], rewards[t_ini:t+1]+[reward] ]
                    new_transf_ob = self._transformation(obs_a_r) if self._transformation else obs_a_r

                # update learning parameters in policy accordingly
                updatePolicy(transf_ob,new_transf_ob,old_action,action,reward,done)

                ############################ <--

                # keep screens as line numpy arrays if asked
                if keep_screen:
                    rnd=env.render(mode='rgb_array')
                    if new_winsize:
                        rnd=scipy.misc.imresize(rnd, (new_winsize[0],new_winsize[1],3))                   
                    screens[t] = rnd
                # keeps the rest of parameters for recording                
                rewards[t] = reward
                actions[t] = action
                episode_starts[t] = new_episode
                reward_ep[i]+=reward
                # message if implemented
                self.message(env, reward, done, i, j)
                # increase counters and save old variables
                t+=1; j+=1
                transf_ob=new_transf_ob;  old_action=action;  new_episode = done
                
                # print expected computation time
                if (t==100 and env._max_episode_steps > 100) or (env._max_episode_steps <= 100 and t==env._max_episode_steps-1) :
                    t2= time.time()
                    div = 100 if (env._max_episode_steps > 100) else (env._max_episode_steps-1)
                    # maximum expected computation time is  : mean_step_time * max_remaining_steps = (t2-t1)/div * (n_ep*env._max_episode_steps -t+1)
                    if env._max_episode_steps > 100:
                        print('Maximum expected computation time: '+str( (t2-t1)/div * ( n_ep*env._max_episode_steps -t+1 ) ) +' s')                
                    elif env._max_episode_steps <= 100:
                        print('Maximum expected computation time: '+str( (t2-t1)/div * (n_ep*env._max_episode_steps -t+1)) +' s')

        t2= time.time()
        print('Final computation time: '+str(t2-t1) +' s\n')
        # create stats returned dict. [0:t] to only return used cells
        stats = {'reward_ep':reward_ep,
                'rewards':rewards[0:t],
                'observations':observations[0:t],
                'actions':actions[0:t],                
                'num_actions':env.num_actions,
                'episode_starts':episode_starts[0:t],
                'screens':screens[0:t],
                'ground_truth':ground_truth[0:t]
                }

        # if actions are discrete, add actions with integer coding
        if env.controller.discrete:
            stats['discrete_actions'] = True       
            ## create a converted list of actions coded as integers from 0 to num_actions instead of tuples or other
            if type(actions[0]) != type(np.int64()) and type(actions[0]) != type(int()):
                # sometimes actions are not of type int, so try to create a copy of actions with this type :
                actions_mapping = env.actions_mapping            
                stats['actions_int'] = tools.convert_actions_to_integers(actions[0:t], actions_mapping)
            else:
                # actions are already of type int
                stats['actions_int'] = actions[0:t]            
        # if actions are continuous, set discrete_actions to False, and add centered and reduced actions
        else:
            stats['discrete_actions'] = False
            stats['centered_reduced_actions'] = env.controller.center_reduce_actions(actions[0:t])

        return stats


class RoundBotAgent(Agent):
    """
    Agent designed for Round Bot Environment
    """
    def __init__(self, policy, transformation=lambda X: X):
        super(RoundBotAgent,self).__init__(policy, transformation)
        self.cum_reward=0
        
    def act(self, observation, reward, done, previous_action=None):
        action=self._policy(observation)

    def message(self, env, reward, done, episode, step):
        """
        Send message to Round_Bot environnement
        """
        self.cum_reward += reward
        message = 'Rew : %d ep : %d %d' % (self.cum_reward, episode, step)
        env.unwrapped.message(message)
    

