#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Cressot Loic & Merckling Astrid
    ISIR - CNRS / Sorbonne Université
    02/2018

    This module rewrite the gym.make function so as to return wrapper containing the number of actions of the environment
"""

import gym
from gym.wrappers.time_limit import TimeLimit


class ActionWrapper(TimeLimit):
    """
    Environment Wrapper for dealing with actions.
    Implementation of properties : num_actions, actions_mapping
    """

    def __init__(self, env):
        """ Copy constructor
            Param : TimeLimit env
        """
        self.env = env.env
        self.action_space = env.action_space 
        self.observation_space = env.observation_space 
        self.reward_range = env.reward_range 
        self.metadata = env.metadata 
        self._max_episode_seconds = env._max_episode_seconds 
        self._max_episode_steps = env._max_episode_steps 
        self._elapsed_steps = env._elapsed_steps 
        self._episode_started_at = env._episode_started_at
        self._warn_double_wrap()
    
    @property
    def num_actions(self):
        """
        Get number of actions of gym space
        """
        name=type(self.action_space).__name__
        if name == 'Discrete':
            return self.action_space.n# +1 to count the 0
        elif name =='MultiDiscrete':
            prod=1
            for d in self.action_space.nvec:
                prod *= d
            return int(prod)
        elif name=='Box': # if actions are continous, return the number of dimensions of the action space
            return self.action_space.shape[0]
        else:
            raise Exception('Unknown action_space class name. Should be either Discrete, MultiDiscrete or Box')

    @property
    def actions_mapping(self):
        name=type(self.action_space).__name__
        if name == 'Discrete':
            return {i:i for i in range(self.env.action_space.n)}
        elif name =='MultiDiscrete':
            r=[[]]
            for x in self.env.action_space.nvec:
                t = []
                for y in list(range(x)):
                    for i in r:
                        t.append(i+[y])
                r = t
            return {tuple(r[i]): i for i in range(len(r))}

    @property
    def max_steps(self):
        return self._max_episode_steps

    @property
    def controller(self):
        return self.env.controller

    @max_steps.setter
    def max_steps(self, ms):
        self._max_episode_steps = ms

    def message(self,message):
        self.env.message(message)


