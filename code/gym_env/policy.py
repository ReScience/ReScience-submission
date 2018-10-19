#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Cressot Loic & Merckling Astrid
    ISIR - CNRS / Sorbonne Université
    02/2018

    This module defines policy classes to control the agent
""" 

from collections import defaultdict
import random

class Policy(object):
    """
    Abstract Policy class for the agent : action = p(observation)
    """
    def __init__(self, action_space):
        self.action_space=action_space
        if type(self) is Policy:
            raise NotImplementedError('Cannot instantiate this abstract class')
    def __call___(self,observation):
        raise NotImplemented
    def update(self, observation, new_observation, old_action, action, reward, done):
        """Default and constant policis have no updates"""
        return
    def convert_action(action):
        """
        Convert actions to either int or tuple of ints when they are of other similar types
        """    
        type_action = type(action)
        if type_action!=type(int()) and type_action!=type(tuple()):
            try :
                action=int(action) # try converting to int first
            except:
                action=tuple(action) # convert lists to tuple for controller
        return action
    def sample(self):
        return Policy.convert_action( self.action_space.sample() )

class Plug_policy(Policy):
    """
    Used to plug a function to the Policy interface
    """
    def __init__(self, plug_policy, action_space=None):
        super(Plug_policy,self).__init__(action_space)
        self.plug_policy=plug_policy
    def __call__(self, observation):
        return Policy.convert_action( self.plug_policy(observation) )

class Compose_policy(Plug_policy):
    """
    Compose a policy with a transformation: action = p( transformation( observation ) )
    """
    def __init__(self, policy, transformation, action_space=None):
        super(Compose_policy,self).__init__(policy, action_space)
        self.transformation = transformation
    def __call__(self, observation):
        return super(Compose_policy,self).__call__( self.transformation(observation)  )

class Mapping_policy(Plug_policy):
    """
    Maps existing policy's actions to another type of actions (ex : int to tuples) : action = mapping( p( observation ) )
    """
    def __init__(self, policy, mapping, action_space=None):
        super(Mapping_policy,self).__init__(policy, action_space)
        self.mapping = mapping
    def __call__(self, observation):        
        return self.mapping[ super(Mapping_policy,self).__call__(observation) ]


class Random_policy(Policy):
    def __init__(self,action_space):
        super(Random_policy,self).__init__(action_space)
    def __call__(self,observation):
        return Policy.convert_action( self.action_space.sample() )
 

class Greedy_policy(Policy):
    def __init__(self,action_space):
        super(Greedy_policy,self).__init__(action_space)
        self.Q = defaultdict(lambda: np.zeros(action_space.n))
    def __call__(self,observation):
        return Policy.convert_action( np.argmax(self.Q[observation]) )


class Epsilon_greedy_policy(Greedy_policy):
    def __init__(self,action_space,epsilon):
        super(Epsilon_greedy_policy,self).__init__(action_space)
        self.epsilon=epsilon
    def __call__(self,observation):
        """ random sample epsilon % of the time """        
        if random.random() > self.epsilon:
            return super(Epsilon_greedy_policy,self).__call__(observation)
        else:
            return Policy.convert_action( self.action_space.sample() )


class Constant_policy(Policy):
    def __init__(self,constant_action, action_space=None):
        super(Constant_policy,self).__init__(action_space)
        self.constant_action = constant_action
    def __call__(self, observation):
        return Policy.convert_action( self.constant_action )


class Epsilon_Plug_policy(Plug_policy):
    """ plug a policy 1-epsilon of the time and random sample epsilon of the time
    """
    def __init__(self, policy, action_space, epsilon):
        super(Epsilon_Plug_policy,self).__init__(policy, action_space)
        self.epsilon=epsilon
    def __call__(self, observation):
        """ random sample epsilon % of the time """
        if random.random() > self.epsilon:
            return super(Epsilon_Plug_policy,self).__call__(observation)
        else:
            return Policy.convert_action( self.action_space.sample() )
