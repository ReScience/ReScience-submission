#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Cressot Loic & Merckling Astrid
    ISIR - CNRS / Sorbonne Université
    02/2018

    This simple module allows to run a policy in the round_bot environment
"""
import sys, os

import gym
import gym_round_bot
from gym_round_bot.envs import round_bot_env

path = os.path.abspath(__file__)
dir_path = os.path.dirname(path)
sys.path.insert(0, dir_path)

from action_wrapper import ActionWrapper

import agent as agent
import policy as policy


def make_round_bot_env( controller,
                        world={'name':'square','size':[20,20]},
                        texture='minecraft+',
                        obssize=(16,16),
                        winsize=None,
                        global_pov=None,
                        visible=False,
                        perspective = True,
                        multiview=False,
                        focal=65.0,
                        max_step=100,
                        crash_stop=False,
                        random_start=True,
                        reward_count_stop = False,
                        reward_stop=False,
                        normalize_observations=True,
                        normalize_rewards=True,
                        position_observations=False,
                        observation_transformation = None,
                        distractors=False,
                        sandboxes=False,
                        trigger_button=False,
                        robot_diameter=2,
                        ):

    # set loading variables before creating env with make method
    round_bot_env.set_metadata(
                world=world,
                texture=texture,
                obssize=obssize,
                winsize=winsize,
                controller=controller,
                global_pov=global_pov,
                visible=visible,
                perspective = perspective,
                multiview=multiview,
                random_start=random_start,
                focal=focal,
                crash_stop=crash_stop,
                reward_count_stop = reward_count_stop,
                reward_stop=reward_stop,
                normalize_observations=normalize_observations,
                normalize_rewards=normalize_rewards,
                position_observations=position_observations,
                observation_transformation = observation_transformation,
                distractors=distractors,
                sandboxes=sandboxes,
                trigger_button=trigger_button,
                robot_diameter=robot_diameter,
              )

    # create env
    env = gym.make('RoundBot-v0')

    # wrapp it in action_wrapper
    env = ActionWrapper(env)
    # exceptionnaly change private attribute of TimeLimit env object
    env._max_episode_steps = max_step
    
    return env


def run_round_bot(  policy,
                    controller,
                    transformation=None,
                    n_ep=1,
                    winsize=None,
                    obssize=(16,16),
                    visible=False,
                    global_pov=None,
                    focal=65.0,
                    ortho_pers=False,
                    multiview=None,
                    seed=None,
                    world={'name':'square','size':[20,20]},
                    texture='minecraft+',
                    max_step = 100,
                    crash_stop=False,
                    reward_count_stop = False,
                    reward_stop=False,
                    n_seq=1,
                    normalize_observations=True,
                    normalize_rewards=True,
                    position_observations=False,
                    observation_transformation = None,
                    robot_diameter=2,
                  ):
    """
    Runs a policy in the round_bot environment

    Parameters
    ----------
    - winsize (tuple of (int,int)) : window's size
    - global_pov (int) : (height of global point of view) or None
    - controller : Controller dict. ex: {'name':'Theta','speed':args.speed,'dtheta':args.dtheta},
    - n_seq : max number of previous observations to take for current state representation 
    """
    
    env = make_round_bot_env( world=world,
                            texture=texture,
                            obssize=obssize,
                            winsize=winsize,
                            controller=controller,
                            global_pov=global_pov,
                            visible=visible,
                            perspective = not ortho_pers,
                            multiview=multiview,
                            focal=focal,
                            max_step=max_step,
                            crash_stop=crash_stop,
                            reward_count_stop = reward_count_stop,
                            reward_stop=reward_stop,
                            normalize_observations=normalize_observations,
                            normalize_rewards=normalize_rewards,
                            position_observations=position_observations,
                            observation_transformation = observation_transformation,
                            robot_diameter=robot_diameter
                            )

    # TODO : check is Policy implements the policy.Policy interface
    # if not issubclass(type(policy), policy.Policy):
    #     # if not, create a plug policy object to use the given policy in agent
    #     policy=policy.Plug_policy(policy)    

    # create an agent with this policy and transformation
    #   recall policy is learnt from states, whereas transformation is learnt from observation
    #   so they are two different things
    rb_agent = agent.RoundBotAgent(policy, transformation)
    # train the agent in the env
    stats = rb_agent.run_in_env(env, n_ep, seed, n_seq=n_seq)

    # return stats dictionnary
    return stats
    