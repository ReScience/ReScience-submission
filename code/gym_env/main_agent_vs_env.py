#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Cressot Loic & Merckling Astrid
    ISIR - CNRS / Sorbonne Université
    02/2018

    This script allows to connect an agent with a given policy in a given environment
    If the policy is a learning policy, it will be updated, if not, it will stay constant
    Possibility to save all observation, actions and rewards, and also command line arguments
        ( the saved command line arguments allow to keep the used controler and environnement )

    Usage example (to update in README)
    -------------

    ```Bash
    python3 main_agent_vs_env.py -c XZ -xzt 2 2 1 --world square --texture graffiti -ms 50 -nep 20 -nr 0.1 -obs 16 16 --speed 5 -f 70 -r ../../data/ 
    ```
    generation of observations for : 
    - XZ controler
    - x range : 2 
    - z range : 2
    - world : square
    - texture : graffiti
    - episode maximum steps : 100
    - number of episodes : 20
    - observation size = 16 x 16
    - speed : 5 units
    - focal angle : 70
    - noise ratio : 10 %
    - record directory : ../../data
    - other arguments set to default value  

"""

import argparse

import gym
import gym_round_bot
from gym_round_bot.envs import round_bot_env

import numpy as np
import matplotlib.pyplot as plt
import os

import agent
import policy
from action_wrapper import ActionWrapper


if __name__ == '__main__':

    #################### PARSE COMMAND LINE OPTIONS ###########################

    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('-id','--env_id', default='RoundBot-v0', help='Select the environment to run')
    parser.add_argument('--name', default='', help='suffix name to automatic name')
    parser.add_argument('-c','--controller', default='XZ', choices=['XZ', 'Theta', 'Theta2','XZF','XZc','XZca'], help='Select the agent\'s controller')
    parser.add_argument('--world_name', default='square', choices=['square', 'square_1wall'], help='Select the world')
    parser.add_argument('--texture', default='minecraft', choices=['graffiti', 'minecraft', 'minecraft+', 'colours'], help='Select the bricks texture')
    parser.add_argument('-p','--policy', default='Random_policy', choices=['Random_policy','Constant_policy'], help='Select the policy of the agent to run')
    parser.add_argument('-xzt','--xztrange', nargs=3, type=int, default=(1,1,1), help='xrange, zrange and theta range')
    parser.add_argument('-ms','--max_step',  type=int,default=200, help='Number of steps in a trajectory')
    parser.add_argument('-f','--focal',  type=float,default=65.0, help='Focal angle of camera')
    parser.add_argument('-w','--winsize', nargs=2, type=int, metavar=('w','h'), default=None, help='Size of rendering window')
    parser.add_argument('-nr','--noise_ratio', type=float, default=0.0, help='Ratio of speed in additive gaussian noise stdv')
    parser.add_argument('--fixed_point', nargs=2, type=int, default=[0,0], help='Fixed point for XZF controller')
    parser.add_argument('-obs','--obssize', nargs=2, type=int, metavar=('w','h'), default=(16,16), help='Size of environment observations')
    parser.add_argument('-ml','--multiview', nargs='*', type=float, default=None, help='Multi-view angles, nargs = *')
    parser.add_argument('-nep','--n_ep', type=int,default=26, help='number of episodes/trajetories')
    parser.add_argument('-rt','--recordto', default='', help='Path without .npz to record oupputs (params will be added to that name)')
    parser.add_argument('-gp','--global_pov', type=float, default=None, help='global_pov height (centered), default is None')
    parser.add_argument('-agp','--auto_global_pov', action='store_true', default=False, help='automatic computation of global_pov')
    parser.add_argument('-v','--verbose', action='store_true', default=False, help='verbose')
    parser.add_argument('-vis','--visible', action='store_true', default=False, help='set window visible, 10 times slower')
    parser.add_argument('-ort','--orthogonal', action='store_true', default=False, help='render with orthogonal perspective')
    parser.add_argument('-dis','--distractors', action='store_true', default=False, help='add visual distractors')
    parser.add_argument('-sb','--sandboxes', action='store_true', default=False, help='add sandboxes on the ground that slow down the robot')
    parser.add_argument('-tb','--trigger_button', action='store_true', default=False, help='add a trigger button on the ground that triggers some effect')
    parser.add_argument('--world_size', nargs=2, type=float, default=[20,20], help='width an depth of the world')
    parser.add_argument('--speed', type=float, default=3.0, help='agent\'s speed')
    parser.add_argument('--dtheta', type=float, default=7.0, help='rotation angle')
    parser.add_argument('-rd','--robot_diameter', type=float, default=2.0, help='robot diameter')
    args = parser.parse_args()



    #################### READ AND VERIFY COMAND LINE OPTIONS ###########################

    if args.auto_global_pov:
        args.global_pov = True
        
    if args.recordto and not os.path.isdir(os.path.dirname(args.recordto)):
            raise Exception('Error: Directory '+os.path.dirname(args.recordto)+' does not exist !')


    n_ep = int(args.n_ep) #trajectories of size max_step
    max_step = int(args.max_step)
    winsize = args.winsize
    obssize = args.obssize



    ####################  CONFIGURE ENVIRONMENT ###########################

    if args.env_id=='RoundBot-v0':
               
        from gym_round_bot.envs import round_bot_controller as rbc

        if args.global_pov!=True:
            global_pov = (0.0,args.global_pov,0.0) if args.global_pov and args.global_pov!=0.0 else None
        else:
            global_pov = True
        
        controller=rbc.make(name=args.controller,
                            speed=args.speed,
                            dtheta=args.dtheta,
                            fixed_point=list(args.fixed_point),
                            xzrange=args.xztrange[0:2],
                            thetarange=args.xztrange[2],
                            noise_ratio=args.noise_ratio)
                
        # set env metadata
        round_bot_env.set_metadata(
                world={'name':args.world_name,'size':tuple(args.world_size)},
                robot_diameter=args.robot_diameter,
                texture=args.texture,
                obssize=args.obssize,
                winsize=args.winsize,
                controller=controller,
                global_pov=global_pov,
                visible=args.visible,
                perspective =not args.orthogonal,
                multiview=args.multiview,
                focal=args.focal,
                crash_stop=False,
                reward_stop=False,
                random_start= (args.policy=='Random_policy'),
                reward_count_stop = False,
                distractors = args.distractors,
                sandboxes= args.sandboxes,
                trigger_button= args.trigger_button,
              )
        
    # create env 
    env = gym.make(args.env_id)
    # wrap in in ActioWrapper
    env = ActionWrapper(env)
    # set manually the maximum episodes steps in the env
    env.max_steps = max_step


    ####################  RUN AGENT IN ENVIRONNEMENT ###########################

    # create a policy and an agent
    if args.policy=='Random_policy':
        policy=policy.Random_policy(env.action_space)
    elif args.policy=='Constant_policy':
        #policy=policy.Constant_policy((0, 1)) # constant rotation policy  for Theta2
        policy=policy.Constant_policy((0, 0)) # constant rotation policy for Theta

    if args.env_id=='RoundBot-v0':
        myagent=agent.RoundBotAgent(policy)
    else:
        myagent=agent.Agent(policy)

    # run (train) agent with this policy
    stats=myagent.run_in_env(env,n_ep,keep_screen=True,new_winsize=None)


    ################  SAVE STATS IF ASKED ################

    if args.recordto:     

        # create the name of the data file :

        dtheta_str = '_dt'+ str( int(args.dtheta)) if 'Theta' in args.controller else ''
        gp_str = '_gp-' + str(args.global_pov) if args.global_pov else ''
        ml_str = '_ml' + str(len(args.multiview)) if args.multiview else ''
        tb_str = '_tb' if args.trigger_button else ''
        sb_str = '_sb' if args.sandboxes else ''
        dis_str = '_dis' if args.distractors else ''

        data_name = args.name +\
                    '_' + args.controller +\
                    '_'+ str( int(args.xztrange[0])) + 'x' + str( int(args.xztrange[1])) + 'x' + str( int(args.xztrange[2])) +\
                    '_'+ str( args.obssize[0]) +'x'+ str( args.obssize[1]) +\
                    '_'+ str( args.max_step) + 'x' + str( args.n_ep ) +\
                    '_s'+ str( int(args.speed)) +\
                    dtheta_str +\
                    '_f'+ str( int(args.focal) ) +\
                    '_'+ args.world_name +\
                    '_ws'+ str( int(args.world_size[0])) +'x'+ str( int(args.world_size[1])) +\
                    gp_str +\
                    '_' + args.texture +\
                    '_noise' + '{:1.2f}'.format(args.noise_ratio) +\
                    ml_str +\
                    sb_str +\
                    tb_str +\
                    dis_str
                    #'_ort-'+ str( args.orthogonal

        print('recording to : '+ args.recordto + data_name + '.npz')

        if len(np.asarray(stats['actions']).shape)==1:            
            stats['actions']=np.reshape(stats['actions'],(-1,1))
            if stats['discrete_actions']:
            # if stats['discrete_actions'] is True, then actions are discrete and stats['actions_int'] exists
                stats['actions_int']=np.reshape(stats['actions_int'],(-1,1))
        # save stats and also command line args to be able to recreate a controler and an environnement
        np.savez( args.recordto+ data_name+'.npz', args=vars(args), **stats )

        print('Done.')
