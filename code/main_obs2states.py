#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Cressot Loic
    ISIR - CNRS / Sorbonne Université
    10/2018
""" 

"""
    main_obs2states : main script for computing states from observations

    This script trains and computes states representations from image observations, rewards and actions,
    using the Jonschkowski and Brock method from 2015 paper 'Learning State Representations with Robotic Priors'
"""

# System
import os, warnings, argparse
import random, datetime

# maths
import numpy as np

# repo
import tools
import qlearning as ql
import gym_env.agent as agent
import gym_env.policy as policy
import gym_env.run_round_bot as rrb
# import jonschkowski_priors2015 later to set the random seed before any keras module importation

warnings.filterwarnings('ignore') # ignore warnings


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=None)

    required = parser.add_argument_group('Required arguments')
    required.add_argument('-trd','--training_data', default='', help='Select training data file to load', required=True)

    parser.add_argument('-ted','--test_data', default='', help='Select testing data file to load')
    parser.add_argument('-ql','--qlearning',  action='store_true', default=False, help='Perform Q-learning state evaluation')
    parser.add_argument('-r','--recordto', default='', help='Select a file for recording computed states')
    parser.add_argument('-ne','--num_epochs', type=int, default=10, help='Number of training epochs')
    parser.add_argument('-lr','--learning_rate', type=float, default=1e-4, help='Number of training epochs')
    parser.add_argument('-reg','--l1_reg', type=float, default=1e-3, help='l1 regularizer')
    parser.add_argument('-sd','--state_dim', type=int, default=2, help='State dimensions')
    parser.add_argument('-bs','--batch_size', type=int, default=256, help='Batch size')
    parser.add_argument('-rs','--seed', type=int, default=None, help='Seed for random')
    parser.add_argument('-tu','--tanh_units', type=int, default=None, help='hidden tanh units, default is None (linear model)')
    parser.add_argument('-w','--weights', nargs=4, type=float, metavar=('wt','wp','wc','wr'), default=(1.0,5.0,1.0,5.0), help='Weights of loss components')
    parser.add_argument('-dis','--display', action='store_true', default=False, help='display plots')

    args = parser.parse_args()

    if args.recordto and args.recordto[-1]!='/':
        args.recordto = args.recordto+'/'

    if args.recordto:
        os.makedirs(args.recordto, exist_ok=True)

    if not (args.recordto or args.display) :
        raise Exception('\nPlease either record (-r) or display (-dis) results ! (or both) \n') 

    if args.qlearning and not args.recordto :
        raise Exception('\nPlease give a record folder (-r) for saving results of qlearning \n') 

    # init seed before importing anything from keras
    random.seed(args.seed if args.seed else datetime.datetime.now())
    import jonschkowski_priors2015

    # load data
    training_data = tools.load_data(args.training_data)      
    # create a model implementing the paper method
    jp_model = jonschkowski_priors2015.Priors_model(
                            obs_shape = list(training_data['observations'][0].shape), 
                            state_dim = args.state_dim,
                            learning_rate=args.learning_rate, 
                            l1_reg = args.l1_reg, 
                            loss_weights=args.weights,
                            noise_stddev=1e-6,
                            tanh_units = args.tanh_units
                           )

    # if args.qlearning is False, then run a simple learning and plotting of the representations
    #      you can also specify a test dataset to test the learnt representation afterward
    if not args.qlearning:
        # learn the model directly for every iterations
        tools.learn_states( training_data=training_data,
                            model = jp_model,                                      
                            num_epochs=args.num_epochs,
                            batch_size = args.batch_size,
                            recordto=args.recordto,
                            display=args.display,
                        )    

        if args.test_data:
            tools.compute_states( test_data=args.test_data,
                                model=jp_model,
                                recordto=args.recordto,
                                display=args.display,
                              )
    # else, we will run the qlearning experiment of the paper :
    #   after each learning iteration, we run 10 q fitted iteration learnings
    #       for each q learning, we test it on 20 episodes of 25 steps and average the sum of rewards
    #   plot representations every 5 learnings and average reward over each 20
    else:
        # learning iteration loop
        for learning_epoch in range(args.num_epochs):

            states = tools.learn_states(
                            training_data=training_data,
                            model = jp_model,                                      
                            num_epochs=1,
                            batch_size = args.batch_size,
                            recordto='',
                            display=args.display,
                        ) 

            # plot representations every 5 learnings
            if learning_epoch%5 == 0:
                tools.plot_representation(
                    states[1:],
                    # offset reward of 1 to match observations, and don't show previous rewards for episode_start steps
                    training_data['rewards'][:-1]*(training_data['episode_starts'][1:]==False),
                    name='Observation-State-Mapping for ' + str(learning_epoch+1) + 'learning epoch',
                    add_colorbar=True, 
                    state_dim=min(jp_model.state_dim,3),
                    plot_name='_train_'+str(learning_epoch+1),
                    recordto=args.recordto,
                    display=False,
                    )

            # qlearning loop
            for q_learning in range(10):

                # Perform fitted Q iterations and get states policy (perform on a subsampled list of states):
                qit = ql.Fitted_QIteration(n_rbf=args.nuclei_rbf, n_actions=n_actions)

                state_policy = qit.fit_sk( states, training_data['actions'], training_data['rewards'], 0.9, 40, recompute_mapping=True)

                state_policy = policy.Epsilon_Plug_policy(state_policy, env.controller.action_space_int, epsilon=0.1)

                rb_agent = agent.RoundBotAgent(state_policy)
                # train the agent in the env
                stats = rb_agent.run_in_env(env, args.n_tests, seed)
                # record all episodes rewards for this evaluation run and this model and this learning
                test_performance[i,cur_learning,cur_evaluation,:,:] = np.reshape(np.array(stats['rewards']),[args.n_tests, env_args['max_step']])
                # recod all ground truth
                ground_truth_as_array = np.array([ np.concatenate([np.array(s[0]).flatten(),np.array(s[1]).flatten()],0) for s in stats['ground_truth']])
                all_ground_truth[i,cur_learning,cur_evaluation,:,:,:] = np.reshape( ground_truth_as_array[:,0:gt_shape[1]], [args.n_tests, env_args['max_step'], gt_shape[1] ])
                # record all observation (observed states)
                if not exp==3: # in exp 3, observations are raw and too big to be kept
                    all_observations[exp][cur_learning,cur_evaluation,:,:,:] = np.reshape( np.array(stats['observations']),[args.n_tests, env_args['max_step'], int(np.prod(np.array(input_shapes[exp]),0)) ])

                if verbose:
                    print( 'Mean reward over '+ str(args.n_tests) + ' episodes for ' + experiences[exp] + ' : ' +\
                        str( np.mean( stats['reward_ep'].flatten() ) )+' \n' )








