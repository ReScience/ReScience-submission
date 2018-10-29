#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Cressot Loic
    ISIR - CNRS / Sorbonne Université
    10/2018
""" 

"""
    main : main script for computing states from observations and learn policies with q fitted iteration

    This script trains and computes states representations from image observations, rewards and actions,
    using the Jonschkowski and Brock method from 2015 paper 'Learning State Representations with Robotic Priors'
    You can test the learned representation in RL with our q fitted iteration implementation, as in the paper
"""

# System
import os, warnings, argparse
import random, time

# maths
import numpy as np
# plotting
import matplotlib.pyplot as plt

# repo
import tools
import jonschkowski_priors2015


warnings.filterwarnings('ignore') # ignore warnings



parser = argparse.ArgumentParser(description=None)

required = parser.add_argument_group('Required arguments')
required.add_argument('-trd','--training_data', default='', help='Select training data file to load', required=True)

parser.add_argument('-ted','--test_data', default='', help='Select testing data file to load')
parser.add_argument('-ql','--qlearning',  action='store_true', default=False, help='Perform Q-learning state evaluation')
parser.add_argument('-r','--recordto', default='', help='Select a file for recording computed states')
parser.add_argument('-ne','--num_epochs', type=int, default=25, help='Number of training epochs')
parser.add_argument('-lr','--learning_rate', type=float, default=1e-4, help='Number of training epochs')
parser.add_argument('-reg','--l1_reg', type=float, default=1e-3, help='l1 regularizer')
parser.add_argument('-sd','--state_dim', type=int, default=2, help='State dimensions')
parser.add_argument('-bs','--batch_size', type=int, default=256, help='Batch size')
parser.add_argument('-rs','--seed', type=int, default=None, help='Seed for random')
parser.add_argument('-tu','--tanh_units', type=int, default=None, help='hidden tanh units, default is None (linear model)')
parser.add_argument('-w','--weights', nargs=4, type=float, metavar=('wt','wp','wc','wr'), default=(1.0,5.0,1.0,5.0), help='Weights of loss components')
parser.add_argument('-dis','--display', action='store_true', default=False, help='display plots')
parser.add_argument('-vis','--visible_train', action='store_true', default=False, help='Set the RL to visible')
parser.add_argument('-v','--verbose', action='store_true', default=False, help='Verbose')

args = parser.parse_args()

if args.recordto and args.recordto[-1]!='/':
    args.recordto = args.recordto+'/'

if args.recordto:
    os.makedirs(args.recordto, exist_ok=True)

if not (args.recordto or args.display) :
    raise Exception('\nPlease either record (-r) or display (-dis) results ! (or both) \n') 

if args.qlearning and not args.recordto :
    raise Exception('\nPlease give a record folder (-r) for saving results of qlearning \n') 

# init seed (with np.random, not random) before creating anything with keras, this allows to debug with deterministic seed in args.seed
np.random.seed(args.seed if args.seed else int(time.time()) )

total_loss_record = None
temporalcoherence_loss_record = None
proportionality_loss_record = None
causality_loss_record = None
repeatability_loss_record = None

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
    _ ,history = tools.learn_states(training_data=training_data,
                                    model = jp_model,                                      
                                    num_epochs=args.num_epochs,
                                    batch_size = args.batch_size,
                                    recordto=args.recordto,
                                    display=args.display,
                                    )

    total_loss_record = history.history['loss']
    temporalcoherence_loss_record = history.history['weighted_temporalcoherence_loss']
    proportionality_loss_record = history.history['weighted_proportionality_loss']
    causality_loss_record = history.history['weighted_causality_loss']
    repeatability_loss_record = history.history['weighted_repeatability_loss']

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

    import fqiteration as fqi
    import gym_env.agent as agent
    import gym_env.policy as policy
    import gym_env.run_round_bot as run_round_bot

    # gym round bot
    from gym_round_bot.envs import round_bot_controller

    # those are the argument parameters used for setting the gym environment where we got the training data
    # we can use these arguments to test a learned policy in this same environment
    env_args = training_data['args'].item(0)
    # integer version of the actions
    actions_int = training_data['actions_int'].item(0)
    n_actions = training_data['num_actions'] # number of actions

    # n_qlearnings = 10 # as in the paper, the number of different q iteration learnings
    # n_test_episodes = 20 # as in the paper, the number of episodes to test the learnt policy
    # n_test_steps = 25 # as in the paper, the number of steps per env policy test iteration
    # n_rbf = 100 # as in the paper, the number of rbf kernels used in q fitted iteration

    n_qlearnings = 1 # as in the paper, the number of different q iteration learnings
    n_test_episodes = 20 # as in the paper, the number of episodes to test the learnt policy
    n_test_steps = 25 # as in the paper, the number of steps per env policy test iteration
    n_rbf = 100 # as in the paper, the number of rbf kernels used in q fitted iteration

    # create contained array to save all test performances :
    #   array size : args.num_epochs learning steps X 10 q learnings X 20 test episodes X 25 steps per episode
    test_performance = np.zeros([args.num_epochs, n_qlearnings, n_test_episodes, n_test_steps])

    ### create the gym environment for testing the learnt policies (the same as the training_data env)
    # the controller
    controller=round_bot_controller.make(name=env_args['controller'],
                                speed=env_args['speed'],
                                dtheta=env_args['dtheta'],
                                fixed_point=list(env_args['fixed_point']),
                                xzrange=env_args['xztrange'][0:2],
                                thetarange=env_args['xztrange'][2],
                                noise_ratio=env_args['noise_ratio'],
                                int_actions=True)

    # create a function for computing states form the observations given by the env
    mean_obs = np.mean(training_data['observations'], axis=0, keepdims=True)
    std_obs = np.std(training_data['observations'], ddof=1)
    obs2states = lambda X: jp_model.phi(X) # centering and scaling are done inside phi

    # also retrieve global point of view parameter
    if env_args['auto_global_pov']:
        env_args['global_pov'] = True
    
    if env_args['global_pov']!=True:
        global_pov = (0.0,env_args['global_pov'],0.0) if env_args['global_pov'] and env_args['global_pov']!=0.0 else None
    else:
        global_pov = True

    # finally the env itself
    env = run_round_bot.make_round_bot_env(
            world={'name':env_args['world_name'], 'size':env_args['world_size']},
            texture=env_args['texture'],
            obssize=env_args['obssize'],
            winsize=[200,200] if args.visible_train else None,
            controller=controller,
            global_pov=global_pov,
            visible=args.visible_train,
            perspective = not env_args['orthogonal'],
            multiview=env_args['multiview'],
            focal=env_args['focal'],
            max_step=n_test_steps,
            random_start= True,
            distractors = env_args['distractors'],
            observation_transformation = obs2states,
            )


    ### learning iteration loop

    #    record loss for plotting
    total_loss_record = np.zeros(args.num_epochs)
    temporalcoherence_loss_record = np.zeros(args.num_epochs)
    proportionality_loss_record = np.zeros(args.num_epochs)
    causality_loss_record = np.zeros(args.num_epochs)
    repeatability_loss_record = np.zeros(args.num_epochs)

    for learning_epoch in range(args.num_epochs):

        states, history = tools.learn_states(
                            training_data=training_data,
                            model = jp_model,                                      
                            num_epochs=1,
                            batch_size = args.batch_size,
                            recordto='',
                            display=args.display,
                        ) 

        # record the loss history
        total_loss_record[learning_epoch] = history.history['loss'][0]
        temporalcoherence_loss_record[learning_epoch] = history.history['weighted_temporalcoherence_loss'][0]
        proportionality_loss_record[learning_epoch] = history.history['weighted_proportionality_loss'][0]
        causality_loss_record[learning_epoch] = history.history['weighted_causality_loss'][0]
        repeatability_loss_record[learning_epoch] = history.history['weighted_repeatability_loss'][0]

        # plot representations every learning
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


        # ### qlearning loop
        # for q_learning in range(n_qlearnings):

        #     # Perform fitted Q iterations and get states policy
        #     qit = fqi.Fitted_QIteration(n_rbf=n_rbf, n_actions=n_actions)
        #     # train a policy with q fitted iteration using an integer representation of the actions 'actions_int'
        #     state_policy = qit.fit_sk( states, actions_int, training_data['rewards'], 0.9, 5, recompute_mapping=True)
        #     # plug this policy into our policy class module
        #     state_policy = policy.Plug_policy(state_policy, env.controller.action_space_int)
        #     # create an agent with this policy
        #     rb_agent = agent.RoundBotAgent(state_policy)
        #     # test the agent in the env
        #     stats = rb_agent.run_in_env(env, n_test_episodes, seed=None) # 20 episodes as in the paper
        #     # record all episodes rewards for this evaluation run and this model and this learning
        #     test_performance[learning_epoch, q_learning, :,:] = np.reshape(np.array(stats['rewards']),[n_test_episodes, n_test_steps])            

        #     if args.verbose:
        #         print( 'Q fitted iteration test number : '+ str(q_learning) +'. Mean reward over'+ str(n_test_episodes)+' episodes : ' +\
        #             str( np.mean( stats['reward_ep'].flatten() ) )+' \n' )


if args.display :
    ### plotting 

    causality_loss_record = np.array(causality_loss_record)
    repeatability_loss_record = np.array(repeatability_loss_record)
    proportionality_loss_record = np.array(proportionality_loss_record)
    temporalcoherence_loss_record = np.array(temporalcoherence_loss_record)

    # plot the stacked losses' histories
    fig=plt.figure('Loss')
    #plt.plot(np.array(total_loss_record)/4.0)    
    plt.plot(causality_loss_record)
    plt.plot(repeatability_loss_record + causality_loss_record)
    plt.plot(proportionality_loss_record + repeatability_loss_record + causality_loss_record)
    plt.plot(temporalcoherence_loss_record + proportionality_loss_record + repeatability_loss_record + causality_loss_record)
    #plt.plot(temporalcoherence_loss_record + proportionality_loss_record + repeatability_loss_record + causality_loss_record)
    plt.title('Model loss')
    plt.ylabel('Loss')
    plt.xlabel('Epoch')
    plt.legend(['causality_loss','repeatability_loss','proportionality_loss','temporalcoherence_loss','total_loss',], loc='upper left')
    plt.show()

    input('Press any key to exit plotting')






