import numpy as np
import matplotlib.pyplot as plt
from CoreFunctions import *

### First we do the simulations for the Figure 4A
# The best action has a probability for a positive outcome of 0.75, and the other action of 0.25.
ProbabilityOfPositiveOutcome = np.array([.25, .75]) 

# The task lasts for 800 trials and there are 5,000 iterations (or runs) of 800 trials each.
NumberOfTrials = 800
NumberOfIterations = 5000

# The temperature parameter beta is fixed at 0.3 for all models.
beta = 0.3

# For the optimistic learner, the positive and negative learning rates are fixed at respectively 0.4 and 0.1
# For the rational learner, both positive and negative learning rates are fixed at 0.1
# For the pessimistic learner, the positive and negative learning rates are fixed at respectively 0.1 and 0.4
Alphas = np.matrix([[.4, .1], [.1, .1], [.1, .4]])
# For the meta-learner, the parameter w set at 0.1
w = .1

# modelIdx refers to the identity of the learner.
for modelIdx in range(4):

    # We will compute the performance on each trial for all the iterations and models
    performance = np.nan * np.ones((NumberOfTrials, NumberOfIterations))      

    for iterationIdx in range(NumberOfIterations):

        # The initial Q-values are fixed at 0. (The first element is the Q-value associated with the worst action, and the second with the best action.)
        Q = np.zeros(2)
        # We need to keep track of the reward history for the meta-learner
        Rewards = np.nan * np.ones(NumberOfTrials)

        # t is the trial index
        for t in range(NumberOfTrials):

            choice = randomChoice(Q, beta)
            performance[t, iterationIdx] = choice
            Rewards[t] = randomReward(ProbabilityOfPositiveOutcome[choice])

            if modelIdx < 3: # The update rule for the rational learners is the rule for the asymmetric model with PosAlpha = NegAlpha
                Q[choice] += updateAsymmetricLearner(Q[choice], Rewards[t], Alphas[modelIdx,0], Alphas[modelIdx,1])
            else: # The update rule for the meta-learner
                Q[choice] += updateMetaLearner(Q[choice], Rewards[t], w, t, Rewards[:t])

    if modelIdx == 0:
        meanPerformanceOptimistic  = np.mean(performance, axis=1)
    elif modelIdx == 1:
        meanPerformanceRational    = np.mean(performance, axis=1)
    elif modelIdx == 2:
        meanPerformancePessimistic = np.mean(performance, axis=1)
    else:
        meanPerformanceMetaLearner = np.mean(performance, axis=1)

# We now draw the figure 4A:
plt.figure(4)

plt.subplot(3,2,1)
plt.plot(meanPerformanceOptimistic, 'g', linewidth=1)
plt.plot(meanPerformanceRational, 'b', linewidth=1)
plt.plot(meanPerformancePessimistic, 'r', linewidth=1) 
plt.plot(meanPerformanceMetaLearner, 'm', linewidth=1)

plt.axis([0, 800, .45, 1])

plt.yticks([.5, 1.0])
plt.ylabel('P(best choice)')

plt.xticks([100, 800])
plt.xlabel('Episodes')

# Here we insert the texts around the figure:
plt.text(100, 0.68,'p(reward|arm#1) = 0.25')
plt.text(100, 0.55,'p(reward|arm#2) = 0.75')
plt.text(0.01, 0.9, 'A', fontsize=14, weight='bold', transform=plt.gcf().transFigure)
plt.text(0.01, 0.37, 'B', fontsize=14, weight='bold', transform=plt.gcf().transFigure)
plt.text(0.62, 0.85, 'Meta-learner', fontsize=12, color = 'm', transform=plt.gcf().transFigure)
plt.text(0.62, 0.78, 'Optimistic', fontsize=12, color = 'g', transform=plt.gcf().transFigure)
plt.text(0.62, 0.71, 'Rational', fontsize=12, color = 'b', transform=plt.gcf().transFigure)
plt.text(0.62, 0.64, 'Pessimistic', fontsize=12, color = 'r', transform=plt.gcf().transFigure)



### Then we start the simulations for Figure 4B:
# Now there are 2 different tasks, in which the model can choose between three possible actions.
# In the "low-reward" task, the best action has a probability for a positive outcome of 0.2, the 'middle' action of 0.15, and the worst action of 0.1. 
# In the "high-reward" task, the best action has a probability for a positive outcome of 0.9, the 'middle' action of 0.85, and the worst action of 0.8. 
ProbabilityOfPositiveOutcome = np.matrix([[.1, .15, .2], [.8, .85, .9]])

for modelIdx in range(4):

    performance = np.nan * np.ones((NumberOfTrials, NumberOfIterations, len(ProbabilityOfPositiveOutcome)))

    # ConditionIdx refers to the task: 0 for the low-reward task, 1 for the high-reward task.
    for conditionIdx in range(len(ProbabilityOfPositiveOutcome)):

        for iterationIdx in range(NumberOfIterations):
            
            # The initial Q-values for the three possible actions are fixed at 0.
            # (The first element is the Q-value associated with the worst action, and the second with the middle action and the third with the best action.)
            Q = np.zeros(3)
            # Again we keep track of the reward history for the meta-learner
            Rewards = np.nan * np.ones(NumberOfTrials)
            
            for t in range(NumberOfTrials):
                
                choice = randomChoice(Q, beta)
                performance[t, iterationIdx, conditionIdx] = 1 * (choice == 2)
                Rewards[t] = randomReward(ProbabilityOfPositiveOutcome[conditionIdx, choice])

                if modelIdx < 3: # The update rule for the rational learners is the rule for the asymmetric model with PosAlpha = NegAlpha
                    Q[choice] += updateAsymmetricLearner(Q[choice], Rewards[t], Alphas[modelIdx,0], Alphas[modelIdx,1])
                else: # The update rule for the meta-learner
                    Q[choice] += updateMetaLearner(Q[choice], Rewards[t], w, t, Rewards[:t])
                               
    if modelIdx == 0:
        meanPerformanceOptimisticLow,  meanPerformanceOptimisticHigh  = np.mean(performance[:,:,0], axis=1), np.mean(performance[:,:,1], axis=1)
    elif modelIdx == 1:
        meanPerformanceRationalLow,    meanPerformanceRationalHigh    = np.mean(performance[:,:,0], axis=1), np.mean(performance[:,:,1], axis=1)
    elif modelIdx == 2:
        meanPerformancePessimisticLow, meanPerformancePessimisticHigh = np.mean(performance[:,:,0], axis=1), np.mean(performance[:,:,1], axis=1)
    else:
        meanPerformanceMetaLearnerLow, meanPerformanceMetaLearnerHigh = np.mean(performance[:,:,0], axis=1), np.mean(performance[:,:,1], axis=1)

# The left panel of figure 4B:
plt.subplot(3,2,5)
plt.plot(meanPerformanceOptimisticLow, 'g', linewidth=1)
plt.plot(meanPerformanceRationalLow, 'b', linewidth=1)
plt.plot(meanPerformancePessimisticLow, 'r', linewidth=1) 
plt.plot(meanPerformanceMetaLearnerLow, 'm', linewidth=1)

plt.title('Low reward')
plt.axis([0, 800, .3, .6])

plt.yticks([.3, .4, .5, .6])
plt.ylabel('P(best choice)')

plt.xticks([100, 800])
plt.xlabel('Episodes')

# The right panel of figure 4B:
plt.subplot(3,2,6)
plt.plot(meanPerformanceOptimisticHigh, 'g', linewidth=1)
plt.plot(meanPerformanceRationalHigh, 'b', linewidth=1)
plt.plot(meanPerformancePessimisticHigh, 'r', linewidth=1) 
plt.plot(meanPerformanceMetaLearnerHigh, 'm', linewidth=1)

plt.title('High reward')
plt.axis([0, 800, .3, .6])

plt.yticks([.3, .4, .5, .6], ('','','',''))
plt.xticks([100, 800])
plt.xlabel('Episodes')

plt.savefig('Figure4.png', dpi=150)
plt.show()
