import numpy as np
import matplotlib.pyplot as plt
from CoreFunctions import *

# In the "low-reward" task, the best action has a probability for a positive outcome of 0.2, and the other action of 0.1. 
# In the "high-reward' task, these two probabilities are respectively 0.9 and 0.8.
ProbabilityOfPositiveOutcome = np.matrix([[.1, .2], [.8, .9]])

# Each task lasts for 800 trials and there are 5,000 iterations (or runs) of 800 trials each.
NumberOfTrials = 800
NumberOfIterations = 5000

# The temperature parameter beta is fixed at 0.3 for all models (the rational learners and the meta-learner).
beta = 0.3

# There are three rational learners, each with a different learning rates (0.01, 0.1 or 0.4)
Alphas = np.array([.01, .1, .4])
# And a meta-learner with a parameter w set at 0.1 (w replaces the standard learning rate parameter)
w = .1

# modelIdx refers to the identity of the model.
for modelIdx in range(4):

    # We will compute the performance on each trial for all the iterations and models
    performance = np.nan * np.ones((NumberOfTrials, NumberOfIterations, len(ProbabilityOfPositiveOutcome)))

    # ConditionIdx refers to the task: 0 for the low-reward task, 1 for the high-reward task.
    for conditionIdx in range(len(ProbabilityOfPositiveOutcome)):      

        for iterationIdx in range(NumberOfIterations):
            
            # The initial Q-values are fixed at 0. (The first element is the Q-value associated with the worst action, and the second with the best action.)
            Q = np.zeros(2)
            # We now need to keep track of the reward history for the meta-learner
            Rewards = np.nan * np.ones(NumberOfTrials)
            
            # t is the trial index
            for t in range(NumberOfTrials):

                choice = randomChoice(Q, beta)
                performance[t, iterationIdx, conditionIdx] = choice
                Rewards[t] = randomReward(ProbabilityOfPositiveOutcome[conditionIdx, choice])

                if modelIdx < 3: # The update rule for the rational learners is the asymmetric model rule with PosAlpha = NegAlpha
                    Q[choice] += updateAsymmetricLearner(Q[choice], Rewards[t], Alphas[modelIdx], Alphas[modelIdx])
                else: # The update rule for the meta-learner
                    Q[choice] += updateMetaLearner(Q[choice], Rewards[t], w, t, Rewards[:t])

    if modelIdx == 0:
        meanPerformanceRational1Low,   meanPerformanceRational1High   = np.mean(performance[:,:,0], axis=1), np.mean(performance[:,:,1], axis=1)
    elif modelIdx == 1:
        meanPerformanceRational2Low,   meanPerformanceRational2High   = np.mean(performance[:,:,0], axis=1), np.mean(performance[:,:,1], axis=1)
    elif modelIdx == 2:
        meanPerformanceRational3Low,   meanPerformanceRational3High   = np.mean(performance[:,:,0], axis=1), np.mean(performance[:,:,1], axis=1)
    else:
        meanPerformanceMetaLearnerLow, meanPerformanceMetaLearnerHigh = np.mean(performance[:,:,0], axis=1), np.mean(performance[:,:,1], axis=1)

plt.figure(3)

# The left panel of figure 3:
plt.subplot(3,2,1)
plt.plot(meanPerformanceRational1Low, '#008080', linewidth=1) # hex triplet for teal color
plt.plot(meanPerformanceRational2Low, '#1E32FF', linewidth=1) # hex triplet for blue color
plt.plot(meanPerformanceRational3Low, '#191970', linewidth=1) # hex triplet for dark violet color
plt.plot(meanPerformanceMetaLearnerLow, 'm', linewidth=1)

plt.title('low reward')
plt.axis([0, 800, .45, .8])

plt.yticks([.5, .55, .6, .65, .7, .75, .8], ('0.5', '', '0.6', '', '0.7', '', '0.8'))
plt.ylabel('p(best choice)')

plt.xticks([0, 100, 200, 300, 400, 500, 600, 700, 800], ('0','','','','400','','','','800'))
plt.xlabel('trials')


# The right panel of figure 3:
plt.subplot(3,2,2)
plt.plot(meanPerformanceRational1High, '#008080', linewidth=1) # hex triplet for teal color
plt.plot(meanPerformanceRational2High, '#1E32FF', linewidth=1) # hex triplet for blue color
plt.plot(meanPerformanceRational3High, '#191970', linewidth=1) # hex triplet for dark violet color
plt.plot(meanPerformanceMetaLearnerHigh, 'm', linewidth=1)

plt.title('high reward')
plt.axis([0, 800, .45, .8])

cur_axes = plt.gca()
cur_axes.axes.get_yaxis().set_visible(False)

plt.xticks([0, 100, 200, 300, 400, 500, 600, 700, 800], ('0','','','','400','','','','800'))
plt.xlabel('trials')

# The text around the figure:
legendToPlot = 'NR'
for indexLegend in range(2):
    plt.text(.5 , .855 - .11 * indexLegend, legendToPlot[indexLegend], fontsize=12, transform=plt.gcf().transFigure, horizontalalignment='center')
    plt.text(.93, .855 - .05 * indexLegend, legendToPlot[indexLegend], fontsize=12, transform=plt.gcf().transFigure, horizontalalignment='center')
    
plt.savefig('Figure3.png', dpi=150)
plt.show()
