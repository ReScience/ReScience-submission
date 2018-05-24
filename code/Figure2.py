import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from CoreFunctions import *

# In the "low-reward" task, the best action has a probability for a positive outcome of 0.2, and the other action of 0.1. 
# In the "high-reward' task, these two probabilities are respectively 0.9 and 0.8.
ProbabilityOfPositiveOutcome = np.matrix([[.1, .2], [.8, .9]])

# Each task lasts for 800 trials and there are 5,000 iterations (or runs) of 800 trials each.
NumberOfTrials = 800
NumberOfIterations = 5000

# The temperature parameter beta is fixed at 0.3 for all three models (optimistic, rational and pessimistic).
beta = 0.3

# For the optimistic learner, the positive and negative learning rates are fixed at respectively 0.4 and 0.1
# For the rational learner, both positive and negative learning rates are fixed at 0.1
# For the pessimistic learner, the positive and negative learning rates are fixed at respectively 0.1 and 0.4
Alphas = np.matrix([[.4, .1], [.1, .1], [.1, .4]])

# modelIdx refers to the identity of the model: 0 for optimistic, 1 for rational and 2 for pessimistic
for modelIdx in range(3):

	# We will compute the performance on each trial for all the iterations and models, and the proportion of choice switch after 800 trials.
	performance = np.nan * np.ones((NumberOfTrials, NumberOfIterations, len(ProbabilityOfPositiveOutcome)))
	switch = np.nan * np.ones((NumberOfIterations, len(ProbabilityOfPositiveOutcome)))

	# ConditionIdx refers to the task: 0 for the low-reward task, 1 for the high-reward task.
	for conditionIdx in range(len(ProbabilityOfPositiveOutcome)):

		for iterationIdx in range(NumberOfIterations):
            
            # The initial Q-values are fixed at 0. (The first element is the Q-value associated with the worst action, and the second with the best action.)
			Q = np.zeros(2)        

			# t is the trial index
			for t in range(NumberOfTrials):
				choice = randomChoice(Q, beta)
				performance[t, iterationIdx, conditionIdx] = choice
				reward = randomReward(ProbabilityOfPositiveOutcome[conditionIdx, choice])
				Q[choice] += updateAsymmetricLearner(Q[choice], reward, Alphas[modelIdx,0], Alphas[modelIdx,1])

			# Finally we compute if the model will switch choices or not after 800 trials
			nextChoice = randomChoice(Q, beta)
			switch[iterationIdx, conditionIdx] = 1 * (nextChoice != choice)

	if modelIdx == 0:
		meanPerformanceOptimisticLow,  meanPerformanceOptimisticHigh  = np.mean(performance[:,:,0], axis=1), np.mean(performance[:,:,1], axis=1)
		switchOptimistic  = np.mean(switch, axis=0)
	elif modelIdx == 1:
		meanPerformanceRationalLow,    meanPerformanceRationalHigh    = np.mean(performance[:,:,0], axis=1), np.mean(performance[:,:,1], axis=1)
		switchRational    = np.mean(switch, axis=0)
	else:
		meanPerformancePessimisticLow, meanPerformancePessimisticHigh = np.mean(performance[:,:,0], axis=1), np.mean(performance[:,:,1], axis=1)
		switchPessimistic = np.mean(switch, axis=0)

plt.figure(2)

# The left panel of figure 2A:
plt.subplot(3,2,1)
plt.plot(meanPerformanceOptimisticLow, 'g', linewidth=1)
plt.plot(meanPerformanceRationalLow, 'b', linewidth=1)
plt.plot(meanPerformancePessimisticLow, 'r', linewidth=1)

plt.title('low reward')
plt.axis([0, 800, .45, .8])

plt.yticks([.5, .55, .6, .65, .7, .75, .8], ('0.5', '', '0.6', '', '0.7', '', '0.8'))
plt.ylabel('p(best choice)')

plt.xticks([0, 100, 200, 300, 400, 500, 600, 700, 800], ('0','','','','400','','','','800'))
plt.xlabel('trials')

# The right panel of figure 2A:
plt.subplot(3,2,2)
plt.plot(meanPerformanceOptimisticHigh, 'g', linewidth=1)
plt.plot(meanPerformanceRationalHigh, 'b', linewidth=1)
plt.plot(meanPerformancePessimisticHigh, 'r', linewidth=1)

plt.title('high reward')
plt.axis([0, 800, .45, .8])

cur_axes = plt.gca()
cur_axes.axes.get_yaxis().set_visible(False)

plt.xticks([0, 100, 200, 300, 400, 500, 600, 700, 800], ('0','','','','400','','','','800'))
plt.xlabel('trials')

# The figure 2B:
gs = gridspec.GridSpec(3, 4)
plt.subplot(gs[2, 1:-1])

plt.bar([1,6], switchOptimistic, width = 1, align = 'center', edgecolor = 'k', color='g')
plt.bar([2,7], switchRational, width = 1, align = 'center', edgecolor = 'k', color='b')
plt.bar([3,8], switchPessimistic, width = 1, align = 'center', edgecolor = 'k', color='r')

plt.axis([.5, 8.5, 0, 1])
plt.locator_params(axis='y',nbins=6)
plt.ylabel('p(switch)')
plt.xticks([])
plt.xlabel('low reward             high reward')

legendToPlot = 'ORP'
for indexLegend in range(3):
	plt.text(indexLegend+1, 0.7,legendToPlot[indexLegend], horizontalalignment='center', verticalalignment='center')
for indexLegend in range(3):
	plt.text(indexLegend+6, 0.7,legendToPlot[indexLegend], horizontalalignment='center', verticalalignment='center')

# Here we insert the text and legends around the figure:
plt.text(0.01, 0.9, 'A', fontsize=14, weight='bold', transform=plt.gcf().transFigure)
plt.text(0.01, 0.35, 'B', fontsize=14, weight='bold', transform=plt.gcf().transFigure)
for indexLegend in range(3):
	plt.text(.5, .85-.075*indexLegend, legendToPlot[indexLegend], fontsize=12, transform=plt.gcf().transFigure, horizontalalignment='center')
for indexLegend in range(3):
	plt.text(.93, .72+.06*indexLegend, legendToPlot[indexLegend], fontsize=12, transform=plt.gcf().transFigure, horizontalalignment='center')
    
plt.savefig('Figure2.png', dpi=150)
plt.show()
