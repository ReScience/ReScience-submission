import numpy as np
import matplotlib.pyplot as plt
from CoreFunctions import *

# There are four different probabilities p of getting a positive outcome (r = 1). The probability of getting a negative outcome (r = -1) is therefore 1-p.
ProbabilityOfPositiveOutcome = np.array([.1, .2, .8, .9]) 

# The task lasts for 800 trials for each probability and there are 5,000 iterations (or runs) of 800 trials each.
NumberOfTrials = 800
NumberOfIterations = 5000

# The positive and negative learning rates used will depend on the model identity:
# For the optimistic learner, the positive and negative learning rates are fixed at respectively 0.4 and 0.1
# For the rational learner, both positive and negative learning rates are fixed at 0.1
# For the pessimistic learner, the positive and negative learning rates are fixed at respectively 0.1 and 0.4
Alphas = np.matrix([[.4, .1], [.1, .1], [.1, .4]])

# modelIdx refers to the identity of the model: 0 for optimistic, 1 for rational and 2 for pessimistic
for modelIdx in range(3):

	# The initial Q-values are fixed at 0.
	Q = np.zeros((NumberOfIterations, len(ProbabilityOfPositiveOutcome)))

	# ConditionIdx refers to the different possible probabilities for positive outcomes (p = 0.1, 0.2, 0.8 or 0.9)
	for conditionIdx in range(len(ProbabilityOfPositiveOutcome)):

		for iterationIdx in range(NumberOfIterations):

			# t is the trial index
			for t in range(NumberOfTrials):
				reward = randomReward(ProbabilityOfPositiveOutcome[conditionIdx])
				Q[iterationIdx, conditionIdx] += updateAsymmetricLearner(Q[iterationIdx, conditionIdx], reward, Alphas[modelIdx,0], Alphas[modelIdx,1])

	if modelIdx == 0:
		meanQoptimistic,  varQoptimistic  = np.mean(Q, axis=0), np.var(Q,  axis=0)
	elif modelIdx == 1:
		meanQrational,    varQrational    = np.mean(Q, axis=0), np.var(Q,  axis=0)
	else:
		meanQpessimistic, varQpessimistic = np.mean(Q, axis=0), np.var(Q,  axis=0)

# We now draw the figure 1.
plt.figure(1)
plt.subplot(2,2,1)

TrueValuesOfQ = 2 * ProbabilityOfPositiveOutcome - 1
for i in range(len(ProbabilityOfPositiveOutcome)):
	plt.errorbar([0,1,2],[meanQpessimistic[i], meanQrational[i], meanQoptimistic[i]], 
             [varQpessimistic[i], varQrational[i], varQoptimistic[i]], marker='o', color = 'k', capsize = 4, antialiased=True)
	plt.plot([-.5, 2.5], [TrueValuesOfQ[i], TrueValuesOfQ[i]], 'k:', antialiased=True)

plt.axis([-.5, 2.5, -1, 1])

plt.locator_params(axis='y',nbins=5)
plt.ylabel(r'Q$\infty$',rotation=0)

plt.xticks([0, 1, 2], ('0.25', '1', '4'))
plt.xlabel(r'$\alpha^+/ \alpha^-$')

plt.savefig('Figure1.png', dpi=150)
plt.show()
