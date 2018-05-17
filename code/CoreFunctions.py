import numpy as np

def randomReward(PosProba):
    """This function returns a random reward, either positive (r = 1) or negative (r = -1), depending on the probability of positive outcome parameter."""
    return 2 * np.random.choice(2, p = np.array([1 - PosProba, PosProba])) - 1

def randomChoice(Q, beta):
    """This function returns a random choice, depending on the Q-values and the beta parameter.
    The probabilities of choosing the different actions are computed from the Q-values throught a softmax or normalized exponential function."""
    return np.random.choice(len(Q), p = np.exp(Q/beta) / np.sum(np.exp(Q/beta)))

def updateAsymmetricLearner(Q, r, alphaPos, alphaNeg):
    """This function returns the Q-value update that needs to be done on each trial according to the asymmetric learner.
    Note that the rational (or symmetric) learner is implemented if alphaPos and alphaNeg are set to be equal.
    This update depends on the prior Q-value (Q), the reward received on this trial (r = +1 or -1), and the positive and negative learning rates (between 0 and 1)."""
    return alphaPos * (r - Q) * (r == +1) + alphaNeg * (r - Q) * (r == -1)

def updateMetaLearner(Q, r, w, t, h):
    """This function returns the Q-value update that needs to be done on each trial according to the meata-learner.
    This update depends on the prior Q-value (Q), the reward received on this trial (r = +1 or -1), the parameter w, the trial index (t), and the reward history (h)."""
    if t == 0:
        ProbaReward = .5
    else:
        ProbaReward = (np.mean(h) + 1)/2
    return w * (1 - ProbaReward) * (r - Q) * (r == +1) + w * ProbaReward * (r - Q) * (r == -1)
