# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2014, Nicolas P. Rougier
# Distributed under the (new) BSD License.
#
# Contributors: Nicolas P. Rougier (Nicolas.Rougier@inria.fr)
# -----------------------------------------------------------------------------
# References:
#
# * Interaction between cognitive and motor cortico-basal ganglia loops during
#   decision making: a computational study. M. Guthrie, A. Leblois, A. Garenne,
#   and T. Boraud. Journal of Neurophysiology, 109:3025–3040, 2013.
# -----------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from model import *

# --- Parameter
ms         = 0.001
settling   = 500*ms
trial      = 2500*ms
dt         = 1*ms

debug      = False
threshold  = 40
alpha_c    = 0.025
alpha_LTP  = 0.004
alpha_LTD  = 0.002
Wmin, Wmax = 0.25, 0.75
tau        = 0.01
clamp      = Clamp(min=0, max=1000)
sigmoid    = Sigmoid(Vmin=0, Vmax=20, Vh=16, Vc=3)

CTX = AssociativeStructure(
                 tau=tau, rest=- 3.0, noise=0.01, activation=clamp )
STR = AssociativeStructure(
                 tau=tau, rest=  0.0, noise=0.01, activation=sigmoid )
STN = Structure( tau=tau, rest=-10.0, noise=0.01, activation=clamp )
GPI = Structure( tau=tau, rest=+10.0, noise=0.03, activation=clamp )
THL = Structure( tau=tau, rest=-40.0, noise=0.01, activation=clamp )
structures = (CTX, STR, STN, GPI, THL)

def weights(shape):
    Wmin, Wmax = 0.25, 0.75
    N = np.random.normal(0.5, 0.005, shape)
    N = np.minimum(np.maximum(N, 0.0),1.0)
    return (Wmin+(Wmax-Wmin)*N)


# These weights will change (learning)
W = weights(4)

connections = [
    OneToOne( CTX.cog.V, STR.cog.Isyn, W,            gain=+1.0 ),
    OneToOne( CTX.mot.V, STR.mot.Isyn, weights(4),   gain=+1.0 ),
    OneToOne( CTX.ass.V, STR.ass.Isyn, weights(4*4), gain=+1.0 ),
    CogToAss( CTX.cog.V, STR.ass.Isyn, weights(4),   gain=+0.2 ),
    MotToAss( CTX.mot.V, STR.ass.Isyn, weights(4),   gain=+0.2 ),
    OneToOne( CTX.cog.V, STN.cog.Isyn, np.ones(4),   gain=+1.0 ),
    OneToOne( CTX.mot.V, STN.mot.Isyn, np.ones(4),   gain=+1.0 ),
    OneToOne( STR.cog.V, GPI.cog.Isyn, np.ones(4),   gain=-2.0 ),
    OneToOne( STR.mot.V, GPI.mot.Isyn, np.ones(4),   gain=-2.0 ),
    AssToCog( STR.ass.V, GPI.cog.Isyn, np.ones(4),   gain=-2.0 ),
    AssToMot( STR.ass.V, GPI.mot.Isyn, np.ones(4),   gain=-2.0 ),
    OneToAll( STN.cog.V, GPI.cog.Isyn, np.ones(4),   gain=+1.0 ),
    OneToAll( STN.mot.V, GPI.mot.Isyn, np.ones(4),   gain=+1.0 ),
    OneToOne( GPI.cog.V, THL.cog.Isyn, np.ones(4),   gain=-0.5 ),
    OneToOne( GPI.mot.V, THL.mot.Isyn, np.ones(4),   gain=-0.5 ),
    OneToOne( THL.cog.V, CTX.cog.Isyn, np.ones(4),   gain=+1.0 ),
    OneToOne( THL.mot.V, CTX.mot.Isyn, np.ones(4),   gain=+1.0 ),
    OneToOne( CTX.cog.V, THL.cog.Isyn, np.ones(4),   gain=+0.4 ),
    OneToOne( CTX.mot.V, THL.mot.Isyn, np.ones(4),   gain=+0.4 ),
]

cues_value = np.ones(4) * 0.5
cues_reward = np.array([3.0,2.0,1.0,0.0])/3.0

def set_trial(cues_cog, cues_mot):
    c1,c2 = cues_cog
    m1,m2 = cues_mot
    v = 7
    noise = 0.0001
    CTX.mot.Iext = 0
    CTX.cog.Iext = 0
    CTX.ass.Iext = 0
    CTX.mot.Iext[m1]  = v + np.random.normal(0,v*noise)
    CTX.mot.Iext[m2]  = v + np.random.normal(0,v*noise)
    CTX.cog.Iext[c1]  = v + np.random.normal(0,v*noise)
    CTX.cog.Iext[c2]  = v + np.random.normal(0,v*noise)
    CTX.ass.Iext[c1*4+m1] = v + np.random.normal(0,v*noise)
    CTX.ass.Iext[c2*4+m2] = v + np.random.normal(0,v*noise)


def clip(V, Vmin, Vmax):
    return np.minimum(np.maximum(V, Vmin), Vmax)


def iterate(dt):
    global connections, structures

    # Flush connections
    for connection in connections:
        connection.flush()

    # Propagate activities
    for connection in connections:
        connection.propagate()

    # Compute new activities
    for structure in structures:
        structure.evaluate(dt)


def reset():
    global cues_values, structures
    for structure in structures:
        structure.reset()


def learn(time, cues_cog, cues_mot, debug=True):
    # A motor decision has been made
    c1, c2 = cues_cog
    m1, m2 = cues_mot
    mot_choice = np.argmax(CTX.mot.V)
    cog_choice = np.argmax(CTX.cog.V)

    # The actual cognitive choice may differ from the cognitive choice
    # Only the motor decision can designate the chosen cue
    if mot_choice == m1:
        choice = c1
    else:
        choice = c2

    if choice == min(c1,c2):
        P.append(1)
    else:
        P.append(0)

    # Compute reward
    reward = np.random.uniform(0,1) < cues_reward[choice]

    # Compute prediction error
    #error = cues_reward[choice] - cues_value[choice]
    error = reward - cues_value[choice]

    # Update cues values
    cues_value[choice] += error* alpha_c

    # Learn
    lrate = alpha_LTP if error > 0 else alpha_LTD
    dw = error * lrate * STR.cog.V[choice]
    W[choice] = W[choice] + dw * (W[choice]-Wmin)*(Wmax-W[choice])

    if not debug: return

    # Just for displaying ordered cue
    oc1,oc2 = min(c1,c2), max(c1,c2)
    if choice == oc1:
        print "Choice:          [%d] / %d  (good)" % (oc1,oc2)
    else:
        print "Choice:           %d / [%d] (bad)" % (oc1,oc2)
    print "Reward (%3d%%) :   %d" % (int(100*cues_reward[choice]),reward)
    print "Mean performance: %.3f" % np.array(P).mean()
    print "Mean reward:      %.3f" % np.array(R).mean()
    print "Response time:    %d ms" % (time)




if 0:
    R = np.load("250-simulations.npy")

else:
    R = np.zeros((250,120))
    # All combinations of cues or positions
    Z = [[0,1], [0,2], [0,3], [1,2], [1,3], [2,3]]
    # 20 x all cues combinations
    C = np.repeat(np.arange(6),20)
    # 20 x all cues positions
    M = np.repeat(np.arange(6),20)

    # 250 experiments
    for k in range(250):

        np.random.shuffle(C)
        np.random.shuffle(M)
        cues_value = np.ones(4) * 0.5
        W[...] = weights(4)
        P = []

        # 120 trials
        for j in range(120):
            reset()

            cues_cog = Z[C[j]]
            cues_mot = Z[M[j]]

            # Settling phase (500ms)
            i0 = 0
            i1 = i0+int(settling/dt)
            for i in xrange(i0,i1):
                iterate(dt)

            # Trial setup
            set_trial(cues_cog, cues_mot)

            # Learning phase (2500ms)
            i0 = int(settling/dt)
            i1 = i0+int(trial/dt)
            decision = False
            for i in xrange(i0,i1):
                iterate(dt)
                # Test if a decision has been made
                if CTX.mot.delta > threshold:
                    decision = True
                    learn(time=i-500, cues_cog=cues_cog, cues_mot=cues_mot, debug=debug)
                    break
            if not decision:
                P.append(0)

        R[k] = P
        print "Experiment %d: %.3f" % (k+1,np.array(P).mean())
    np.save("250-simulations.npy", R)


# Display
P = np.zeros((250,120))
for i in range(250):
    for j in range(120):
        P[i,j] = np.mean(R[i,j])

plt.figure(figsize=(12,5.5))
ax = plt.subplot(111)
ax.patch.set_facecolor("w")
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.yaxis.set_ticks_position('left')
ax.yaxis.set_tick_params(direction="in")
ax.xaxis.set_ticks_position('bottom')
ax.xaxis.set_tick_params(direction="in")

X = 1+np.arange(120)
plt.plot(X, P.mean(axis=0), c='b', lw=2)
plt.plot(X, P.mean(axis=0)+P.var(axis=0), c='b',lw=.5)
plt.plot(X, P.mean(axis=0)-P.var(axis=0), c='b',lw=.5)
plt.fill_between(X, P.mean(axis=0)+P.var(axis=0),
                    P.mean(axis=0)-P.var(axis=0), color='b', alpha=.1)
plt.xlabel("Trial number", fontsize=16)
plt.ylabel("Performance", fontsize=16)
plt.ylim(0,1.0)
plt.xlim(1,120)
plt.savefig("figure-2.pdf")
plt.show()

