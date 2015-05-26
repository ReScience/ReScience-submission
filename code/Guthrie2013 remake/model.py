# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2015, Meropi Topalidou
# Distributed under the (new) BSD License.
#
# Contributors: Meropi Topalidou (Meropi.Topalidou@inria.fr)
#				Nicolas Rougier  (Nicolas.Rougier@inria.fr)
# -----------------------------------------------------------------------------
# References:
#
# * Interaction between cognitive and motor cortico-basal ganglia loops during
#   decision making: a computational study. M. Guthrie, A. Leblois, A. Garenne,
#   and T. Boraud. Journal of Neurophysiology, 109:3025–3040, 2013.
# -----------------------------------------------------------------------------
from c_dana import *
from parameters import *

clamp   = Clamp(min=0, max=1000)
sigmoid = Sigmoid(Vmin=Vmin, Vmax=Vmax, Vh=Vh, Vc=Vc)

CTX = AssociativeStructure(
                 tau=tau, rest=CTX_rest, noise=Cortex_N, activation=clamp )
STR = AssociativeStructure(
                 tau=tau, rest=STR_rest, noise=Striatum_N, activation=sigmoid )
STN = Structure( tau=tau, rest=STN_rest, noise=STN_N, activation=clamp )
GPI = Structure( tau=tau, rest=GPI_rest, noise=GPi_N, activation=clamp )
THL = Structure( tau=tau, rest=THL_rest, noise=Thalamus_N, activation=clamp )

structures = (CTX, STR, STN, GPI, THL)

CUE = np.zeros(4, dtype=[("mot", float),
                         ("cog", float),
                         ("value" , float),
                         ("reward", float)])

choices  = np.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
cues_cog = choices.repeat(n_trials/6, axis = 0)
np.random.shuffle(cues_cog)


CUE["mot"]    = 0,1,2,3
CUE["cog"]    = 0,1,2,3
CUE["value"]  = 0.5
#CUE["reward"] = rewards

def weights(shape):
    N = np.random.normal(0.5, 0.005, shape)
    N = np.minimum(np.maximum(N, 0.0),1.0)
    return (Wmin+(Wmax-Wmin)*N)

W1 = (2*np.eye(4) - np.ones((4,4))).ravel()
W2 = (2*np.eye(16) - np.ones((16,16))).ravel()

connections = {
    "CTX.cog -> STR.cog" : OneToOne( CTX.cog.V, STR.cog.Isyn, weights(4)  ), # plastic (RL)
    "CTX.mot -> STR.mot" : OneToOne( CTX.mot.V, STR.mot.Isyn, weights(4)  ),
    "CTX.ass -> STR.ass" : OneToOne( CTX.ass.V, STR.ass.Isyn, weights(4*4)),
    "CTX.cog -> STR.ass" : CogToAss( CTX.cog.V, STR.ass.Isyn, weights(4)  ),
    "CTX.mot -> STR.ass" : MotToAss( CTX.mot.V, STR.ass.Isyn, weights(4)  ),
    "CTX.cog -> STN.cog" : OneToOne( CTX.cog.V, STN.cog.Isyn, np.ones(4)  ),
    "CTX.mot -> STN.mot" : OneToOne( CTX.mot.V, STN.mot.Isyn, np.ones(4)  ),
    "STR.cog -> GPI.cog" : OneToOne( STR.cog.V, GPI.cog.Isyn, np.ones(4)  ),
    "STR.mot -> GPI.mot" : OneToOne( STR.mot.V, GPI.mot.Isyn, np.ones(4)  ),
    "STR.ass -> GPI.cog" : AssToCog( STR.ass.V, GPI.cog.Isyn, np.ones(4)  ),
    "STR.ass -> GPI.mot" : AssToMot( STR.ass.V, GPI.mot.Isyn, np.ones(4)  ),
    "STN.cog -> GPI.cog" : OneToAll( STN.cog.V, GPI.cog.Isyn, np.ones(4)  ),
    "STN.mot -> GPI.mot" : OneToAll( STN.mot.V, GPI.mot.Isyn, np.ones(4)  ),
    "THL.cog -> CTX.cog" : OneToOne( THL.cog.V, CTX.cog.Isyn, np.ones(4)  ),  # changed
    "THL.mot -> CTX.mot" : OneToOne( THL.mot.V, CTX.mot.Isyn, np.ones(4)  ),  # changed
    "CTX.cog -> THL.cog" : OneToOne( CTX.cog.V, THL.cog.Isyn, np.ones(4)  ),  # changed
    "CTX.mot -> THL.mot" : OneToOne( CTX.mot.V, THL.mot.Isyn, np.ones(4)  ),  # changed
    "GPI.cog -> THL.cog" : OneToOne( GPI.cog.V, THL.cog.Isyn, np.ones(4), ), # changed
    "GPI.mot -> THL.mot" : OneToOne( GPI.mot.V, THL.mot.Isyn, np.ones(4), ), # changed
}
for name,gain in gains.items():
    connections[name].gain = gain

# -----------------------------------------------------------------------------
def set_trial(n=2, cog_shuffle=True, mot_shuffle=True, noise=noise, trial = 0, protocol = 'Guthrie', familiar = 'True'):

    if protocol == 'Guthrie':
		temp = cues_cog[trial,:]
		np.random.shuffle(temp)
		CUE["cog"][0], CUE["cog"][1] = temp[0], temp[1]
		if cog_shuffle:
			np.random.shuffle(CUE["cog"][:n])
    elif protocol == 'Piron':#if protocol == 'Piron':
		if familiar:
			CUE["cog"][0], CUE["cog"][1] = 0, 1
		else:

			CUE["cog"][0], CUE["cog"][1] = 2, 3
		if cog_shuffle:
			np.random.shuffle(CUE["cog"][:n])
    if mot_shuffle:
		np.random.shuffle(CUE["mot"])

    CTX.mot.Iext = 0
    CTX.cog.Iext = 0
    CTX.ass.Iext = 0
    for i in range(n):
        c, m = CUE["cog"][i], CUE["mot"][i]

        CTX.mot.Iext[m]     = V_cue #+ np.random.uniform(-noise/2,noise/2)
        CTX.cog.Iext[c]     = V_cue #+ np.random.uniform(-noise/2,noise/2)
        CTX.ass.Iext[c*4+m] = V_cue #+ np.random.uniform(-noise/2,noise/2)


def iterate(dt):
    # Flush connections
    for connection in connections.values():
        connection.flush()

    # Propagate activities
    for connection in connections.values():
        connection.propagate()

    # Compute new activities
    for structure in structures:
        structure.evaluate(dt)


def reset(protocol = 'Guthrie', W_STR = None):
    CUE["mot"]    = 0,1,2,3
    CUE["cog"]    = 0,1,2,3
    CUE["value"]  = 0.5
    if protocol == 'Guthrie':
    	CUE["reward"] = rewards_Guthrie
    elif protocol == 'Piron':
    	CUE["reward"] = rewards_Piron
    np.random.shuffle(cues_cog)
    np.random.shuffle(CUE["mot"][:n])
    # CUE["reward"] = rewards
    if W_STR is not None:
		connections["CTX.cog -> STR.cog"].weights = W_STR
    else:
		connections["CTX.cog -> STR.cog"].weights = weights(4)
    reset_activities()

def reset_activities():
    for structure in structures:
        structure.reset()
def history():
	history = np.zeros(duration, dtype=dtype)
	history["CTX"]["mot"] = CTX.mot.history[:duration]
	history["CTX"]["cog"] = CTX.cog.history[:duration]
	history["CTX"]["ass"] = CTX.ass.history[:duration]
	history["STR"]["mot"] = STR.mot.history[:duration]
	history["STR"]["cog"] = STR.cog.history[:duration]
	history["STR"]["ass"] = STR.ass.history[:duration]
	history["STN"]["mot"] = STN.mot.history[:duration]
	history["STN"]["cog"] = STN.cog.history[:duration]
	history["GPI"]["mot"] = GPI.mot.history[:duration]
	history["GPI"]["cog"] = GPI.cog.history[:duration]
	history["THL"]["mot"] = THL.mot.history[:duration]
	history["THL"]["cog"] = THL.cog.history[:duration]
	return history
def reset_history():
	CTX.mot.history[:duration] = 0
	CTX.cog.history[:duration] = 0
	CTX.ass.history[:duration] = 0
	STR.mot.history[:duration] = 0
	STR.cog.history[:duration] = 0
	STR.ass.history[:duration] = 0
	STN.mot.history[:duration] = 0
	STN.cog.history[:duration] = 0
	GPI.mot.history[:duration] = 0
	GPI.cog.history[:duration] = 0
	THL.mot.history[:duration] = 0
	THL.cog.history[:duration] = 0

def process(n=2, learning=True, P = [], D = [], AP = np.zeros(n), RP = np.zeros(n), reverse = False, reverse_all = True):
    # A motor decision has been made
    # The actual cognitive choice may differ from the cognitive choice
    # Only the motor decision can designate the chosen cue
    mot_choice = np.argmax(CTX.mot.U)
    cog_choice = np.argmax(CTX.cog.U)

    # The actual cognitive choice may differ from the cognitive choice
    # Only the motor decision can designate the chosen cue
    for i in range(n):
        #print mot_choice, CUE["mot"][:n][i]
        if mot_choice == CUE["mot"][:n][i]:
            choice = int(CUE["cog"][:n][i])
    if not reverse:
		if choice == min(CUE["cog"][:n][0],CUE["cog"][:n][1]):
			P.append(1)
		else:
			P.append(0)
    else:
		if reverse_all:
			if choice == max(CUE["cog"][:n][0],CUE["cog"][:n][1]):
				P.append(1)
			else:
				P.append(0)
		else:
			c1, c2 = np.sort(CUE["cog"][:n])
			if np.array_equal([c1, c2], [1,2]):
				if choice == max(CUE["cog"][:n][0],CUE["cog"][:n][1]):
					P.append(1)
				else:
					P.append(0)
			elif choice == min(CUE["cog"][:n][0],CUE["cog"][:n][1]):
				P.append(1)
			else:
					P.append(0)

    D.append(0 if cog_choice == choice else 1)
    AP[choice] += 1
    if learning:
		# Compute reward
		reward = float(np.random.uniform(0,1) < CUE["reward"][choice])
		RP[choice] += reward

		# Compute prediction error
		error = reward - CUE["value"][choice]

		# Update cues values
		CUE["value"][choice] += error* alpha_CUE

        # Reinforcement learning
		lrate = alpha_LTP if error > 0 else alpha_LTD
		dw = error * lrate * STR.cog.V[choice]
		W = connections["CTX.cog -> STR.cog"].weights
		W[choice] = W[choice] + dw * (Wmax-W[choice]) * (W[choice]-Wmin)
		connections["CTX.cog -> STR.cog"].weights = W
    return choice


def debug_learning(Wstr, cues_value, f = None):
		print "Cues Values			: ", cues_value, '\n'
		print "Striatal Weights		: ", Wstr
		if f is not None:
			f.write("\nCues Values			: "+ str(cues_value))
			f.write("\nStriatal Weights		: "+ str(Wstr))

def debug(f = None, cgchoice = None, c1 = None, c2 = None, m1 = None, m2 = None, reverse = False, reverse_all = True, P = [], reward = [], RT = [], R = [], D = [], RP = None, AP = None, mBc = [], ABC = [], NoMove = []):

	if cgchoice is not None:
		print "Choice:         ",
		if cgchoice == c1:
			print " 	[%d]" % c1,
		else:
			print " 	%d" % c1,
		if cgchoice == c2:
			print " [%d]" % c2,
		else:
			print " %d" % c2,
		if not reverse:
			if cgchoice == np.minimum(c1,c2):
				print " (good)"
			else:
				print " (bad)"
		else:
			if reverse_all:
				if cgchoice == np.maximum(c1,c2):
					print " (good)"
				else:
					print " (bad)"
			else:
				temp1, temp2 = np.sort([c1,c2])
				if np.array_equal([temp1, temp2],[1,2]):
					if cgchoice == np.maximum(c1,c2):
						print " (good)"
					else:
						print " (bad)"
				elif cgchoice == np.minimum(c1,c2):
					print " (good)"

				else:
					print " (bad)"

	if m1 is not None:
		print "Positions:         	 %d   %d" %(m1, m2)
	if NoMove:
		print "Mean No move trials		: %.3f %%" % (len(NoMove)/float(n_trials))
	if P:
		print "Mean performance	 	: %.3f %%" % (np.array(P).mean()*100)
	if D:
		print "Mean Different choices 	 	: %.3f %%" % (np.array(D).mean()*100)
	if mBc:
		print "Motor decision before Cognitive	: %.3f %%" % (np.array(mBc).mean()*100)
	if ABC:
		print "Activity before Cues		: %.3f %%" % (np.array(ABC).mean()*100)


	if RT:
		print "Mean Response time		: %.3f ms" % (np.array(RT).mean())
	if reward:
		print "Reward	  		 	: %d" % (reward)
		print "Mean reward		 	: %.3f %%" % (np.array(R).mean()*100)
	if RP is not None:
		print "Reward Probabilities		: ", RP/AP*100
	if AP is not None:
		print "Number of Chosen		: ", AP # Number of chosen cues
	if f:
		if cgchoice is not None:
			f.writ( "Choice:         ")
			if cgchoice == c1:
				f.write(" 	[%d]" % c1)
			else:
				f.write(" 	%d" % c1)
			if cgchoice == c2:
				f.write(" [%d]" % c2)
			else:
				f.write(" %d" % c2)
			if not reverse:
				if cgchoice == np.minimum(c1,c2):
					f.write(" (good)")
				else:
					f.write(" (bad)")
			else:
				temp1, temp2 = np.sort([c1,c2])
				if np.array_equal([temp1, temp2],[1,2]):
					if cgchoice == np.maximum(c1,c2):
						f.write(" (good)")
					else:
						f.write(" (bad)")
				elif cgchoice == np.minimum(c1,c2):
					f.write(" (good)")

				else:
					f.write(" (bad)")
				if 0:
					if cgchoice == np.maximum(c1,c2):
						f.write(" (good)")
					else:
						f.write(" (bad)")
			if NoMove:
				f.write("\nMean No move trials		: %.3f %%" % (len(NoMove)/float(n_trials)))
			if P:
				f.write("\nMean performance	 	: %.3f %%" % (np.array(P).mean()*100))
			if D:
				f.write("\nMean Different choices 	 	: %.3f %%" % (np.array(D).mean()*100))
			if mBc:
				f.write("\nMotor decision before Cognitive	: %.3f %%" % (np.array(mBc).mean()*100))
			if ABC:
				f.write("\nActivity before Cues		: %.3f %%" % (np.array(ABC).mean()*100))


			if RT:
				f.write("\nMean Response time		: %.3f ms" % (np.array(RT).mean()))
			if reward:
				f.write("\nReward	  		 	: %d" % (reward))
				f.write("\nMean reward		 	: %.3f %%" % (np.array(R).mean()*100))
			if RP is not None:
				f.write("\nReward Probabilities		: "+ str(RP/AP*100))
			if AP is not None:
				f.write("\nNumber of Chosen		: "+ str(AP))


def debug_total(P, D, ABC, NoMove, AP, RP = None, CV = None, Wstr = None, f = None):

	print "Mean Performance		: " , (P.mean(axis=1)).mean(axis = 0)*100, '%'
	print "Trials with diff move		: " , (D.mean()*100), '%'
	print "Mean trials with activity \nbefore cue			: " , (ABC.mean()*100), '%'
	print "Mean Trials with no move	: " , (NoMove.mean()*100), '%'
	if RP is not None:
		print "Mean Reward Probabilities	:" + str((RP/AP*100).mean(axis = 0))
	if CV is not None:
		print "Mean Cues Values		:" + str(CV.mean(axis = 0))
		print 'Mean Striatal Weights		: ' + str(Wstr[:,-1].mean(axis = 0))
	if f is not None:
		f.write("\nMean Performance		: " + str((P.mean(axis=0)*100)) + '%')
		f.write("\nTrials with diff move		: " + str((D.mean()*100)) + '%')
		f.write("\nMean trials with activity \nbefore cue			: " + str((ABC.mean()*100)) + '%')
		f.write("\nMean Trials with no move	: " + str((NoMove.mean()*100)) + '%')
		if RP is not None:
			f.write("\nMean Reward Probabilities	:" + str((RP/AP*100).mean(axis = 0)))
		if CV is not None:
			f.write("\nMean Cues Values		:" + str(CV.mean(axis = 0)))
			f.write('\nMean Striatal Weights		: ' + str(Wstr[:,-1].mean(axis = 0)))
