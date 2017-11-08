# ----------------------------------------------------------------------------
# Contributors: Renan O. Shimoura
#               Nilton L. Kamiji
#               Rodrigo F. O. Pena
#               Vinicius L. Cordeiro
#               Cesar C. Ceballos
#               Cecilia Romaro
#               Antonio C. Roque
# ----------------------------------------------------------------------------
# References:
#
# *The Cell-Type Specific Cortical Microcircuit: Relating Structure and Activity
# in a Full-Scale Spiking Network Model*,
# Tobias C. Potjans and Markus Diesmann,
# Cerebral Cortex, 24(3):785-806, 2014.
# ----------------------------------------------------------------------------
# File description:
#
# Function responsible to make connections between different populations of
# neurons and to connect stimulation (background or DC currents).
# ----------------------------------------------------------------------------

from netParams import *

def PDnet(NeuronGroup, stim, bg_type, w_ex, g):

	std_w_ex = 0.1*w_ex        # standard deviation weigth

	# Background number per layer
	if bg_type == 0:
		# layer-specific
		bg_layer = bg_layer_specific
	else:
		# layer-independent
		bg_layer = bg_layer_independent

	pop = [] # Stores NeuronGroups, one for each population
	for r in range(0, 8):
		pop.append(NeuronGroup[nn_cum[r]:nn_cum[r+1]])

		# DC-current normalized by population
		if (stim == 1 or stim == 2):
			NeuronGroup.Iext[pop[r]] = 0.3512*pA*bg_layer[r]

	###########################################################################
	# Creating synapse connections
	###########################################################################

	syn_model = '''
				w:amp			# synaptic weight
				'''

	# equations executed only when presynaptic spike occurs:
	# for excitatory connections
	pre_eq_exc1 = '''
			I_post += w
			'''
	# synaptic weight from L4e to L2/3e is doubled
	pre_eq_exc2 = '''
			I_post += 2*w
			'''
	# for inhibitory connections
	pre_eq_inh = '''
			I_post -= g*w
			'''

	con = [] # Stores connections

	###########################################################################
	# Connecting neurons
	###########################################################################
	pre_index = []
	post_index = []

	for c in range(0, 8):
		for r in range(0, 8):

			#nsyn = int(n_layer[c]*n_layer[r]*table[r][c])
			nsyn = int(log(1.0-table[r][c])/log(1.0 - (1.0/float(n_layer[c]*n_layer[r]))))
			pre_index = randint(n_layer[c], size=nsyn)
			post_index = randint(n_layer[r], size=nsyn)

			if nsyn<1:
				pass
			else:
				# Excitatory connections
				if (c % 2) == 0:
					# Synaptic weight from L4e to L2/3e is doubled
					if c == 2 and r == 0:
						con.append(Synapses(pop[c], pop[r], model=syn_model, on_pre=pre_eq_exc2))
						con[-1].connect(i = pre_index, j = post_index)
					else:
						con.append(Synapses(pop[c], pop[r], model=syn_model, on_pre=pre_eq_exc1))
						con[-1].connect(i = pre_index, j = post_index)
					con[-1].w = 'clip((w_ex + std_w_ex*randn()),w_ex*0.0, w_ex*inf)'
					con[-1].delay = 'clip(d_ex + std_d_ex*randn(), 0.1*ms, d_ex*inf)'

				# Inhibitory connections
				else:
					con.append(Synapses(pop[c], pop[r], model=syn_model, on_pre=pre_eq_inh))
					con[-1].connect(i = pre_index, j = post_index)
					con[-1].w = 'clip((w_ex + std_w_ex*randn()),w_ex*0.0, w_ex*inf)'
					con[-1].delay = 'clip(d_in + std_d_in*randn(), 0.1*ms, d_in*inf)'


	###########################################################################
	# Creating poissonian background inputs
	###########################################################################
	bg_in  = []
	if (stim==0 or stim==2):
		for r in range(0, 8):
			bg_in.append(PoissonInput(pop[r], 'I', bg_layer[r], 8*Hz, weight=w_ex))

	###########################################################################
	# Creating spike monitors
	###########################################################################
	smon_net = SpikeMonitor(NeuronGroup)

	return pop, con, bg_in, smon_net
