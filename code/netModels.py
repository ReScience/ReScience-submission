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

def PDnet(NeuronGroup, stim, bg_type, w_ex, g, bg_freq, nsyn_type, thal):

	w_ex = w_ex*pA		   	# excitatory synaptic weight
	std_w_ex = 0.1*w_ex     # standard deviation weigth

	# Background number per layer
	if bg_type == 0:
		# layer-specific:
		bg_layer = bg_layer_specific
	elif bg_type == 1:
		# layer-independent:
		bg_layer = bg_layer_independent
	elif bg_type == 2:
		bg_layer = zeros(8)
		#layer-independent-random:
		for i in range(0,8,2):
		    # range of the number of inputs given to an excitatory population:
		    exc_bound_A = bg_layer_specific[i]
		    exc_bound_B = bg_layer_independent[i]
		    diff_exc = abs(exc_bound_A-exc_bound_B)

		    # randomly choosing a number for the external input to an excitatory population:
		    if exc_bound_A<=exc_bound_B: exc_input = exc_bound_A + rand()*diff_exc
		    elif exc_bound_A>exc_bound_B: exc_input = exc_bound_B + rand()*diff_exc

		    # range of the number of inputs given to an inhibitory population:
		    if i!=6: inh_bound_A = ((1-0.1)/(1+0.1))*exc_input		# eq. 4 from the article
		    else: inh_bound_A = ((1-0.2)/(1+0.2))*exc_input		# eq. 4 from the article
		    inh_bound_B = exc_input
		    diff_inh = abs(inh_bound_A-inh_bound_B)

		    # randomly choosing a number for the external input to an inhibitory population:
		    if inh_bound_A<=inh_bound_B: inh_input = inh_bound_A + rand()*diff_inh
		    else: inh_input = inh_bound_B + rand()*diff_inh

		    # array created to save the values:
		    bg_layer[i] = int(exc_input)
		    bg_layer[i+1] = int(inh_input)

	pop = [] # Stores NeuronGroups, one for each population
	for r in range(0, 8):
		pop.append(NeuronGroup[nn_cum[r]:nn_cum[r+1]])

		# DC-current normalized by population
		if (stim == 1):
			NeuronGroup.Iext[pop[r]] = 0.3512*pA*bg_layer[r]

	###########################################################################
	# Creating synapse connections
	###########################################################################

	syn_model = '''
				w:amp			# synaptic weight
				'''

	# equations executed only when presynaptic spike occurs:
	# for excitatory connections
	pre_eq = '''
			I_post += w
			'''

	con = [] # Stores connections

	###########################################################################
	# Connecting neurons
	###########################################################################
	pre_index = []
	post_index = []

	for c in range(0, 8):
		for r in range(0, 8):

			if (nsyn_type==0):
				# number of synapses calculated with equation 3 from the article
				nsyn = int(log(1.0-table[r][c])/log(1.0 - (1.0/float(n_layer[c]*n_layer[r]))))
			elif (nsyn_type==1):
				# number of synapses calculated with equation 5 from the article
				nsyn = int(n_layer[c]*n_layer[r]*table[r][c])

			pre_index = randint(n_layer[c], size=nsyn)
			post_index = randint(n_layer[r], size=nsyn)

			if nsyn<1:
				pass
			else:
				# Excitatory connections
				if (c % 2) == 0:
					# Synaptic weight from L4e to L2/3e is doubled
					if c == 2 and r == 0:
						con.append(Synapses(pop[c], pop[r], model=syn_model, on_pre=pre_eq))
						con[-1].connect(i = pre_index, j = post_index)
						con[-1].w = '2.0*clip((w_ex + std_w_ex*randn()),w_ex*0.0, w_ex*inf)'
					else:
						con.append(Synapses(pop[c], pop[r], model=syn_model, on_pre=pre_eq))
						con[-1].connect(i = pre_index, j = post_index)
						con[-1].w = 'clip((w_ex + std_w_ex*randn()),w_ex*0.0, w_ex*inf)'
					con[-1].delay = 'clip(d_ex + std_d_ex*randn(), 0.1*ms, d_ex*inf)'

				# Inhibitory connections
				else:
					con.append(Synapses(pop[c], pop[r], model=syn_model, on_pre=pre_eq))
					con[-1].connect(i = pre_index, j = post_index)
					con[-1].w = '-g*clip((w_ex + std_w_ex*randn()),w_ex*0.0, w_ex*inf)'
					con[-1].delay = 'clip(d_in + std_d_in*randn(), 0.1*ms, d_in*inf)'

	###########################################################################
	# Creating poissonian background inputs
	###########################################################################
	bg_in  = []
	if (stim==0):
		for r in range(0, 8):
			bg_in.append(PoissonInput(pop[r], 'I', bg_layer[r], bg_freq*Hz, weight=w_ex))

	###############################################################################
	# Creating thalamic neurons as poissonian inputs
	###############################################################################
	thal_con = []
	thal_input = []
	if thal=="ON":
		thal_input = PoissonGroup(n_layer[8], rates=120.0*Hz)	#from PD paper: rates=15Hz
		for r in range(0,8):
			thal_con.append(Synapses(thal_input, pop[r], model=syn_model, on_pre=pre_eq))
			thal_con[-1].connect(p=table[r][8])
			thal_con[-1].w = 0.0

	###########################################################################
	# Creating spike monitors
	###########################################################################
	smon_net = SpikeMonitor(NeuronGroup)

	return pop, con, bg_in, smon_net, thal_input ,thal_con
