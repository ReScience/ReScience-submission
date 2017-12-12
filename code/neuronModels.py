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
# Neuron model equations and parameters.
# ----------------------------------------------------------------------------

from brian2 import *

def LIFparams():
	tau_m   = 10.0*ms		# membrane time constant
	tau_ref = 2.0*ms		# absolute refractory period
	Cm      = 250.0*pF		# membrane capacity
	v_r     = -65.0*mV		# reset potential
	v_th    = -50.0*mV		# fixed firing threshold
	return tau_m, tau_ref, Cm, v_r, v_th

# Leaky integrate-and-fire model equations
# dv/dt: equation 1 from the article
# dI/dt: equation 2 from the article
LIFmodel = '''
	dv/dt = (-v + v_r)/tau_m + (I+Iext)/Cm : volt (unless refractory)
	dI/dt = -I/tau_syn : amp
	Iext : amp
	'''
# Reset condition
resetLIF = '''
	v = v_r
	'''
