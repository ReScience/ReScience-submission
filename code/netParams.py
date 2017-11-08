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
# Networks structure parameters.
# ----------------------------------------------------------------------------

from brian2 import *

###############################################################################
# Network parameters
###############################################################################
# Population size per layer
#          2/3e   2/3i   4e    4i    5e    5i    6e     6i    Th
n_layer = [20683, 5834, 21915, 5479, 4850, 1065, 14395, 2948, 902]

# Total cortical Population
N = sum(n_layer[:-1])

# Number of neurons accumulated
nn_cum = [0]
nn_cum.extend(cumsum(n_layer))

# Prob. connection table
table = array([[0.101,  0.169, 0.044, 0.082, 0.032, 0.,     0.008, 0.,     0.    ],
               [0.135,  0.137, 0.032, 0.052, 0.075, 0.,     0.004, 0.,     0.    ],
               [0.008,  0.006, 0.050, 0.135, 0.007, 0.0003, 0.045, 0.,     0.0983],
               [0.069,  0.003, 0.079, 0.160, 0.003, 0.,     0.106, 0.,     0.0619],
               [0.100,  0.062, 0.051, 0.006, 0.083, 0.373,  0.020, 0.,     0.    ],
               [0.055,  0.027, 0.026, 0.002, 0.060, 0.316,  0.009, 0.,     0.    ],
               [0.016,  0.007, 0.021, 0.017, 0.057, 0.020,  0.040, 0.225,  0.0512],
               [0.036,  0.001, 0.003, 0.001, 0.028, 0.008,  0.066, 0.144,  0.0196]])

# Synapses parameters

d_ex = 1.5*ms      	# Excitatory delay
std_d_ex = 0.75*ms 	# Std. Excitatory delay
d_in = 0.80*ms      # Inhibitory delay
std_d_in = 0.4*ms  	# Std. Inhibitory delay
tau_syn = 0.5*ms    # Post-synaptic current time constant

# Layer-specific background input
bg_layer_specific = array([1600, 1500 ,2100, 1900, 2000, 1900, 2900, 2100])

# Layer-independent background input
bg_layer_independent = array([2000, 1850 ,2000, 1850, 2000, 1850, 2000, 1850])
