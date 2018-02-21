# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Mario Senden mario.senden@maastrichtuniversity.nl
# -----------------------------------------------------------------------------
# References:
#
# Gancarz, G., Grossberg, S. "A Neural Model of the Saccade Generator in the Reticular Formation." 
# Neural Networks 11, no. 7-8 (October 1998): 1159-74. doi:10.1016/S0893-6080(98)00096-3.
# -----------------------------------------------------------------------------
# File description:
# 
# Sets up the SG model developed by Gancarz & Grossberg (1998) in PyNEST for subsequent simulation.
# -----------------------------------------------------------------------------

import nest
import numpy as np


###########################################
#### general setup & NEST initialization ##
###########################################

# simulation parameters 
dt     		=  0.05		# time step for numerical integration (in ms)
delay  		=  0.05		# conductance delay (in ms) - setting this equal to time step implies no delay
tau    		= 50.00		# time constant (in ms) 
sigma  		=  0.00		# scaling factor of additive noise


# NEST kernel initialization
nest.ResetKernel()	 
nest.SetKernelStatus({'resolution': dt, 'use_wfr': False})


###########################################
#### creation of neuron & device objects ##
###########################################

# neuron parameters
# alpha				: upper limit for a piecewise linear unit (2nd threshold for threshold_lin_rate)
# g_in				: gain in multiplicative coupling for inhibitory input
# lambda			: passive decay rate 
# linear_summation		: if false, gain function is applied to input before summation; if true, it is applied afterwards
# mult_coupling 		: if false, neuron does not exhibit multiplicative coupling; if true, it is does
# rectify_output		: if false, activations updated according to numerical integration; if true, activations are rectified after each update
# theta 			: lower limit for threshold (piecewise) linear unit (1st thresold for threshold_lin_rate)
# theta_ex			: offset in multiplicative coupling for excitatory input
# theta_in			: offset in multiplicative coupling for inhibitory input	

Params_llbn = {'tau': tau,'std': sigma,
			   'lambda':1.3,
			   'rectify_output': True,
			   'linear_summation': True}
Params_ebn  = {'tau': tau,'std': sigma,
			   'lambda':3.5,
			   'mult_coupling': True,
			   'theta_ex':2.,'theta_in':1.,
			   'rectify_output': True,
			   'linear_summation': True}
Params_ibn  = {'tau': tau,'std': sigma,
			   'lambda':2.4,
			   'rectify_output': True,
			   'linear_summation': True}
Params_opn  = {'tau': tau,'std': sigma,
			   'lambda':0.2,
			   'mult_coupling': True,
			   'theta_ex':1.,'theta_in':0.4,
			   'g_in':3.5, 
			   'rectify_output': True,
			   'linear_summation': True}
Params_tn   = {'tau': tau,'std': sigma,
			   'lambda':0.,'rate':.5, 
			   'rectify_output': False,
			   'linear_summation': True}
Params_sc 	= {'tau': tau,'std': sigma}
Params_bias = {'tau': dt, 'std': sigma,
			   'mean': 1.}
Params_ext 	= {'tau': dt, 'std': sigma,
			   'mean': 0.}
Params_gs 	= {'theta':0.,'alpha':1.}

# neurons 
LLBN 		= [None]*4
EBN  		= [None]*4
IBN  		= [None]*4
TN   		= [None]*4
OPN 		= nest.Create('lin_rate_ipn', 			            
			params = Params_opn)
SC  		= nest.Create('lin_rate_ipn', 
			params = Params_sc)

for i in range(4):
    LLBN[i]     = nest.Create('lin_rate_ipn', 						
				params = Params_llbn)
    EBN[i]      = nest.Create('lin_rate_ipn',						
				params = Params_ebn)
    IBN[i]      = nest.Create('lin_rate_ipn',						 
				params = Params_ibn)
    TN[i]       = nest.Create('lin_rate_ipn',						
				params = Params_tn)

# bias unit (sends constant activity)
Bias 		= nest.Create('lin_rate_ipn', 									
			params = Params_bias)	

# unit representing external stimulation to OPN
EXT 		= nest.Create('lin_rate_ipn', 									
			params = Params_ext)

# auxiliary units (apply nonlinearities)
gS  		= nest.Create('rate_transformer_threshold_lin',
			params = Params_gs)				                
gP  		= nest.Create('rate_transformer_sigmoid_gg_1998') 				     
gL  		= [None]*4
for i in range(0,4):
    gL[i] 	= nest.Create('rate_transformer_sigmoid_gg_1998') 		
    

###########################################
#### 		connections		 ##
###########################################

k 		= [1,0,3,2]
for i in range(0,4):
# to LLBNs
    nest.Connect(IBN[i], LLBN[i], 'all_to_all', {
                'model': 'rate_connection_instantaneous', 'weight': -2.0})
# to EBNs
    nest.Connect(LLBN[i], EBN[i], 'all_to_all', {
                'model': 'rate_connection_instantaneous', 'weight': 5.0})
    nest.Connect(Bias, EBN[i], 'all_to_all', {
                'model': 'rate_connection_instantaneous', 'weight': 1.0})
    nest.Connect(LLBN[k[i]], EBN[i], 'all_to_all', {
                'model': 'rate_connection_instantaneous', 'weight': -10.0})
    nest.Connect(gP, EBN[i], 'all_to_all', {
                'model': 'rate_connection_instantaneous', 'weight': -20.0}) 
# to IBNs
    nest.Connect(EBN[i], IBN[i], 'all_to_all', {
                'model': 'rate_connection_instantaneous', 'weight': 3.0})
# to TNs
    nest.Connect(EBN[i], TN[i], 'all_to_all', {
                'model': 'rate_connection_instantaneous', 'weight': .1})
    nest.Connect(EBN[k[i]], TN[i], 'all_to_all', {
                'model': 'rate_connection_instantaneous', 'weight': -.1})
# to OPN
    nest.Connect(gL[i], OPN, 'all_to_all', {
                'model': 'rate_connection_instantaneous', 'weight': -1.0})

# to auxiliary units
    nest.Connect(LLBN[i], gL[i], 'all_to_all', {
                'model': 'rate_connection_instantaneous', 'weight': 1.0})

# to OPN (cont'd)
nest.Connect(Bias, OPN, 'all_to_all', {
                'model': 'rate_connection_instantaneous', 'weight': 1.2})
nest.Connect(EXT, OPN, 'all_to_all', {
                'model': 'rate_connection_instantaneous', 'weight': 1.})

# to auxiliary units (cont'd)
nest.Connect(OPN, gP, 'all_to_all', {
                'model': 'rate_connection_instantaneous', 'weight': 1.0})
nest.Connect(SC, gS, 'all_to_all', {
                'model': 'rate_connection_instantaneous', 'weight': 1.0})
