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
# Simulates reimplemented saccade generation (SG) model of Gancarz & Grossberg (1998).
# Generates Fig. 4: Saccadic staircase (figure 7 in original publication). 
# Inputs to the horizontal and vertical circuits were held at (0.2, 0.33) for 250 ms
# -----------------------------------------------------------------------------

import pylab as pl
from matplotlib import rcParams
from scipy.signal import argrelextrema

###########################################
#### 		set up model				 ##
###########################################

execfile('setup_model.py')


###########################################
#### 			auxiliary				 ##
###########################################

# additional variables
cm2inch		= .394 # inch/cm
t_relax 	= 100  # ms
g_pos		= 260. # gain eye position
sr 			= 40   # sampling rate

# figure setup
rcParams.update({'figure.autolayout': True})

fig_name	= '../article/figures/fig4.eps'
fig_size 	= np.multiply([8.5,11.6],cm2inch)
ppi			= 1200
face	 	= 'white'
edge	 	= 'white'
fig 		= pl.figure(facecolor = face, edgecolor = edge, figsize = fig_size)
ax 			= fig.add_subplot(1,1,1)

ax.set_xlim([0.,10.])
ax.set_ylim([0.,15.5])
ax.set_xlabel('horizontal eye position (deg)')
ax.set_ylabel('vertical eye position (deg)')


###########################################
#### 		set up experiment			 ##
###########################################

# input & external electric stimulation
I_horizontal = .20
I_vertical   = .33 
J   		 = .0

# timing protocol (in ms)
preStim  	 =   0
Stim     	 = 250
postStim 	 =  50

# time vector T
t_start 	 = 0
t_end   	 = preStim + Stim + postStim
t_steps 	 = int((t_end-t_start)/dt)-1
T       	 = np.linspace(t_start,t_end,t_steps)


###########################################
#### 	set up recording devices 		 ##
###########################################

multimeter	 = nest.Create('multimeter')
nest.SetStatus(multimeter, {'interval': dt, 'record_from': ['rate']})


###########################################
#### 			relaxation				 ##
###########################################

# let system reach equilibrium
# in the absence of input and stimulation
nest.Simulate(t_relax)

###########################################
#### 	connect recording devices  		 ##
###########################################

nest.Connect(multimeter, TN[1], syn_spec = {'delay': dt})
nest.Connect(multimeter, TN[3], syn_spec = {'delay': dt})


###########################################
#### 			simulation				 ##
###########################################

# pre-stimulus period
nest.Simulate(preStim)

# stimulus period
nest.SetStatus(LLBN[1],{'mean': I_horizontal})
nest.SetStatus(LLBN[3],{'mean': I_vertical})
nest.SetStatus(OPN,{'mean': J})
nest.Simulate(Stim)

# post-stimulus period
nest.SetStatus(LLBN[1],{'mean': 0.})
nest.SetStatus(LLBN[3],{'mean': 0.})
nest.Simulate(postStim)


###########################################
#### 			create figure			 ##
###########################################

# gather data from recording device
data 	 	= nest.GetStatus(multimeter)
senders  	= data[0]['events']['senders']
voltages 	= data[0]['events']['rate']

# compute output variables (horizontal and vertical eye position)
theta_h 	= g_pos*(voltages[np.where(senders == TN[1])]-.5)
theta_v 	= g_pos*(voltages[np.where(senders == TN[3])]-.5)

# plot
ax.plot(theta_h[::sr],theta_v[::sr],'k.')

pl.savefig(fig_name, format='eps', dpi=ppi, bbox_inches='tight')
pl.show()



