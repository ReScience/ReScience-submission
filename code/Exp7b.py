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
# As Exp7a with the difference that generation of the input is adjusted to exactly reproduce
# the curve shown in the original publication (by reducing the time constant after stimulus offset). 
# -----------------------------------------------------------------------------

import pylab as pl
from matplotlib import rcParams
from itertools import cycle


###########################################
#### 		set up model				 ##
###########################################

execfile('setup_model.py')


###########################################
#### 			auxiliary				 ##
###########################################

# additional variables
cm2inch		= .394	# inch/cm
t_relax 	= 100   # ms

# figure setup
rcParams.update({'figure.autolayout': True})

fig_name	= '../article/figures/fig7b.eps'
fig_size 	= np.multiply([11.6,17.6],cm2inch)
fig_rows 	= 5
fig_cols 	= 1
fig_plots	= fig_rows*fig_cols
ppi			= 1200
face	 	= 'white'
edge	 	= 'white'

ax 		 	= [None]*fig_plots
fig 		= pl.figure(facecolor = face, edgecolor = edge, figsize = fig_size)
for i in range(0,fig_plots):
	ax[i] 	= fig.add_subplot(fig_rows,fig_cols,i+1)
	ax[i].get_xaxis().set_visible(False)
	ax[i].get_yaxis().set_ticks([])
	ax[i].set_ylim([-.5,2.])
	ax[i].tick_params(right='off')
	ax[i].tick_params(top='off')

lines 	   	= ['k--','k-']
linecycler 	= cycle(lines)


###########################################
#### 		set up experiment			 ##
###########################################
# input & external electrict stimulation
F			= [3.,1.3]
W 			= 2.
J   		= 0.
nest.Connect(gS, LLBN[0], 'all_to_all', {'model': 'rate_connection_instantaneous', 'weight': W})

# timing protocol (in ms)
preStim  	= 50
Stim		= [68,117]
postStim 	= [132,83]

# time vector T
t_start 	= 0
t_end   	= preStim + Stim[0] + postStim[0]
t_steps 	= int((t_end-t_start)/dt)-1
T       	= np.linspace(t_start,t_end,t_steps)


###########################################
#### 	set up recording devices 		 ##
###########################################

MM 			= [None]*2
for s in range(0,2):
	MM[s] 		= nest.Create('multimeter')
	nest.SetStatus(MM[s], {'interval': dt, 'record_from': ['rate']})


###########################################
#### 			relaxation				 ##
###########################################

# let system reach equilibrium
# in the absence of input and stimulation
	nest.Simulate(t_relax)


###########################################
#### 	connect recording devices  		 ##
###########################################

	nest.Connect(MM[s], gS, syn_spec	  = {'delay': dt})
	nest.Connect(MM[s], OPN, syn_spec	  = {'delay': dt})
	nest.Connect(MM[s], LLBN[0], syn_spec = {'delay': dt})
	nest.Connect(MM[s], EBN[0], syn_spec  = {'delay': dt})
	nest.Connect(MM[s], TN[0], syn_spec   = {'delay': dt})


###########################################
#### 			simulation				 ##
###########################################

# pre-stimulus period
	nest.SetStatus(SC,{'mean': 0.})
	nest.Simulate(preStim)

# stimulus period
	nest.SetStatus(SC,{'mean': F[s]})
	nest.SetStatus(OPN,{'mean': J})
	nest.Simulate(Stim[s])

# post-stimulus period
	nest.SetStatus(SC,{'mean': 0.,'tau':25.})
	nest.Simulate(postStim[s])

# reset rates for next simulation
	nest.SetStatus(SC,{'rate': 0.,'tau':50.})
	nest.SetStatus(TN[0],{'rate': .5})
	nest.SetStatus(TN[1],{'rate': .5})
	

###########################################
#### 			create figure			 ##
###########################################

# gather data from recording device
	data 	 = nest.GetStatus(MM[s])
	senders  = data[0]['events']['senders']
	voltages = data[0]['events']['rate']

	l 		 = next(linecycler)	
	
	ax[0].plot(T,voltages[np.where(senders == gS)],l)
	ax[0].set_ylabel('Input')
	ax[1].plot(T,voltages[np.where(senders == LLBN[0])],l)
	ax[1].set_ylabel('LLBN')

	ax[2].plot(T,voltages[np.where(senders == EBN[0])],l)
	ax[2].set_ylabel('EBN')

	ax[3].plot(T,voltages[np.where(senders == OPN)],l)
	ax[3].set_ylabel('OPN')

	ax[4].plot(T,voltages[np.where(senders == TN[0])],l)
	ax[4].get_xaxis().set_visible(True)
	ax[4].set_xlabel('time (ms)')
	ax[4].set_ylabel('TN')
	ax[4].set_ylim([.3,.7])

ax[0].text(-0.075, 1.1, 'B', transform=ax[0].transAxes, size=16, weight='bold')
pl.savefig(fig_name, format='eps', dpi=ppi,bbox_inches='tight')
pl.show()



