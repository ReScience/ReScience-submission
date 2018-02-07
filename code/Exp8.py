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
# Generates Fig. 8: Smooth staircase eye movement (figure 11 in original publication).  
# Input I to the left side of the SG was set equal to 3 for 300 ms
# -----------------------------------------------------------------------------

import pylab as pl
from matplotlib import rcParams


###########################################
#### 		set up model		 ##
###########################################

execfile('setup_model.py')


###########################################
#### 		auxiliary		 ##
###########################################

# additional variables
cm2inch		= .394	# inch/cm
t_relax 	= 100  	# ms

# figure setup
rcParams.update({'figure.autolayout': True})

fig_name	= '../article/figures/fig8.eps'
fig_size 	= np.multiply([17.6,11.6],cm2inch)
fig_rows 	= 5
fig_cols 	= 1
fig_plots	= fig_rows*fig_cols
ppi		= 1200
face	 	= 'white'
edge	 	= 'white'

ax 		= [None]*fig_plots
fig 		= pl.figure(facecolor = face, edgecolor = edge, figsize = fig_size)
for i in range(0,fig_plots):
	ax[i] 	= fig.add_subplot(fig_rows,fig_cols,i+1)
	ax[i].get_xaxis().set_visible(False)
	ax[i].get_yaxis().set_ticks([])
	ax[i].set_ylim([-1.,1.5])
	ax[i].tick_params(right='off')
	ax[i].tick_params(top='off')


###########################################
#### 	     set up experiment		 ##
###########################################

# input & external electrict stimulation
I  		= 3.
J   		= 0.

# timing protocol (in ms)
preStim  	= 100
Stim     	= 300
postStim 	=  50

# time vector T
t_start 	= 0
t_end   	= preStim + Stim + postStim
t_steps 	= int((t_end-t_start)/dt)-1
T       	= np.linspace(t_start,t_end,t_steps)

###########################################
#### 	set up recording devices 	 ##
###########################################

multimeter	= nest.Create('multimeter')
nest.SetStatus(multimeter, {'interval': dt, 'record_from': ['rate']})


###########################################
#### 		relaxation		 ##
###########################################

# let system reach equilibrium
# in the absence of input and stimulation
nest.Simulate(t_relax)


###########################################
#### 	connect recording devices  	 ##
###########################################

nest.Connect(multimeter, OPN, syn_spec	   = {'delay': dt})
nest.Connect(multimeter, LLBN[0], syn_spec = {'delay': dt})
nest.Connect(multimeter, EBN[0], syn_spec  = {'delay': dt})
nest.Connect(multimeter, IBN[0], syn_spec  = {'delay': dt})
nest.Connect(multimeter, TN[0], syn_spec   = {'delay': dt})


###########################################
#### 		simulation		 ##
###########################################

# pre-stimulus period
nest.Simulate(preStim)

# stimulus period
nest.SetStatus(LLBN[0],{'mean': I})
nest.SetStatus(OPN,{'mean': J})
nest.Simulate(Stim)

# post-stimulus period
nest.SetStatus(LLBN[0],{'mean': 0.})
nest.Simulate(postStim)


###########################################
#### 		create figure		 ##
###########################################

# gather data from recording device
data		= nest.GetStatus(multimeter)
senders		= data[0]['events']['senders']
voltages	= data[0]['events']['rate']

Input		= np.zeros(t_end-t_start)
Input[preStim+1:preStim+Stim] = I

# plot
ax[0].plot(range(t_start,t_end),Input,'k')
ax[0].set_ylabel('Input')
ax[0].set_ylim([-.5,3.5])

ax[1].plot(T,voltages[np.where(senders == LLBN[0])],'k')
ax[1].set_ylabel('LLBN')

ax[2].plot(T,voltages[np.where(senders == EBN[0])],'k')
ax[2].set_ylabel('EBN')

ax[3].plot(T,voltages[np.where(senders == OPN)],'k')
ax[3].set_ylabel('OPN')

ax[4].plot(T,voltages[np.where(senders == TN[0])],'k')
ax[4].get_xaxis().set_visible(True)
ax[4].set_xlabel('time (ms)')
ax[4].set_ylabel('TN')
ax[i].set_ylim([.4,1.25])

pl.savefig(fig_name, format='eps', dpi=ppi, bbox_inches='tight')
pl.show()


