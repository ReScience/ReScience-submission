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
# Generates Fig. 1: Activity profiles in left and right SG (figure 3 in original publication). 
# Input I to the left side of the SG was set equal to 1 for 265 ms
# -----------------------------------------------------------------------------

import pylab as pl
from matplotlib import rcParams


###########################################
#### 		set up model		 ##
###########################################

execfile('setup_model.py')


###########################################
#### 			auxiliary	 ##
###########################################

# additional variables
cm2inch		= .394 # inch/cm
t_relax 	= 100  # ms

# figure setup
rcParams.update({'figure.autolayout': True})

fig_name	= '../article/figures/fig1.eps'
fig_size 	= np.multiply([17.6,11.6],cm2inch)
fig_rows 	= 6
fig_cols 	= 2
fig_plots	= fig_rows*fig_cols
ppi			= 1200
face	 	= 'white'
edge	 	= 'white'

ax 			= [None]*fig_plots
fig 		= pl.figure(facecolor = face, edgecolor = edge, figsize = fig_size)
for i in range(0,fig_plots):
	col 		= np.mod(i,fig_cols)
	ax[i] 		= fig.add_subplot(fig_rows,fig_cols,i+1)
	ax[i].get_xaxis().set_visible(False)
	ax[i].get_yaxis().set_ticks([])
	ax[i].set_ylim([-1.,1.5])
	ax[i].tick_params(right='off')
	ax[i].tick_params(top='off')
	if (col):
		ax[i].get_yaxis().set_visible(False)

ax[0].text(-0.075, 1.25, 'A', transform=ax[0].transAxes, size=16, weight='bold')
ax[1].text(-0.075, 1.25, 'B', transform=ax[1].transAxes, size=16, weight='bold')

###########################################
#### 		set up experiment	 		 ##
###########################################

# input & external electric stimulation
I  			= 1.
J   		= 0.

# timing protocol (in ms)
preStim  	=  50
Stim     	= 265
postStim 	=  85

# time vector T
t_start 	= 0
t_end   	= preStim + Stim + postStim
t_steps 	= int((t_end-t_start)/dt)-1
T       	= np.linspace(t_start,t_end,t_steps)


###########################################
#### 	set up recording devices 		 ##
###########################################

multimeter	= nest.Create('multimeter')
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

nest.Connect(multimeter, OPN, syn_spec={'delay': dt})
for i in range(0,2):
    nest.Connect(multimeter, LLBN[i], syn_spec = {'delay': dt})
    nest.Connect(multimeter, EBN[i],  syn_spec = {'delay': dt})
    nest.Connect(multimeter, IBN[i],  syn_spec = {'delay': dt})
    nest.Connect(multimeter, TN[i],   syn_spec = {'delay': dt})


###########################################
#### 			simulation				 ##
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
#### 			create figure			 ##
###########################################

# gather data from recording device
data		= nest.GetStatus(multimeter)
senders 	= data[0]['events']['senders']
voltages	= data[0]['events']['rate']

Input 		= np.zeros((2,t_end-t_start))
Input[0][preStim+1:preStim+Stim] = I

# plot
ax[0].plot(range(t_start,t_end),Input[0],'k',linewidth=2)
ax[0].set_ylabel('Input')

ax[1].plot(range(t_start,t_end),Input[1],'k',linewidth=2)

ax[2].plot(T,voltages[np.where(senders == LLBN[0])],'k',linewidth=2)
ax[2].set_ylabel('LLBN')

ax[3].plot(T,voltages[np.where(senders == LLBN[1])],'k',linewidth=2)

ax[4].plot(T,voltages[np.where(senders == EBN[0])],'k',linewidth=2)
ax[4].set_ylabel('EBN')

ax[5].plot(T,voltages[np.where(senders == EBN[1])],'k',linewidth=2)

ax[6].plot(T,voltages[np.where(senders == IBN[0])],'k',linewidth=2)
ax[6].set_ylabel('IBN')

ax[7].plot(T,voltages[np.where(senders == IBN[1])],'k',linewidth=2)

ax[8].plot(T,voltages[np.where(senders == OPN)],'k',linewidth=2)
ax[8].set_ylabel('OPN')

ax[9].plot(T,voltages[np.where(senders == OPN)],'k',linewidth=2)

ax[10].plot(T,voltages[np.where(senders == TN[0])],'k',linewidth=2)
ax[10].get_xaxis().set_visible(True)
ax[10].set_xticks(range(t_start,t_end+1,100))
ax[10].set_xlabel('time (ms)')
ax[10].set_ylabel('TN')
ax[10].set_ylim([.25,.75])

ax[11].plot(T,voltages[np.where(senders == TN[1])],'k',linewidth=2)
ax[11].get_xaxis().set_visible(True)
ax[11].set_xticks(range(t_start,t_end+1,100))
ax[11].set_xlabel('time (ms)')
ax[11].set_ylim([.25,.75])

pl.savefig(fig_name, format='eps', dpi=ppi, bbox_inches='tight')
pl.show()


