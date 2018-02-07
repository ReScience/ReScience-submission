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
# Generates Fig. 6: Effect of stimulation frequency on saccade amplitude, duration, and velocity (figure 9 in original publication).  
# SC stimulation frequency F was varied between 1 and 2.4 at increments of 0.2. The weight W was set equal to 2, and stimulation duration was 125 ms.
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
cm2inch		= .394			# inch/cm
t_relax 	= 100  			# ms
g_pos		= 260.			# gain eye position	
threshold	= 30. 			# eye velocity marking onset of saccade		
Amplitude	= np.zeros(8)	# saccade amplitude
Duration	= np.zeros(8)	# saccade duration
Velocity	= np.zeros(8)	# saccade velocity

# figure setup
rcParams.update({'figure.autolayout': True})

fig_name	= '../article/figures/fig6.eps'
fig_size 	= np.multiply([17.6,8.5],cm2inch)
fig_rows 	= 1
fig_cols 	= 3
fig_plots	= fig_rows*fig_cols
ppi			= 1200
face	 	= 'white'
edge	 	= 'white'

ax 		 	= [None]*fig_plots
fig 		= pl.figure(facecolor = face, edgecolor = edge, figsize = fig_size)
for i in range(0,fig_plots):
	ax[i] 	= fig.add_subplot(fig_rows,fig_cols,i+1)
	ax[i].set_xlabel('frequency')
	ax[i].set_xlim([.75,2.75])
	ax[i].set_xticks([1.0,1.5,2.0,2.5])
	ax[i].locator_params(axis='y',nbins=5)
	ax[i].tick_params(right='off')
	ax[i].tick_params(top='off')
	ax[i].spines["right"].set_visible(False)
	ax[i].spines["top"].set_visible(False)

ax[0].text(-0.075, 1.1, 'A', transform=ax[0].transAxes, size=16, weight='bold')
ax[1].text(-0.075, 1.1, 'B', transform=ax[1].transAxes, size=16, weight='bold')
ax[2].text(-0.075, 1.1, 'C', transform=ax[2].transAxes, size=16, weight='bold')
###########################################
#### 		set up experiment			 ##
###########################################

# input & external electric stimulation
F			= np.linspace(1.,2.4,8)
W 			= 2.
J   		= 0.
nest.Connect(gS, LLBN[1], 'all_to_all', {'model': 'rate_connection_instantaneous', 'weight': W})

# timing protocol (in ms)
preStim  	=   0
Stim     	= 125 
postStim 	=  20

# time vector T
t_start 	= 0
t_end   	= preStim + Stim + postStim
t_steps 	= int((t_end-t_start)/dt)-1
T       	= np.linspace(t_start,t_end,t_steps)


###########################################
#### 	set up recording devices 		 ##
###########################################

MM 			= [None]*8
for s in range(0,8):
	MM[s] = nest.Create('multimeter')
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

	nest.Connect(MM[s], TN[1], syn_spec = {'delay': dt})


###########################################
#### 			simulation				 ##
###########################################

# pre-stimulus period
	nest.Simulate(preStim)

# stimulus period
	nest.SetStatus(SC,{'mean': F[s]})
	nest.SetStatus(OPN,{'mean': J})
	nest.Simulate(Stim)

# post-stimulus period
	nest.SetStatus(SC,{'mean': 0.})
	nest.Simulate(postStim)

# reset rates for next simulation
	nest.SetStatus(SC,{'rate': 0.})
	nest.SetStatus(TN[0],{'rate': .5})
	nest.SetStatus(TN[1],{'rate': .5})


###########################################
#### 			create figure			 ##
###########################################

# gather data from recording device
	data 	 = nest.GetStatus(MM[s])
	senders  = data[0]['events']['senders']
	voltages = data[0]['events']['rate']

# compute output variables
	a   	 = g_pos*(voltages[np.where(senders == TN[1])]-.5)  # amplitude (degree)	
	v   	 = np.diff(a)/(1e-3*dt)						   		# velocity (degree/sec)

	greater  = np.where(v>threshold)
	lesser   = np.where(v<threshold)
	on  	 = greater[0][0]							   		# saccade onset
	ID 		  = np.where(lesser[0][:]>greater[0][0])						 
	if (ID[0].size==0):								 	   		# if velocity did not drop below threhsold
		off 	  = argrelextrema(v, np.less)			   		# find local minima to identify offset
		off 	  = off[0][(len(off[0]))-1]	
	else:
		
		off 	  = lesser[0][ID[0][0]] 				   		# offset is determined by first moment velocity if drops below threshold	    
	Amplitude[s] = a[off]-a[on]
	Duration[s]  = T[off]-T[on]	
	Velocity[s]  = np.max(v)

ax[0].plot(F,Amplitude,'ko-',linewidth=2)
ax[0].set_ylabel('amplitude (deg)')
ax[0].set_ylim([20.,50.])

ax[1].plot(F,Duration,'ko-',linewidth=2)
ax[1].set_ylabel('duration (ms)')
ax[1].set_ylim([50.,120.])

ax[2].plot(F,Velocity,'ko-',linewidth=2)
ax[2].set_ylabel('peak velocity (deg/s)')
ax[2].set_ylim([450.,650.])

pl.savefig(fig_name, format='eps', dpi=ppi, bbox_inches='tight')
pl.show()


