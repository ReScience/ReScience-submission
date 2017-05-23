import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


# Helper functions
def calcPowerSpectrum(meg,dt,sim_time):
	# fourier sample rate
  	fs = 1. / dt	

	tn = np.linspace(0,sim_time,int(sim_time/dt))
  	npts = len(meg)

  	pxx,freqs=mlab.psd(meg,NFFT=npts,Fs=fs,noverlap=0,window=mlab.window_none)
  	pxx[0] = 0.0
	
	return pxx,freqs


def getSingleSpikeTimes(neuron,dt):
	spike_times = []
	old = 0.0
	for i,n in enumerate(neuron):
		
		# if theta passes (2l-1)*pi, l integer, with dtheta/dt>0 then the neuron spikes (see Boergers and Kopell, 2003)
		if (n%(2*np.pi))>np.pi and (old%(2*np.pi))<np.pi:
			spike_time = i*dt
			spike_times.append(spike_time)
		old = n

	return spike_times

def getSpikeTimes(data,dt):
	nx,ny = data.shape
	spike_times = [None]*nx
	for i in range(nx):
		spike_times[i] = getSingleSpikeTimes(data[i,:],dt)

	return spike_times





# Parameters
sim_time = 500
s = 2**13
dt = float(sim_time)/float(s)
time = np.linspace(0,sim_time,int(sim_time/dt))



# Flags
savefig = 0	# set savefig to 1 if you want to store the figures
showfig = 1	# set plotfi to 1 if you want to show the figures

########### Figure 7 - Exploration of the influence of noise strength
# load data
ctrl_avg_meg_20_A_0 = np.load('../data/Data-Average/sims_ctrl_avg_20Hz_A0.0.npy')
schiz_avg_meg_40_A_0 = np.load('../data/Data-Average/sims_schiz_avg_40Hz_A0.0.npy')

ctrl_avg_meg_20_A_0_25 = np.load('../data/Data-Average/sims_ctrl_avg_20Hz_A0.25.npy')
schiz_avg_meg_40_A_0_25 = np.load('../data/Data-Average/sims_schiz_avg_40Hz_A0.25.npy')

ctrl_avg_meg_20_A_0_375 = np.load('../data/Data-Average/sims_ctrl_avg_20Hz_A0.375.npy')
schiz_avg_meg_40_A_0_375 = np.load('../data/Data-Average/sims_schiz_avg_40Hz_A0.375.npy')
# plot
f,((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,sharex=True,sharey=True,figsize=[19.0,15.0])

ax1.plot(time,schiz_avg_meg_40_A_0,'k',linewidth=1.5)
ax1.set_xlabel('Time (ms)',fontsize=18)
ax1.set_ylabel('Simulated MEG',fontsize=18)
ax1.annotate('0% \n Noise ',xy=(0,0.5),xytext=(-ax1.yaxis.labelpad-5,0),xycoords=ax1.yaxis.label,textcoords='offset points',size=35,ha='right',va='center')
ax2.plot(time,ctrl_avg_meg_20_A_0,'k',linewidth=1.5)
ax2.set_xlabel('Time (ms)',fontsize=18)
ax2.set_ylabel('Simulated MEG',fontsize=18)
ax3.plot(time,schiz_avg_meg_40_A_0_25,'k',linewidth=1.5)
ax3.set_xlabel('Time (ms)',fontsize=18)
ax3.set_ylabel('Simulated MEG',fontsize=18)
ax3.annotate('50% \n Noise',xy=(0,0.5),xytext=(-ax3.yaxis.labelpad-5,0),xycoords=ax3.yaxis.label,textcoords='offset points',size=35,ha='right',va='center')
ax4.plot(time,ctrl_avg_meg_20_A_0_25,'k',linewidth=1.5)
ax4.set_xlabel('Time (ms)',fontsize=18)
ax4.set_ylabel('Simulated MEG',fontsize=18)
ax5.plot(time,schiz_avg_meg_40_A_0_375,'k',linewidth=1.5)
ax5.set_xlabel('Time (ms)',fontsize=18)
ax5.set_ylabel('Simulated MEG',fontsize=18)
ax5.annotate('75% \n Noise',xy=(0,0.5),xytext=(-ax5.yaxis.labelpad-5,0),xycoords=ax5.yaxis.label,textcoords='offset points',size=35,ha='right',va='center')
ax5.annotate('Schizophrenia 40 Hz Drive',xy=(0.5,0),xytext=(0,-75),xycoords='axes fraction',textcoords='offset points',size=35,ha='center',va='bottom')
ax6.plot(time,ctrl_avg_meg_20_A_0_375,'k',linewidth=1.5)
ax6.set_xlabel('Time (ms)',fontsize=18)
ax6.set_ylabel('Simulated MEG',fontsize=18)
ax6.annotate('Control 20 Hz Drive',xy=(0.5,0),xytext=(0,-75),xycoords='axes fraction',textcoords='offset points',size=35,ha='center',va='bottom')

plt.setp(ax1.get_xticklabels(),visible=True,fontsize=18)
plt.setp(ax2.get_xticklabels(),visible=True,fontsize=18)
plt.setp(ax3.get_xticklabels(),visible=True,fontsize=18)
plt.setp(ax4.get_xticklabels(),visible=True,fontsize=18)
plt.setp(ax5.get_xticklabels(),visible=True,fontsize=18)
plt.setp(ax6.get_xticklabels(),visible=True,fontsize=18)

plt.setp(ax2.get_yticklabels(),visible=True,fontsize=18)
plt.setp(ax4.get_yticklabels(),visible=True,fontsize=18)
plt.setp(ax6.get_yticklabels(),visible=True,fontsize=18)
plt.setp(ax1.get_yticklabels(),visible=True,fontsize=18)
plt.setp(ax3.get_yticklabels(),visible=True,fontsize=18)
plt.setp(ax5.get_yticklabels(),visible=True,fontsize=18)

if savefig:
	plt.savefig('../data/Figures/Noise-Exploration.png',dpi=600)
	plt.savefig('../data/Figures/Noise-Exploration.eps',dpi=600)

if showfig:
	plt.show()
