import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


# Helper functions
def calcPowerSpectrum(meg,dt,sim_time):
	# fourier sample rate
  	fs = 1. / dt	

	tn = np.linspace(0,sim_time,int(sim_time/dt))

  	npts = len(meg)
 	startpt = int(0.2*fs)

  	if (npts - startpt)%2!=0:
      		startpt = startpt + 1

  	meg = meg[startpt:]
  	tn = tn[startpt:]
  	nfft = len(tn)

  	pxx,freqs=mlab.psd(meg,NFFT=nfft,Fs=fs,noverlap=0,window=mlab.window_none)
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

# Load data
ctrl_avg_meg_20 = np.load('../data/Data-Average/sims_ctrl_avg_20Hz.npy')
ctrl_avg_meg_30 = np.load('../data/Data-Average/sims_ctrl_avg_30Hz.npy')
ctrl_avg_meg_40 = np.load('../data/Data-Average/sims_ctrl_avg_40Hz.npy')

schiz_avg_meg_20 = np.load('../data/Data-Average/sims_schiz_avg_20Hz.npy')
schiz_avg_meg_30 = np.load('../data/Data-Average/sims_schiz_avg_30Hz.npy')
schiz_avg_meg_40 = np.load('../data/Data-Average/sims_schiz_avg_40Hz.npy')





# Parameters
sim_time = 500
s = 2**13
dt = float(sim_time)/float(s)
time = np.linspace(0,sim_time,int(sim_time/dt))



# Flags
savefig = 0	# set savefig to 1 if you want to store the figures
showfig = 1	# set plotfi to 1 if you want to show the figures


########### Figure 1 - Replication of Figure 4 from Vierling-Claassen et al.
#plt.figure(1)
f,((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,sharex=True,sharey=True,figsize=[15.0,15.0])

ax1.plot(time,ctrl_avg_meg_40,'k',linewidth=1.5)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Simulated MEG')
ax1.annotate('40 Hz \n Drive',xy=(0,0.5),xytext=(-ax1.yaxis.labelpad-5,0),xycoords=ax1.yaxis.label,textcoords='offset points',size=25,ha='right',va='center')
ax2.plot(time,schiz_avg_meg_40,'k',linewidth=1.5)
ax2.set_xlabel('Time (ms)')
ax2.set_ylabel('Simulated MEG')
ax3.plot(time,ctrl_avg_meg_30,'k',linewidth=1.5)
ax3.set_xlabel('Time (ms)')
ax3.set_ylabel('Simulated MEG')
ax3.annotate('30 Hz \n Drive',xy=(0,0.5),xytext=(-ax3.yaxis.labelpad-5,0),xycoords=ax3.yaxis.label,textcoords='offset points',size=25,ha='right',va='center')
ax4.plot(time,schiz_avg_meg_30,'k',linewidth=1.5)
ax4.set_xlabel('Time (ms)')
ax4.set_ylabel('Simulated MEG')
ax5.plot(time,ctrl_avg_meg_20,'k',linewidth=1.5)
ax5.set_xlabel('Time (ms)')
ax5.set_ylabel('Simulated MEG')
ax5.annotate('20 Hz \n Drive',xy=(0,0.5),xytext=(-ax5.yaxis.labelpad-5,0),xycoords=ax5.yaxis.label,textcoords='offset points',size=25,ha='right',va='center')
ax5.annotate('Control',xy=(0.5,0),xytext=(0,-75),xycoords='axes fraction',textcoords='offset points',size=25,ha='center',va='bottom')
ax6.plot(time,schiz_avg_meg_20,'k',linewidth=1.5)
ax6.set_xlabel('Time (ms)')
ax6.set_ylabel('Simulated MEG')
ax6.annotate('Schizophrenia',xy=(0.5,0),xytext=(0,-75),xycoords='axes fraction',textcoords='offset points',size=25,ha='center',va='bottom')

plt.setp(ax1.get_xticklabels(),visible=True)
plt.setp(ax2.get_xticklabels(),visible=True)
plt.setp(ax3.get_xticklabels(),visible=True)
plt.setp(ax4.get_xticklabels(),visible=True)

plt.setp(ax2.get_yticklabels(),visible=True)
plt.setp(ax4.get_yticklabels(),visible=True)
plt.setp(ax6.get_yticklabels(),visible=True)

if savefig:
	plt.savefig('Replication-Figure4.png',dpi=600)
	plt.savefig('Replication-Figure4.eps',dpi=600)


########### Figure 2 - Replication of Figure 5 from Vierling-Claassen et al.
# Calculate power spectra
pxx_ctrl_20,freqs = calcPowerSpectrum(ctrl_avg_meg_20,dt,sim_time)
pxx_schiz_20,_ = calcPowerSpectrum(schiz_avg_meg_20,dt,sim_time)
pxx_ctrl_30,_ = calcPowerSpectrum(ctrl_avg_meg_30,dt,sim_time)
pxx_schiz_30,_ = calcPowerSpectrum(schiz_avg_meg_30,dt,sim_time)
pxx_ctrl_40,_ = calcPowerSpectrum(ctrl_avg_meg_40,dt,sim_time)
pxx_schiz_40,_ = calcPowerSpectrum(schiz_avg_meg_40,dt,sim_time)

# Plot figure
#plt.figure(2)
f2,((ax11,ax22),(ax33,ax44),(ax55,ax66)) = plt.subplots(3,2,sharex=True,sharey=True,figsize=[15.0,15.0])

ax11.plot(freqs*1000,pxx_ctrl_40,'k',linewidth=2)
ax11.axis(xmin=0, xmax=55)
ax11.set_xlabel('Frequency (Hz)')
ax11.set_ylabel('Power')
ax11.annotate('40 Hz \n Drive',xy=(0,0.5),xytext=(-ax11.yaxis.labelpad-5,0),xycoords=ax11.yaxis.label,textcoords='offset points',size=25,ha='right',va='center')
xticks=ax11.xaxis.get_ticklabels()
xticks[2].set_weight('bold')
xticks[3].set_weight('bold')
xticks[4].set_weight('bold')

ax22.plot(freqs*1000,pxx_schiz_40,'k',linewidth=2)
ax22.set_xlabel('Frequency (Hz)')
ax22.set_ylabel('Power')
xticks=ax22.xaxis.get_ticklabels()
xticks[2].set_weight('bold')
xticks[3].set_weight('bold')
xticks[4].set_weight('bold')

ax33.plot(freqs*1000,pxx_ctrl_30,'k',linewidth=2)
ax33.set_xlabel('Frequency (Hz)')
ax33.set_ylabel('Power')
ax33.annotate('30 Hz \n Drive',xy=(0,0.5),xytext=(-ax33.yaxis.labelpad-5,0),xycoords=ax33.yaxis.label,textcoords='offset points',size=25,ha='right',va='center')
xticks=ax33.xaxis.get_ticklabels()
xticks[2].set_weight('bold')
xticks[3].set_weight('bold')
xticks[4].set_weight('bold')

ax44.plot(freqs*1000,pxx_schiz_30,'k',linewidth=2)
ax44.set_xlabel('Frequency (Hz)')
ax44.set_ylabel('Power')
xticks=ax44.xaxis.get_ticklabels()
xticks[2].set_weight('bold')
xticks[3].set_weight('bold')
xticks[4].set_weight('bold')

ax55.plot(freqs*1000,pxx_ctrl_20,'k',linewidth=2)
ax55.set_xlabel('Frequency (Hz)')
ax55.set_ylabel('Power')
ax55.annotate('20 Hz \n Drive',xy=(0,0.5),xytext=(-ax55.yaxis.labelpad-5,0),xycoords=ax55.yaxis.label,textcoords='offset points',size=25,ha='right',va='center')
ax55.annotate('Control',xy=(0.5,0),xytext=(0,-75),xycoords='axes fraction',textcoords='offset points',size=25,ha='center',va='bottom')
xticks=ax55.xaxis.get_ticklabels()
xticks[2].set_weight('bold')
xticks[3].set_weight('bold')
xticks[4].set_weight('bold')

ax66.plot(freqs*1000,pxx_schiz_20,'k',linewidth=2)
ax66.set_xlabel('Frequency (Hz)')
ax66.set_ylabel('Power')
ax66.annotate('Schizophrenia',xy=(0.5,0),xytext=(0,-75),xycoords='axes fraction',textcoords='offset points',size=25,ha='center',va='bottom')
xticks=ax66.xaxis.get_ticklabels()
xticks[2].set_weight('bold')
xticks[3].set_weight('bold')
xticks[4].set_weight('bold')


plt.setp(ax11.get_xticklabels(),visible=True)
plt.setp(ax22.get_xticklabels(),visible=True)
plt.setp(ax33.get_xticklabels(),visible=True)
plt.setp(ax44.get_xticklabels(),visible=True)

plt.setp(ax22.get_yticklabels(),visible=True)
plt.setp(ax44.get_yticklabels(),visible=True)
plt.setp(ax66.get_yticklabels(),visible=True)

if savefig:
	plt.savefig('Replication-Figure5.png',dpi=600)
	plt.savefig('Replication-Figure5.eps',dpi=600)


########### Figure 3 - Replication of Figure 6 from Vierling-Claassen et al.
# Load data
data_ex  = np.load('../data/Single-Trial-Data-40Hz/sims_ctrl_40Hz0-Ex.npy')
data_inh = np.load('../data/Single-Trial-Data-40Hz/sims_ctrl_40Hz0-Inh.npy')
data_meg = np.load('../data/Single-Trial-Data-40Hz/sims_ctrl_40Hz0-MEG.npy')

# Plot figure
plt.figure(3,figsize=[15.0,15.0])
ax1 = plt.subplot2grid((2,2),(0,0), rowspan=2)
ax2 = plt.subplot2grid((2,2),(0,1))
ax3 = plt.subplot2grid((2,2),(1,1))


# Raster plot
spike_times = getSpikeTimes(data_ex,dt)
for i,times in enumerate(spike_times):
		y = [i]*len(times)
		ax2.plot(times,y,linestyle='None',color='k',marker='|',markersize=10)

spike_times = getSpikeTimes(data_inh,dt)
for i,times in enumerate(spike_times):
		y = [i+20]*len(times)
		ax2.plot(times,y,linestyle='None',color='b',marker='|',markersize=10)
		ax2.axis([0,sim_time,-0.5,30])
ax2.set_xlabel('Time (ms)')
ax2.set_yticks([])
ax2.annotate('E cells',xy=(0,0.5),xytext=(-ax2.yaxis.labelpad,-75),xycoords=ax2.yaxis.label,textcoords='offset points',size=12,ha='right',va='top')
ax2.annotate('I cells',xy=(0,0.5),xytext=(-ax2.yaxis.labelpad,75),xycoords=ax2.yaxis.label,textcoords='offset points',size=12,ha='right',va='bottom')

# MEG plot
ax3.plot(time,data_meg,'k',linewidth=1.5)
ax3.set_xlabel('Time (ms)')
ax3.set_ylabel('Simulated MEG')


# Power spectrum plot
pxx,freqs = calcPowerSpectrum(data_meg,dt,sim_time)

ax1.plot(freqs*1000,pxx,'k',linewidth=2)
ax1.axis(xmin=0, xmax=55)
ax1.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('Power')

if savefig:
	plt.savefig('Replication-Figure6.png',dpi=600)
	plt.savefig('Replication-Figure6.eps',dpi=600)


########### Figure 4 - Replication of Figure 7 from Vierling-Claassen et al.
# Load data
data_ex  = np.load('../data/Single-Trial-Data-40Hz/sims_schiz_40Hz0-Ex.npy')
data_inh = np.load('../data/Single-Trial-Data-40Hz/sims_schiz_40Hz0-Inh.npy')
data_meg = np.load('../data/Single-Trial-Data-40Hz/sims_schiz_40Hz0-MEG.npy')

# Plot figure
plt.figure(4,figsize=[15.0,15.0])
ax1 = plt.subplot2grid((2,2),(0,0), rowspan=2)
ax2 = plt.subplot2grid((2,2),(0,1))
ax3 = plt.subplot2grid((2,2),(1,1))


# Raster plot
spike_times = getSpikeTimes(data_ex,dt)
for i,times in enumerate(spike_times):
		y = [i]*len(times)
		ax2.plot(times,y,linestyle='None',color='k',marker='|',markersize=10)

spike_times = getSpikeTimes(data_inh,dt)
for i,times in enumerate(spike_times):
		y = [i+20]*len(times)
		ax2.plot(times,y,linestyle='None',color='b',marker='|',markersize=10)
		ax2.axis([0,sim_time,-0.5,30])
ax2.set_xlabel('Time (ms)')
ax2.set_yticks([])
ax2.annotate('E cells',xy=(0,0.5),xytext=(-ax2.yaxis.labelpad,-75),xycoords=ax2.yaxis.label,textcoords='offset points',size=12,ha='right',va='top')
ax2.annotate('I cells',xy=(0,0.5),xytext=(-ax2.yaxis.labelpad,75),xycoords=ax2.yaxis.label,textcoords='offset points',size=12,ha='right',va='bottom')

# MEG plot
ax3.plot(time,data_meg,'k',linewidth=1.5)
ax3.set_xlabel('Time (ms)')
ax3.set_ylabel('Simulated MEG')


# Power spectrum plot
pxx,freqs = calcPowerSpectrum(data_meg,dt,sim_time)

ax1.plot(freqs*1000,pxx,'k',linewidth=2)
ax1.axis(xmin=0, xmax=55)
ax1.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('Power')

if savefig:
	plt.savefig('Replication-Figure7.png',dpi=600)
	plt.savefig('Replication-Figure7.eps',dpi=600)


########### Figure 5 - Replication of Figure 10 from Vierling-Claassen et al.
# Load data
data_ex  = np.load('../data/Single-Trial-Data-20Hz/sims_ctrl_20Hz13-Ex.npy')
data_inh = np.load('../data/Single-Trial-Data-20Hz/sims_ctrl_20Hz13-Inh.npy')
data_meg = np.load('../data/Single-Trial-Data-20Hz/sims_ctrl_20Hz13-MEG.npy')

# Plot figure
plt.figure(5,figsize=[15.0,15.0])
ax1 = plt.subplot2grid((2,2),(0,0), rowspan=2)
ax2 = plt.subplot2grid((2,2),(0,1))
ax3 = plt.subplot2grid((2,2),(1,1))


# Raster plot
spike_times = getSpikeTimes(data_ex,dt)
for i,times in enumerate(spike_times):
		y = [i]*len(times)
		ax2.plot(times,y,linestyle='None',color='k',marker='|',markersize=10)

spike_times = getSpikeTimes(data_inh,dt)
for i,times in enumerate(spike_times):
		y = [i+20]*len(times)
		ax2.plot(times,y,linestyle='None',color='b',marker='|',markersize=10)
		ax2.axis([0,sim_time,-0.5,30])
ax2.set_xlabel('Time (ms)')
ax2.set_yticks([])
ax2.annotate('E cells',xy=(0,0.5),xytext=(-ax2.yaxis.labelpad,-75),xycoords=ax2.yaxis.label,textcoords='offset points',size=12,ha='right',va='top')
ax2.annotate('I cells',xy=(0,0.5),xytext=(-ax2.yaxis.labelpad,75),xycoords=ax2.yaxis.label,textcoords='offset points',size=12,ha='right',va='bottom')

# MEG plot
ax3.plot(time,data_meg,'k',linewidth=1.5)
ax3.set_xlabel('Time (ms)')
ax3.set_ylabel('Simulated MEG')


# Power spectrum plot
pxx,freqs = calcPowerSpectrum(data_meg,dt,sim_time)

ax1.plot(freqs*1000,pxx,'k',linewidth=2)
ax1.axis(xmin=0, xmax=55)
ax1.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('Power')

if savefig:
	plt.savefig('Replication-Figure10.png',dpi=600)
	plt.savefig('Replication-Figure10.eps',dpi=600)


########### Figure 6 - Replication of Figure 11 from Vierling-Claassen et al.
# Load data
data_ex  = np.load('../data/Single-Trial-Data-20Hz/sims_schiz_20Hz13-Ex.npy')
data_inh = np.load('../data/Single-Trial-Data-20Hz/sims_schiz_20Hz13-Inh.npy')
data_meg = np.load('../data/Single-Trial-Data-20Hz/sims_schiz_20Hz13-MEG.npy')

# Plot figure
plt.figure(6,figsize=[15.0,15.0])
ax1 = plt.subplot2grid((2,2),(0,0), rowspan=2)
ax2 = plt.subplot2grid((2,2),(0,1))
ax3 = plt.subplot2grid((2,2),(1,1))


# Raster plot
spike_times = getSpikeTimes(data_ex,dt)
for i,times in enumerate(spike_times):
		y = [i]*len(times)
		ax2.plot(times,y,linestyle='None',color='k',marker='|',markersize=10)

spike_times = getSpikeTimes(data_inh,dt)
for i,times in enumerate(spike_times):
		y = [i+20]*len(times)
		ax2.plot(times,y,linestyle='None',color='b',marker='|',markersize=10)
		ax2.axis([0,sim_time,-0.5,30])
ax2.set_xlabel('Time (ms)')
ax2.set_yticks([])
ax2.annotate('E cells',xy=(0,0.5),xytext=(-ax2.yaxis.labelpad,-75),xycoords=ax2.yaxis.label,textcoords='offset points',size=12,ha='right',va='top')
ax2.annotate('I cells',xy=(0,0.5),xytext=(-ax2.yaxis.labelpad,75),xycoords=ax2.yaxis.label,textcoords='offset points',size=12,ha='right',va='bottom')

# MEG plot
ax3.plot(time,data_meg,'k',linewidth=1.5)
ax3.set_xlabel('Time (ms)')
ax3.set_ylabel('Simulated MEG')


# Power spectrum plot
pxx,freqs = calcPowerSpectrum(data_meg,dt,sim_time)

ax1.plot(freqs*1000,pxx,'k',linewidth=2)
ax1.axis(xmin=0, xmax=55)
ax1.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('Power')

if savefig:
	plt.savefig('Replication-Figure11.png',dpi=600)
	plt.savefig('Replication-Figure11.eps',dpi=600)


if showfig:
	plt.show()


