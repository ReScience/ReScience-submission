import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


def getSingleSpikeTimes(neuron,dt):
        '''
           Calculates the spike times from the trace of a single theta neuron
           Parameters:
           neuron: the single neuron trace
	   dt: time step
        '''
	spike_times = []
	#spike_flag = 0
	old = 0.0
	for i,n in enumerate(neuron):
		
		# if theta passes (2l-1)*pi, l integer, with dtheta/dt>0 then the neuron spikes (see Boergers and Kopell, 2003)
		if (n%(2*np.pi))>np.pi and (old%(2*np.pi))<np.pi:
			spike_time = i*dt
			spike_times.append(spike_time)
		old = n

	return spike_times

def getSpikeTimes(data,dt):
        '''
           Calculates the spike times from an array of theta neuron traces
           Parameters:
           data: the traces array
	   dt: time step
        '''
	nx,ny = data.shape
	spike_times = [None]*nx
	for i in range(nx):
		spike_times[i] = getSingleSpikeTimes(data[i,:],dt)

	return spike_times


def rasterPlot(spike_times,sim_time):
        '''
           Plots a raster plot for an array of spike trains
           Parameters:
           spike_times: array of spike trains
           sim_time: duration of the simulation
        '''
	fig = plt.figure()
	ax = fig.add_subplot(111)
	for i,times in enumerate(spike_times):
		y = [i]*len(times)
		ax.plot(times,y,'ko')
		ax.axis([0,sim_time,-0.5,len(spike_times)])

	plt.show()


def plotMEGTrace(meg,sim_time,dt,save,filename):
        '''
           Plots a simulated MEG signal versus time
           Parameters:
           meg: the simulated MEG signal to plot
           sim_time: the duration of the simulation
	   dt: time step
	   save: flag whether to save the plot or not
	   filename: filename to save the plot
        '''
	fig = plt.figure()
	ax = fig.add_subplot(111)
	time = np.linspace(0,sim_time,int(sim_time/dt))
	ax.plot(time,meg,'k')
 
 
	if save:
            filenamepng = filename+'-MEG.png'
            print filenamepng
            plt.savefig(filenamepng,dpi=600)
        

	plt.show()



def calcPowerSpectrum(meg,dt,sim_time):
        '''
           Calculates the power spectral density of a simulated MEG signal
           Parameters:
           meg: the simulated MEG signal
	   dt: time step
           sim_time: the duration of the simulation
        '''
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


def plotPowerSpectrum(pxx,freqs,fmax,save,filename):
 	'''
           Plots the power spectral density of a simulated MEG signal
           Parameters:
           freqs: frequency vector
           pxx: power spectral density vector
           fmax: maximum frequency to display
	   save: flag whether to save the plot or not
	   filename: filename to save the plot
        '''
	fig = plt.figure()
	ax = fig.add_subplot(111)


	ax.plot(freqs*1000,pxx)
	ax.axis(xmin=0, xmax=fmax)
	if save:
		filenamepng = filename+'-PSD.png'
		plt.savefig(filenamepng,dpi=600)
		
	plt.show()

