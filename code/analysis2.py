import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


def getSingleSpikeTimes(neuron,dt):
	spike_times = []
	#spike_flag = 0
	old = 0.0
	for i,n in enumerate(neuron):
		
		# if theta passes (2l-1)*pi, l integer, with dtheta/dt>0 then the neuron spikes (see Boergers and Kopell, 2003)
		if (n%(2*np.pi))>np.pi and (old%(2*np.pi))<np.pi:
			#if spike_flag == 0:
			spike_time = i*dt
			spike_times.append(spike_time)
			#	spike_flag = 1
		#else:
		#	spike_flag = 0
		old = n

	return spike_times

def getSpikeTimes(data,dt):
	nx,ny = data.shape
	spike_times = [None]*nx
	for i in range(nx):
		spike_times[i] = getSingleSpikeTimes(data[i,:],dt)

	return spike_times


def rasterPlot(spike_times,sim_time):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	for i,times in enumerate(spike_times):
		y = [i]*len(times)
		ax.plot(times,y,'ko')
		ax.axis([0,sim_time,-0.5,len(spike_times)])

	#plt.show()

def plotNeuron(data,id,sim_time,dt):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	time = np.linspace(0,sim_time,int(sim_time/dt)+1)	
	ax.plot(time,np.sin(data[id,:]),'r')
    

	#plt.show()

def plotMEGTrace(meg,sim_time,dt,save,filename):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	time = np.linspace(0,sim_time,int(sim_time/dt)+1)
	ax.plot(time,meg,'k')
 
 
	if save:
            filenamepng = filename+'-MEG.png'
            print filenamepng
            plt.savefig(filenamepng,dpi=600)
        

	#plt.show()


def calcAverageFiringRate(spike_times,sim_time):
	avg_firing_rates = np.zeros((len(spike_times),))
	for i,times in enumerate(spike_times):
		avg_firing_rates[i] = len(times)*(sim_time/1000.0)

	return avg_firing_rates

def calcPowerSpectrum(meg,dt,sim_time):
	# fourier sample rate
  	fs = 1. / dt	

	tn = np.linspace(0,sim_time,int(sim_time/dt)+1)

  	npts = len(meg)
 	startpt = int(0.2*fs)

  	if (npts - startpt)%2!=0:
      		startpt = startpt + 1

  	meg = meg[startpt:]
  	tn = tn[startpt:]
  	nfft = len(tn)
  	#overlap = nfft//2


  	pxx,freqs=mlab.psd(meg,NFFT=nfft,Fs=fs,noverlap=0,window=mlab.window_none)
  	pxx[0] = 0.0
	
	return pxx,freqs


def plotPowerSpectrum(pxx,freqs,save,filename):
    fig = plt.figure()
    ax = fig.add_subplot(111)


    ax.plot(freqs*1000,pxx)
    ax.axis(xmin=0, xmax=50)
    if save:
            filenamepng = filename+'-PSD.png'
            print filenamepng
            plt.savefig(filenamepng,dpi=600)
        
	#plt.show()


if __name__ == "__main__":

	data = np.load('test.npy')
	meg = np.load('test_40_schiz_1.npy')
	sim_time = 1000
	dt = 0.01
	
	plotNeuron(data,5,sim_time,dt)

	spike_times = getSpikeTimes(data,dt)
	print spike_times[5]
	print len(spike_times[5])
	#rasterPlot(spike_times,sim_time)
	#avg = calcAverageFiringRate(spike_times,sim_time)
	#print np.mean(avg)
	#pxx,freqs = calcPowerSpectrum(meg,dt)
	#plotPowerSpectrum(pxx,freqs)
	
