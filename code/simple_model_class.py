####################################################################
# Implements the simple model from Vierling-Claassen et al., 
# J Neurophysiol, 2008
#
# @author: Christoph Metzner, 03/02/2017
####################################################################


import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab





class simpleModel(object):
    
    '''
	The simple model from Vierling-Claassen et al. (J Neurophysiol, 2008)

	Attributes:
		n_ex		: number of excitatory cells
		n_inh		: number of inhibitory cells

		eta		: synaptic scaling factor
		tau_R		: 'synaptic rise time'
		tau_ex		: exc. 'synaptic decay time'
		tau_inh		: inh. 'synaptic decay time' 

		g_ee		: E-E weight
		g_ei		: E-I weight
		g_ie		: I-E weight
		g_ii		: I-I weight
		g_de		: Drive-E weight
		g_di		: Drive-I weight 

		dt		: time step

		b_ex		: applied current to excitatory cells
		b_inh		: applied current to inhibitory cells
		drive_frequency : drive frequency (a drive frequency of 0.0 means no drive at all and the network is solely driven by noise)
		

		background_rate : rate of the background noise spike trains
		A		: scaling factor for the background noise strength

		seed		: seed for the random generator
	'''

    def __init__(self,n_ex=20,n_inh=10,eta=5.0,tau_R=0.1,tau_ex=2.0,tau_inh=8.0,g_ee=0.015,g_ei=0.025,g_ie=0.015,g_ii=0.02,g_de=0.3,g_di=0.08,dt=0.05,b_ex=-0.01,b_inh=-0.01,drive_frequency=0.0,background_rate=33.3,A=0.5,seed=12345,filename='default',directory='/'):
        self.n_ex = n_ex
        self.n_inh = n_inh
        self.eta = eta
        self.tau_R = tau_R
        self.tau_ex = tau_ex
        self.tau_inh = tau_inh
        self.g_ee = g_ee
        self.g_ei = g_ei
        self.g_ie = g_ie
        self.g_ii = g_ii
        self.g_de = g_de
        self.g_di = g_di
        self.dt = dt
        self.b_ex = b_ex
        self.b_inh = b_inh
        self.drive_frequency = drive_frequency
        self.background_rate = background_rate
        self.A = A
        self.seed = seed
        self.filename = filename
        self.directory = directory
    
    def run(self,time=100.0,saveMEG=0,saveEX=0,saveINH=0):
        '''
        Runs the model and returns (and stores) the results
               
        Parameters:
        time : the length of the simulation (in ms)
        saveMEG: flag that signalises whether the MEG signal should be stored
        saveEX: flag that signalises whether the exc. population activity should be stored
        saveINH: flag that signalises whether the inh.population activity should be stored
        '''
            
        time_points = np.linspace(0,time,int(time/self.dt)) # number of time steps 
    
    	# Initialisations
        drive_cell  =   np.zeros((len(time_points),))	# the pacemaking drive cell
    
    
        theta_ex = np.zeros((self.n_ex,len(time_points)))		# exc. neurons
        theta_inh = np.zeros((self.n_inh,len(time_points)))		# inh. neurons
    
        s_ee = np.zeros((self.n_ex,self.n_ex,len(time_points))) 	# E-E snyaptic gating variables
        s_ei = np.zeros((self.n_ex,self.n_inh,len(time_points)))	# E-I snyaptic gating variables
        s_ie = np.zeros((self.n_inh,self.n_ex,len(time_points)))	# I-E snyaptic gating variables
        s_ii = np.zeros((self.n_inh,self.n_inh,len(time_points)))	# I-I snyaptic gating variables
        s_de = np.zeros((self.n_ex,len(time_points)))			# Drive-E snyaptic gating variables
        s_di = np.zeros((self.n_inh,len(time_points)))			# Drive-I snyaptic gating variables
        
        N_ex = np.zeros((self.n_ex,len(time_points)))			# Noise to exc. cells
        N_inh = np.zeros((self.n_inh,len(time_points)))			# Noise to inh. cells
        
        S_ex = np.zeros((self.n_ex,len(time_points)))			# Synaptic inputs for exc. cells
        S_inh = np.zeros((self.n_inh,len(time_points)))			# Synaptic inputs for inh. cells
        
        meg = np.zeros((self.n_ex,len(time_points)))			# MEG component for each cell (only E-E EPSCs)
        
        # applied currents
        B_ex    = self.b_ex * np.ones((self.n_ex,))				# applied current for exc. cells
        B_inh   = self.b_inh * np.ones((self.n_inh,))			# applied current for inh. cells
        
        # Frequency = 1000/period(in ms) and b= pi**2 / period**2 (because period = pi* sqrt(1/b); see Boergers and Kopell 2003) 
        period  = 1000.0/self.drive_frequency
        b_drive = np.pi**2/period**2			# applied current for drive cell 
        
        # Seed the random generator
        random.seed( self.seed )
        
        # Noise spike trains
        ST_ex = [None]*self.n_ex
        ST_inh = [None]*self.n_inh
        
        rate_parameter = 1000*(1.0/self.background_rate)
        rate_parameter = 1.0/rate_parameter
        for i in range(self.n_ex):
            template_spike_array = []
            # Produce Poissonian spike train
            total_time = 0.0
            while total_time < time:
                next_time = random.expovariate(rate_parameter)
                total_time = total_time + next_time 
                if total_time < time:
                    template_spike_array.append(total_time)
                    
            ST_ex[i] = template_spike_array
                
        
        for i in range(self.n_inh):
            template_spike_array = []
            # Produce Poissonian spike train
            total_time = 0.0
            while total_time < time:
                next_time = random.expovariate(rate_parameter)
                total_time = total_time + next_time 
                if total_time < time:
                    template_spike_array.append(total_time)
                    
            ST_inh[i] = template_spike_array
                
        a = np.zeros((self.n_ex,1))    
        b = np.zeros((self.n_inh,1)) 
        # Simulation
        for t in range(1,len(time_points)):
            # calculate noise (not done in a very efficient way!)
            for i in range(self.n_ex):
                for tn in ST_ex[i]:
                    N_ex[i,t] = N_ex[i,t] + self._noise(t,tn) 
        
            # calculate noise (not done in a very efficient way!)
            for i in range(self.n_inh):
                for tn in ST_inh[i]:
                    N_inh[i,t] = N_inh[i,t] + self._noise(t,tn) 
            # evolve gating variables
            s_ee[:,:,t] 	= s_ee[:,:,t-1] + self.dt*(-1.0*(s_ee[:,:,t-1]/self.tau_ex) + np.exp(-1.0*self.eta*(1+np.cos(theta_ex[:,t-1])))*((1.0-s_ee[:,:,t-1])/self.tau_R))
            # this seems awfully complicated            
            for k in range(self.n_ex):
                a[k,0]=theta_ex[k,t-1]
            s_ei[:,:,t] 	= s_ei[:,:,t-1] + self.dt*(-1.0*(s_ei[:,:,t-1]/self.tau_ex) + np.exp(-1.0*self.eta*(1+np.cos(a)))*((1.0-s_ei[:,:,t-1])/self.tau_R))
            # this seems awfully complicated            
            for l in range(self.n_inh):
                b[l,0]=theta_inh[l,t-1]
            s_ie[:,:,t] 	= s_ie[:,:,t-1] + self.dt*(-1.0*(s_ie[:,:,t-1]/self.tau_inh) + np.exp(-1.0*self.eta*(1+np.cos(b)))*((1.0-s_ie[:,:,t-1])/self.tau_R))
            s_ii[:,:,t] 	= s_ii[:,:,t-1] + self.dt*(-1.0*(s_ii[:,:,t-1]/self.tau_inh) + np.exp(-1.0*self.eta*(1+np.cos(theta_inh[:,t-1])))*((1.0-s_ii[:,:,t-1])/self.tau_R))
            s_de[:,t]  	= s_de[:,t-1]   + self.dt*(-1.0*(s_de[:,t-1]/self.tau_ex) + np.exp(-1.0*self.eta*(1+np.cos(drive_cell[t-1])))*((1.0-s_de[:,t-1])/self.tau_R))
            s_di[:,t]   	= s_di[:,t-1]   + self.dt*(-1.0*(s_di[:,t-1]/self.tau_ex) + np.exp(-1.0*self.eta*(1+np.cos(drive_cell[t-1])))*((1.0-s_di[:,t-1])/self.tau_R))
             
            # calculate total synaptic input
            S_ex[:,t]	= self.g_ee*np.sum(s_ee[:,:,t-1],axis=0) - self.g_ie*np.sum(s_ie[:,:,t-1],axis=0) + self.g_de*s_de[:,t-1]
            S_inh[:,t] = self.g_ei*np.sum(s_ei[:,:,t-1],axis=0) - self.g_ii*np.sum(s_ii[:,:,t-1],axis=0) + self.g_di*s_di[:,t-1]
             
            meg[:,t]	= self.g_ee*np.sum(s_ee[:,:,t-1],axis=0)  #+ self.g_de*s_de[:,t-1]

            
            # evolve drive cell
            drive_cell[t]  	= drive_cell[t-1]  + self.dt*((1 - np.cos(drive_cell[t-1])) + b_drive*(1 + np.cos(drive_cell[t-1])))
             
            # evolve theta
            theta_ex[:,t]  	= theta_ex[:,t-1]  + self.dt*( (1 - np.cos(theta_ex[:,t-1])) + (B_ex + S_ex[:,t] + N_ex[:,t])*(1 + np.cos(theta_ex[:,t-1])))
            theta_inh[:,t] 	= theta_inh[:,t-1] + self.dt*( (1 - np.cos(theta_inh[:,t-1])) + (B_inh + S_inh[:,t] + N_inh[:,t])*(1 + np.cos(theta_inh[:,t-1])))
    
    
    
        # Sum EPSCs of excitatory cells
    	MEG = np.sum(meg,axis=0)

     
           
        if saveMEG:
            filenameMEG = self.directory  + self.filename + '-MEG.npy'
            np.save(filenameMEG,MEG)

          
          
        if saveEX:
            filenameEX = self.directory  + self.filename + '-Ex.npy'
            np.save(filenameEX,theta_ex)
          
          
        if saveINH:
           filenameINH = self.directory  + self.filename + '-Inh.npy'
           np.save(filenameINH,theta_inh)
          
        return MEG,theta_ex,theta_inh
    
    
    def plotTrace(self,trace,sim_time,save):
        '''
           Plots a trace signal versus time
           Parameters:
           trace: the trace signal to plot
           sim_time: the duration of the simulation
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        time = np.linspace(0,sim_time,int(sim_time/self.dt))
        ax.plot(time,trace,'k')
        
        plt.show()
    
    def plotMEG(self,MEG,sim_time,save):
        '''
           Plots a simulated MEG signal versus time
           Parameters:
           MEG: the simulated MEG signal to plot
           sim_time: the duration of the simulation
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        time = np.linspace(0,sim_time,int(sim_time/self.dt))
        ax.plot(time,MEG,'k')
        
        if save:
            filenamepng = self.directory+self.filename+'-MEG.png'

            plt.savefig(filenamepng,dpi=600)
        
        
        #plt.show()
    
    def rasterPlot(self,data,sim_time,save,name):
        '''
           Plots a raster plot for an array of spike trains
           Parameters:
           data: array of spike trains
           sim_time: duration of the simulation
        '''
        spiketrains = self._getSpikeTimes(data)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i,times in enumerate(spiketrains):
            y = [i]*len(times)
            ax.plot(times,y,linestyle='None',color='k',marker='|',markersize=10)
            ax.axis([0,sim_time,-0.5,len(spiketrains)])
    
        if save:
            filenamepng = self.directory+self.filename+'-'+name+'-raster.png'

            plt.savefig(filenamepng,dpi=600)
        #plt.show()
    
    def calculatePSD(self,meg,sim_time):
        '''
           Calculates the power spectral density of a simulated MEG signal
           Parameters:
           meg: the simulated MEG signal
           sim_time: the duration of the simulation
        '''
        # fourier sample rate
        fs = 1. / self.dt	
        
        tn = np.linspace(0,sim_time,int(sim_time/self.dt)+1)
        
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
    
           
           
    def plotPSD(self,freqs,psd,fmax,save):
        '''
            Plots the power spectral density of a simulated MEG signal
            Parameters:
            freqs: frequency vector
            psd: power spectral density vector
            fmax: maximum frequency to display
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    
        ax.plot(freqs,psd)
        ax.axis(xmin=0, xmax=fmax)
            
        if save:
            filenamepng = self.directory+self.filename+'-PSD.png'

            plt.savefig(filenamepng,dpi=600)
            
        
        return ax
    
    	
    def _getSingleSpikeTimes(self,neuron):
        '''
           Calculates the spike times from the trace of a single theta neuron
           Parameters:
           neuron: the single neuron trace
        '''
        spike_times = []
        old = 0.0
        for i,n in enumerate(neuron):
        		
          # if theta passes (2l-1)*pi, l integer, with dtheta/dt>0 then the neuron spikes (see Boergers and Kopell, 2003)
          if (n%(2*np.pi))>np.pi and (old%(2*np.pi))<np.pi:
              spike_time = i*self.dt
              spike_times.append(spike_time)
          old = n
        
        return spike_times
        
    def _getSpikeTimes(self,data):
        '''
           Calculates the spike times from an array of theta neuron traces
           Parameters:
           data: the traces array
        '''
        nx,ny = data.shape
        spike_times_array = [None]*nx
        for i in range(nx):
            spike_times_array[i] = self._getSingleSpikeTimes(data[i,:])
        
        return spike_times_array
    
    def _noise(self,t,tn):
    	t  = t * self.dt
    	if t-tn>0:
    		value = (self.A*(np.exp(-(t-tn)/self.tau_ex)-np.exp(-(t-tn)/self.tau_R)))/(self.tau_ex-self.tau_R)
    	else:
    		value = 0
    
    	return value









