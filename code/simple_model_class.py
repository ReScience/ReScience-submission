# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2017, Christoph Metzner
# Distributed under the (new) BSD License.
#
# Contributors: Christoph Metzner (c.metzner@herts.ac.uk)
# -----------------------------------------------------------------------------
# References:
#
# * Vierling-Claassen, D., Siekmeier, P., Stufflebeam, S., & Kopell, N. (2008). 
#   Modeling GABA alterations in schizophrenia: a link between impaired 
#   inhibition and altered gamma and beta range auditory entrainment. 
#   Journal of neurophysiology, 99(5), 2656-2671.
# -----------------------------------------------------------------------------
# The main class implementing the model of the replication study.
#
# -----------------------------------------------------------------------------
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab





class simpleModel(object):   
    '''The simple model from Vierling-Claassen et al. (J Neurophysiol, 2008)

     Attributes
    -----------------
    n_ex        : int
        number of excitatory cells
    n_inh        : int
        number of inhibitory cells

    eta        : float
        synaptic scaling factor
    tau_R        : float
        synaptic rise time
    tau_ex        : float
        exc. synaptic decay time
    tau_inh        : inh. 'synaptic decay time' 

    g_ee        : float
        E-E weight
    g_ei        : float
        E-I weight
    g_ie        : float
        I-E weight
    g_ii        : float
        I-I weight
    g_de        : float
         Drive-E weight
    g_di        : float
        Drive-I weight 

    dt        : float
        time step

    b_ex        : float
        applied current to excitatory cells
    b_inh        : float
        applied current to inhibitory cells
    drive_frequency : float
        drive frequency (a drive frequency of 0.0 means no drive at all and the
        network is solely driven by noise)
    
    background_rate : float
        rate of the background noise spike trains
    A        : float
        scaling factor for the background noise strength

    seed        : int
        seed for the random number generator
    '''

    def __init__(self,n_ex=20,n_inh=10,eta=5.0,tau_R=0.1,tau_ex=2.0,tau_inh=8.0,
g_ee=0.015,g_ei=0.025,g_ie=0.015,g_ii=0.02,g_de=0.3,g_di=0.08,dt=0.05,
b_ex=-0.01,b_inh=-0.01,drive_frequency=0.0,background_rate=33.3,A=0.5,
seed=12345,filename='default',directory='/'):
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
        '''Runs the model and returns (and stores) the results
            
        Parameters
        -----------------
        time    : float
            The duration of the simulation.
        saveMEG : int
            A flag that signalises whether the MEG signal should be stored
        saveEX: : int
            A flag that signalises whether the exc. population activity 
           should be stored
        saveINH : int
            A flag that signalises whether the inh.population activity 
            should be stored
        '''
        # number of time steps 
        time_points = np.linspace(0,time,int(time/self.dt)) 
        
        # Initialisations

        # the pacemaking drive cell
        drive_cell  =   np.zeros((len(time_points),))    
        
        # exc. neurons
        theta_ex = np.zeros((self.n_ex,len(time_points)))        
        # inh. neurons
        theta_inh = np.zeros((self.n_inh,len(time_points)))        
        
        # E-E snyaptic gating variables
        s_ee = np.zeros((self.n_ex,self.n_ex,len(time_points))) 
        # E-I snyaptic gating variables    
        s_ei = np.zeros((self.n_ex,self.n_inh,len(time_points)))    
        # I-E snyaptic gating variables
        s_ie = np.zeros((self.n_inh,self.n_ex,len(time_points)))
        # I-I snyaptic gating variables    
        s_ii = np.zeros((self.n_inh,self.n_inh,len(time_points)))    
        # Drive-E snyaptic gating variables
        s_de = np.zeros((self.n_ex,len(time_points)))  
        # Drive-I snyaptic gating variables          
        s_di = np.zeros((self.n_inh,len(time_points)))
        
        # Noise to exc. cells
        N_ex = np.zeros((self.n_ex,len(time_points)))  
        # Noise to inh. cells          
        N_inh = np.zeros((self.n_inh,len(time_points)))            
        
        # Synaptic inputs for exc. cells
        S_ex = np.zeros((self.n_ex,len(time_points)))
        # Synaptic inputs for inh. cells            
        S_inh = np.zeros((self.n_inh,len(time_points)))            
        
        # MEG component for each cell (only E-E EPSCs)
        meg = np.zeros((self.n_ex,len(time_points)))            
        
        # applied current for exc. cells
        B_ex    = self.b_ex * np.ones((self.n_ex,)) 
        # applied current for inh. cells               
        B_inh   = self.b_inh * np.ones((self.n_inh,))            

        # Frequency = 1000/period(in ms) and b= pi**2 / period**2 
        # (because period = pi* sqrt(1/b); see Boergers and Kopell 2003) 
        period  = 1000.0/self.drive_frequency
        # applied current for drive cell
        b_drive = np.pi**2/period**2             
        
        # Seed the random generator
        random.seed( self.seed )
        
        # Noise spike trains
        ST_ex = [None]*self.n_ex
        ST_inh = [None]*self.n_inh
        
        # adjust rate to ms time scale
        rate_parameter = self.background_rate/1000.0 
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

            # evolve gating variable
       
            # E-E connections
            # exponential decay
            ee_decay            = s_ee[:,:,t-1]/self.tau_ex 
            # synaptic input from other cells
            ee_synaptic_input   = np.exp(-1.0*self.eta*(1+np.cos(theta_ex[:,t-1])))*((1.0-s_ee[:,:,t-1])/self.tau_R) 
            s_ee[:,:,t] = s_ee[:,:,t-1]+self.dt*(-1.0*ee_decay+ee_synaptic_input)
             
            # E-I connections          
            for k in range(self.n_ex):
                a[k,0]=theta_ex[k,t-1]
            # exponential decay
            ei_decay            = s_ei[:,:,t-1]/self.tau_ex 
            # synaptic input from other cells
            ei_synaptic_input   = np.exp(-1.0*self.eta*(1+np.cos(a)))*((1.0-s_ei[:,:,t-1])/self.tau_R) 
            s_ei[:,:,t] = s_ei[:,:,t-1]+self.dt*(-1.0*ei_decay+ei_synaptic_input)
                      
            # I-E connections 
            for l in range(self.n_inh):
                b[l,0]=theta_inh[l,t-1]
            # exponential decay
            ie_decay            = s_ie[:,:,t-1]/self.tau_inh 
            # synaptic input from other cells
            ie_synaptic_input   = np.exp(-1.0*self.eta*(1+np.cos(b)))*((1.0-s_ie[:,:,t-1])/self.tau_R) 
            s_ie[:,:,t] = s_ie[:,:,t-1]+self.dt*(-1.0*ie_decay+ie_synaptic_input)

            # I-I connections
            # exponential decay
            ii_decay            = s_ii[:,:,t-1]/self.tau_inh 
            # synaptic input from other cells
            ii_synaptic_input   = np.exp(-1.0*self.eta*(1+np.cos(theta_inh[:,t-1])))*((1.0-s_ii[:,:,t-1])/self.tau_R) 
            s_ii[:,:,t] = s_ii[:,:,t-1]+self.dt*(-1.0*ii_decay+ii_synaptic_input)

            # D-E connections
            # exponential decay
            de_decay            = s_de[:,t-1]/self.tau_ex 
            # synaptic input from drive cell
            de_synaptic_input   = np.exp(-1.0*self.eta*(1+np.cos(drive_cell[t-1])))*((1.0-s_de[:,t-1])/self.tau_R) 
            s_de[:,t] = s_de[:,t-1]+self.dt*(-1.0*de_decay+de_synaptic_input)

            # D-I connections
            # exponential decay
            di_decay            = s_di[:,t-1]/self.tau_ex 
            # synaptic input from drive cell
            di_synaptic_input   = np.exp(-1.0*self.eta*(1+np.cos(drive_cell[t-1])))*((1.0-s_di[:,t-1])/self.tau_R) 
            s_di[:,t] = s_di[:,t-1]+self.dt*(-1.0*di_decay+di_synaptic_input)
             
            # calculate total synaptic input
            excitation = self.g_ee*np.sum(s_ee[:,:,t-1],axis=0)
            inhibition = self.g_ie*np.sum(s_ie[:,:,t-1],axis=0)
            drive = self.g_de*s_de[:,t-1]
            S_ex[:,t] = excitation-inhibition+drive

            excitation = self.g_ei*np.sum(s_ei[:,:,t-1],axis=0)
            inhibition = self.g_ii*np.sum(s_ii[:,:,t-1],axis=0)
            drive = self.g_di*s_di[:,t-1]
            S_inh[:,t] = excitation-inhibition+drive
             
            meg[:,t]    = self.g_ee*np.sum(s_ee[:,:,t-1],axis=0)  

 
            # evolve drive cell
            part_a = (1-np.cos(drive_cell[t-1]))
            part_b = b_drive*(1 + np.cos(drive_cell[t-1])) 
            drive_cell[t] = drive_cell[t-1]+self.dt*(part_a+part_b)
             
            # evolve theta
            part_a = (1 - np.cos(theta_ex[:,t-1]))
            part_b = (B_ex + S_ex[:,t] + N_ex[:,t])*(1 + np.cos(theta_ex[:,t-1]))
            theta_ex[:,t]= theta_ex[:,t-1]  + self.dt*(part_a+part_b)

            part_a = (1 - np.cos(theta_inh[:,t-1]))
            part_b = (B_inh + S_inh[:,t] + N_inh[:,t])*(1 + np.cos(theta_inh[:,t-1]))
            theta_inh[:,t] = theta_inh[:,t-1] + self.dt*(part_a+part_b)
        
        
        
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
        '''Plots a simulated neuron trace versus time.

         Parameters
        -----------------
        trace      : ndarray
            1D array containing the simulated neuron trace.
        sim_time : float
            The duration of the simulation.
        dt       : float
            The time step.
        save     : int
            A flag whether to save the plot or not.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        time = np.linspace(0,sim_time,int(sim_time/self.dt))
        ax.plot(time,trace,'k')
        

    
    def plotMEG(self,MEG,sim_time,save):
        '''Plots a simulated MEG signal versus time.

         Parameters
        -----------------
        meg      : ndarray
            1D array containing the simulated MEG signal.
        sim_time : float
            The duration of the simulation.
        dt       : float
            The time step.
        save     : int
            A flag whether to save the plot or not.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        time = np.linspace(0,sim_time,int(sim_time/self.dt))
        ax.plot(time,MEG,'k')
        
        if save:
            filenamepng = self.directory+self.filename+'-MEG.png'

            plt.savefig(filenamepng,dpi=600)
        
        
    
    def rasterPlot(self,data,sim_time,save,name):
        '''Plots a raster plot for an array of spike trains.

         Parameters
        -----------------
        spike_times   : list
            A list containing lists of spike times.
        sim_time      : float 
            The duration of the simulation.

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

    def calculatePSD(self,meg,sim_time):
        '''Calculates the power spectral density of a simulated MEG 
        signal (using mlab.psd(), which uses Welch's method).

        Parameters
        -----------------
        meg      : ndarray
            1D array containing the simulated MEG signal.
        sim_time : float
            The duration of the simulation.

        Returns
        -----------------
        ndarray,ndarray
            Two 1D arrays containing the power spectral density of 
            the signal and the according frequencies.
        '''
        # fourier sample rate
        fs = 1. / self.dt    
        
        tn = np.linspace(0,sim_time,int(sim_time/self.dt)+1)        
        npts = len(meg)   
        
        pxx,freqs=mlab.psd(meg,NFFT=npts,Fs=fs,noverlap=0,window=mlab.window_none)
        pxx[0] = 0.0
            
        return pxx,freqs
    
           
           
    def plotPSD(self,freqs,psd,fmax,save):
        '''Plots the power spectral density of a simulated MEG signal.

         Parameters
        -----------------
        freqs   : ndarray
            1D array containing the frequencies.
        pxx     : ndarray
            1D array containing the power spectral desnsity.
        fmax    : float 
            Maximal frequency to display.
        save    : int
            Flag whether to save the plot or not.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)   
      
        ax.plot(freqs*1000,psd) # adjust for ms time scale of data
        ax.axis(xmin=0, xmax=fmax)
            
        if save:
            filenamepng = self.directory+self.filename+'-PSD.png'

            plt.savefig(filenamepng,dpi=600)            
        
        return ax
    
        
    def _getSingleSpikeTimes(self,neuron):
        '''Calculates the spike times from the trace of a single theta neuron.

         Parameters
        -----------------
        neuron : ndarray
            1D array containing the single neuron trace.

        Returns
        -----------------
        list
            A list containing the spike times.
        '''
        spike_times = []
        old = 0.0
        for i,n in enumerate(neuron):

            # if theta passes (2l-1)*pi, l integer, with dtheta/dt>0 then 
            # the neuron spikes (see Boergers and Kopell, 2003)
            # therefore we loop through the time series of the single neuron 
            # and look for these events the spike time is then simply the product
            # of the index times the time step 
            if (n%(2*np.pi))>np.pi and (old%(2*np.pi))<np.pi:
                spike_time = i*self.dt
                spike_times.append(spike_time)
            old = n
        
        return spike_times
        
    def _getSpikeTimes(self,data):
        '''Calculates the spike times from an array of theta neuron traces.

         Parameters
        -----------------
        data   : ndarray
            nD array containing the traces.
        dt     : float 
            The time step.

        Returns
        -----------------
        list
            A list containing lists of spike times.
        '''
        nx,ny = data.shape
        spike_times_array = [None]*nx
        for i in range(nx):
            spike_times_array[i] = self._getSingleSpikeTimes(data[i,:])
        
        return spike_times_array
    
    def _noise(self,t,tn):
        '''Calculates the noise EPSP according to the formula from the model 
        description of the article.

         Parameters
        -----------------
        t   : int
            Current time.
        tn  : int 
            The time of the noise spike.

        Returns
        -----------------
        list
            A list containing lists of spike times.
        '''
        t  = t * self.dt
        if t-tn>0:
            exp1 = np.exp(-(t-tn)/self.tau_ex)
            exp2 = np.exp(-(t-tn)/self.tau_R)
            value = (self.A*(exp1-exp2))/(self.tau_ex-self.tau_R)
        else:
            value = 0
    
        return value









