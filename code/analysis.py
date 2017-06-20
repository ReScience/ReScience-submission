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
# A collection of anaylsis methods used in the replication study.
#
# -----------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


def getSingleSpikeTimes(neuron,dt):
    '''Calculates the spike times from the trace of a single theta neuron.

     Parameters
    -----------------
    neuron : ndarray
        1D array containing the single neuron trace.
    dt     : float 
        The time step.

    Returns
    -----------------
    list
        A list containing the spike times.
    '''
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
    spike_times = [None]*nx
    for i in range(nx):
        spike_times[i] = getSingleSpikeTimes(data[i,:],dt)

    return spike_times


def rasterPlot(spike_times,sim_time):
    '''Plots a raster plot for an array of spike trains.

     Parameters
    -----------------
    spike_times   : list
        A list containing lists of spike times.
    sim_time      : float 
        The duration of the simulation.

    '''
           fig = plt.figure()
    ax = fig.add_subplot(111)
    for i,times in enumerate(spike_times):
        y = [i]*len(times)
        ax.plot(times,y,'ko')
        ax.axis([0,sim_time,-0.5,len(spike_times)])

    plt.show()


def plotMEGTrace(meg,sim_time,dt,save,filename):
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
    filename : str
        The filename to save the plot.
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    time = np.linspace(0,sim_time,int(sim_time/dt))
    ax.plot(time,meg,'k')
 
    if save:
        filenamepng = filename+'-MEG.png'
        plt.savefig(filenamepng,dpi=600)
        
    plt.show()



def calcPowerSpectrum(meg,dt,sim_time):
    '''Calculates the power spectral density of a simulated MEG signal (using mlab.psd(), which uses Welch's method).

     Parameters
    -----------------
    meg      : ndarray
        1D array containing the simulated MEG signal.
    dt       : float 
        The time step.
    sim_time : float
        The duration of the simulation.

    Returns
    -----------------
    ndarray,ndarray
        Two 1D arrays containing the power spectral density of the signal and the according frequencies.
    '''
    # fourier sample rate
      fs = 1. / dt    

    tn = np.linspace(0,sim_time,int(sim_time/dt))

      npts = len(meg)


      pxx,freqs=mlab.psd(meg,NFFT=npts,Fs=fs,noverlap=0,window=mlab.window_none)
      pxx[0] = 0.0
    
    return pxx,freqs


def plotPowerSpectrum(pxx,freqs,fmax,save,filename):
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
    filename : str
        The filename to save the plot.
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)


    ax.plot(freqs*1000,pxx)
    ax.axis(xmin=0, xmax=fmax)
    if save:
        filenamepng = filename+'-PSD.png'
        plt.savefig(filenamepng,dpi=600)
        
    plt.show()

