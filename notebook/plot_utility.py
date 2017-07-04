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
# Utility functions for plotting
#
# -----------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from analysis import calcPowerSpectrum,getSingleSpikeTimes,getSpikeTimes

def plot_MEG_summary(time,meg_data,xlabel,ylabel,annotations,savefig,filename,fontsizes):
    f,((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,sharex=True,sharey=True,figsize=[20.0,25.0])

    ax1.plot(time,meg_data[0],'k',linewidth=1.5)
    ax1.set_xlabel(xlabel,fontsize=fontsizes[0])
    ax1.set_ylabel(ylabel,fontsize=fontsizes[0])
    ax1.annotate(annotations[0],xy=(0,0.5),xytext=(-ax1.yaxis.labelpad-15,0),
    xycoords=ax1.yaxis.label,textcoords='offset points',size=fontsizes[1],ha='right',va='center')
    ax2.plot(time,meg_data[1],'k',linewidth=1.5)
    ax2.set_xlabel(xlabel,fontsize=fontsizes[0])
    ax2.set_ylabel(ylabel,fontsize=fontsizes[0])
    ax3.plot(time,meg_data[2],'k',linewidth=1.5)
    ax3.set_xlabel(xlabel,fontsize=fontsizes[0])
    ax3.set_ylabel(ylabel,fontsize=fontsizes[0])
    ax3.annotate(annotations[1],xy=(0,0.5),xytext=(-ax3.yaxis.labelpad-15,0),
    xycoords=ax3.yaxis.label,textcoords='offset points',size=fontsizes[1],ha='right',va='center')
    ax4.plot(time,meg_data[3],'k',linewidth=1.5)
    ax4.set_xlabel(xlabel,fontsize=fontsizes[0])
    ax4.set_ylabel(ylabel,fontsize=fontsizes[0])
    ax5.plot(time,meg_data[4],'k',linewidth=1.5)
    ax5.set_xlabel(xlabel,fontsize=fontsizes[0])
    ax5.set_ylabel(ylabel,fontsize=fontsizes[0])
    ax5.annotate(annotations[2],xy=(0,0.5),xytext=(-ax5.yaxis.labelpad-15,0),
    xycoords=ax5.yaxis.label,textcoords='offset points',size=fontsizes[1],ha='right',va='center')
    ax5.annotate(annotations[3],xy=(0.5,0),xytext=(0,-115),xycoords='axes fraction',
    textcoords='offset points',size=fontsizes[1],ha='center',va='bottom')
    ax6.plot(time,meg_data[5],'k',linewidth=1.5)
    ax6.set_xlabel(xlabel,fontsize=fontsizes[0])
    ax6.set_ylabel(ylabel,fontsize=fontsizes[0])
    ax6.annotate(annotations[4],xy=(0.5,0),xytext=(0,-115),xycoords='axes fraction',
    textcoords='offset points',size=35,ha='center',va='bottom')

    plt.setp(ax1.get_xticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax2.get_xticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax3.get_xticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax4.get_xticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax5.get_xticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax6.get_xticklabels(),visible=True,fontsize=fontsizes[2])

    plt.setp(ax2.get_yticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax4.get_yticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax6.get_yticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax1.get_yticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax3.get_yticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax5.get_yticklabels(),visible=True,fontsize=fontsizes[2])

    if savefig:
        plt.savefig(filename+'.png',dpi=600)
        plt.savefig(filename+'.eps',dpi=600)



def plot_PSD_summary(freqs,psd_data,xlabel,ylabel,annotations,savefig,filename,fontsizes):
    f2,((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,sharex=True,
    sharey=True,figsize=[25.0,25.0])

    ax1.plot(freqs,psd_data[0],'k',linewidth=2)
    ax1.axis(xmin=0, xmax=55)
    ax1.set_xlabel(xlabel,fontsize=fontsizes[0])
    ax1.set_ylabel(ylabel,fontsize=fontsizes[0])
    ax1.annotate(annotations[0],xy=(0,0.5),xytext=(-ax1.yaxis.labelpad-15,0),
    xycoords=ax1.yaxis.label,textcoords='offset points',size=35,ha='right',va='center')
    xticks=ax1.xaxis.get_ticklabels()
    xticks[2].set_weight('bold')
    xticks[3].set_weight('bold')
    xticks[4].set_weight('bold')

    ax2.plot(freqs,psd_data[1],'k',linewidth=2)
    ax2.set_xlabel(xlabel,fontsize=fontsizes[0])
    ax2.set_ylabel(ylabel,fontsize=fontsizes[0])
    xticks=ax2.xaxis.get_ticklabels()
    xticks[2].set_weight('bold')
    xticks[3].set_weight('bold')
    xticks[4].set_weight('bold')

    ax3.plot(freqs,psd_data[2],'k',linewidth=2)
    ax3.set_xlabel(xlabel,fontsize=fontsizes[0])
    ax3.set_ylabel(ylabel,fontsize=fontsizes[0])
    ax3.annotate(annotations[1],xy=(0,0.5),xytext=(-ax3.yaxis.labelpad-15,0),
    xycoords=ax3.yaxis.label,textcoords='offset points',size=fontsizes[1],ha='right',va='center')
    xticks=ax3.xaxis.get_ticklabels()
    xticks[2].set_weight('bold')
    xticks[3].set_weight('bold')
    xticks[4].set_weight('bold')

    ax4.plot(freqs,psd_data[3],'k',linewidth=2)
    ax4.set_xlabel(xlabel,fontsize=fontsizes[0])
    ax4.set_ylabel(ylabel,fontsize=fontsizes[0])
    xticks=ax4.xaxis.get_ticklabels()
    xticks[2].set_weight('bold')
    xticks[3].set_weight('bold')
    xticks[4].set_weight('bold')

    ax5.plot(freqs,psd_data[4],'k',linewidth=2)
    ax5.set_xlabel(xlabel,fontsize=fontsizes[0])
    ax5.set_ylabel(ylabel,fontsize=fontsizes[0])
    ax5.annotate(annotations[2],xy=(0,0.5),xytext=(-ax5.yaxis.labelpad-15,0),
    xycoords=ax5.yaxis.label,textcoords='offset points',size=fontsizes[1],ha='right',va='center')
    ax5.annotate(annotations[3],xy=(0.5,0),xytext=(0,-115),xycoords='axes fraction',
    textcoords='offset points',size=fontsizes[1],ha='center',va='bottom')
    xticks=ax5.xaxis.get_ticklabels()
    xticks[2].set_weight('bold')
    xticks[3].set_weight('bold')
    xticks[4].set_weight('bold')

    ax6.plot(freqs,psd_data[5],'k',linewidth=2)
    ax6.set_xlabel(xlabel,fontsize=fontsizes[0])
    ax6.set_ylabel(ylabel,fontsize=fontsizes[0])
    ax6.annotate(annotations[4],xy=(0.5,0),xytext=(0,-115),xycoords='axes fraction',
    textcoords='offset points',size=fontsizes[1],ha='center',va='bottom')
    xticks=ax6.xaxis.get_ticklabels()
    xticks[2].set_weight('bold')
    xticks[3].set_weight('bold')
    xticks[4].set_weight('bold')


    plt.setp(ax1.get_xticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax2.get_xticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax3.get_xticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax4.get_xticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax5.get_xticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax6.get_xticklabels(),visible=True,fontsize=fontsizes[2])

    plt.setp(ax2.get_yticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax4.get_yticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax6.get_yticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax1.get_yticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax3.get_yticklabels(),visible=True,fontsize=fontsizes[2])
    plt.setp(ax5.get_yticklabels(),visible=True,fontsize=fontsizes[2])

    if savefig:
        plt.savefig(filename+'.png',dpi=600)
        plt.savefig(filename+'.eps',dpi=600)

def plot_single_trial(sim_time,spike_times,freqs,data_psd,time,data_meg,xlabels,ylabels,annotations,savefig,filename,fontsizes):
    plt.figure(3,figsize=[20.0,15.0])
    ax1 = plt.subplot2grid((2,2),(0,0), rowspan=2)
    ax2 = plt.subplot2grid((2,2),(0,1))
    ax3 = plt.subplot2grid((2,2),(1,1))

    # Raster plot
    for i,times in enumerate(spike_times[0]):
        y = [i+10]*len(times)
        ax2.plot(times,y,linestyle='None',color='k',marker='|',markersize=10)


    for i,times in enumerate(spike_times[1]):
        y = [i]*len(times)
        ax2.plot(times,y,linestyle='None',color='b',marker='|',markersize=10)
        ax2.axis([0,sim_time,-0.5,30])
    ax2.set_xlabel(xlabels[0],fontsize=fontsizes[0])
    ax2.set_yticks([])
    ax2.annotate(annotations[0],xy=(0,0.5),xytext=(-ax2.yaxis.labelpad,-75),
    xycoords=ax2.yaxis.label,textcoords='offset points',size=fontsizes[1],ha='right',va='top')
    ax2.annotate(annotations[1],xy=(0,0.5),xytext=(-ax2.yaxis.labelpad,75),
    xycoords=ax2.yaxis.label,textcoords='offset points',size=fontsizes[1],ha='right',va='bottom')
    plt.setp(ax2.get_xticklabels(),fontsize=fontsizes[2])

    # MEG plot
    ax3.plot(time,data_meg,'k',linewidth=1.5)
    ax3.set_xlabel(xlabels[0],fontsize=fontsizes[0])
    ax3.set_ylabel(ylabels[0],fontsize=fontsizes[0])
    plt.setp(ax3.get_xticklabels(),fontsize=fontsizes[2])
    plt.setp(ax3.get_yticklabels(),fontsize=fontsizes[2])

    # Power spectrum plot
    ax1.plot(freqs*1000,data_psd,'k',linewidth=2)
    ax1.axis(xmin=0, xmax=55)
    ax1.set_xlabel(xlabels[1],fontsize=fontsizes[0])
    ax1.set_ylabel(ylabels[1],fontsize=fontsizes[0])
    plt.setp(ax1.get_xticklabels(),fontsize=fontsizes[2])
    plt.setp(ax1.get_yticklabels(),fontsize=fontsizes[2])

    if savefig:
        plt.savefig(filename+'.png',dpi=600)
        plt.savefig(filename+'.eps',dpi=600)

