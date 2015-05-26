# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) 2014, Nicolas P. Rougier
# Distributed under the (new) BSD License.
#
# Contributors: Nicolas P. Rougier (Nicolas.Rougier@inria.fr)
# -----------------------------------------------------------------------------
# References:
#
# * Interaction between cognitive and motor cortico-basal ganglia loops during
#   decision making: a computational study. M. Guthrie, A. Leblois, A. Garenne,
#   and T. Boraud. Journal of Neurophysiology, 109:3025–3040, 2013.
# -----------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

# Record type for display
htype = [ ("CTX", [("mot", float, 4), ("cog", float, 4), ("ass", float, 16)]),
          ("STR", [("mot", float, 4), ("cog", float, 4), ("ass", float, 16)]),
          ("GPI", [("mot", float, 4), ("cog", float, 4)]),
          ("THL", [("mot", float, 4), ("cog", float, 4)]),
          ("STN", [("mot", float, 4), ("cog", float, 4)])]


# Cortex activity display
# -----------------------
def display_ctx(history, duration=3.0, filename=None):
    fig = plt.figure(figsize=(12,5))
    plt.subplots_adjust(bottom=0.15)

    timesteps = np.linspace(0,duration, len(history))

    fig.patch.set_facecolor('.9')
    ax = plt.subplot(1,1,1)

    plt.plot(timesteps, history["CTX"]["cog"][:,0],c='r', label="Cognitive Cortex")
    plt.plot(timesteps, history["CTX"]["cog"][:,1],c='r')
    plt.plot(timesteps, history["CTX"]["cog"][:,2],c='r')
    plt.plot(timesteps, history["CTX"]["cog"][:,3],c='r')
    plt.plot(timesteps, history["CTX"]["mot"][:,0],c='b', label="Motor Cortex")
    plt.plot(timesteps, history["CTX"]["mot"][:,1],c='b')
    plt.plot(timesteps, history["CTX"]["mot"][:,2],c='b')
    plt.plot(timesteps, history["CTX"]["mot"][:,3],c='b')

    plt.xlabel("Time (seconds)")
    plt.ylabel("Activity (Hz)")
    plt.legend(frameon=False, loc='upper left')
    plt.xlim(0.0,duration)
    plt.ylim(0.0,60.0)

    plt.xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
               ['0.0','0.5\n(Trial start)','1.0','1.5', '2.0','2.5\n(Trial stop)','3.0'])

    if filename is not None:
        plt.savefig(filename)
    plt.show()



# All but cortex activity display
# -------------------------------
def display_all(history, duration=3.0, filename=None):
    fig = plt.figure(figsize=(18,12))
    fig.patch.set_facecolor('1.0')

    timesteps = np.linspace(0,duration, len(history))

    def subplot(rows,cols,n, alpha=0.0):
        ax = plt.subplot(rows,cols,n)
        ax.patch.set_facecolor("k")
        ax.patch.set_alpha(alpha)

        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.yaxis.set_ticks_position('left')
        ax.yaxis.set_tick_params(direction="outward")
        return ax

    ax = subplot(5,3,1)
    ax.set_title("Motor", fontsize=24)
    ax.set_ylabel("STN", fontsize=24)
    for i in range(4):
        plt.plot(timesteps, history["STN"]["mot"][:,i], c='k', lw=.5)
    ax.set_xticks([])

    ax = subplot(5,3,2)
    ax.set_title("Cognitive", fontsize=24)
    for i in range(4):
        plt.plot(timesteps, history["STN"]["cog"][:,i], c='k', lw=.5)
    ax.set_xticks([])

    ax = subplot(5,3,3,alpha=0)
    ax.set_title("Associative", fontsize=24)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['left'].set_color('none')


    ax = subplot(5,3,4)
    ax.set_ylabel("Cortex", fontsize=24)
    for i in range(4):
        plt.plot(timesteps, history["CTX"]["mot"][:,i], c='k', lw=.5)
    ax.set_xticks([])

    ax = subplot(5,3,5)
    for i in range(4):
        plt.plot(timesteps, history["CTX"]["cog"][:,i], c='k', lw=.5)
    ax.set_xticks([])

    ax = subplot(5,3,6)
    for i in range(16):
        plt.plot(timesteps, history["CTX"]["ass"][:,i], c='k', lw=.5)
    ax.set_xticks([])

    ax = subplot(5,3,7)
    ax.set_ylabel("Striatum", fontsize=24)
    for i in range(4):
        plt.plot(timesteps, history["STR"]["mot"][:,i], c='k', lw=.5)
    ax.set_xticks([])

    ax = subplot(5,3,8)
    for i in range(4):
        plt.plot(timesteps, history["STR"]["cog"][:,i], c='k', lw=.5)
    ax.set_xticks([])

    ax = subplot(5,3,9)
    for i in range(16):
        plt.plot(timesteps, history["STR"]["ass"][:,i], c='k', lw=.5)
    ax.set_xticks([])

    ax = subplot(5,3,10)
    ax.set_ylabel("GPi", fontsize=24)
    for i in range(4):
        plt.plot(timesteps, history["GPI"]["mot"][:,i], c='k', lw=.5)
    ax.set_xticks([])

    ax = subplot(5,3,11)
    for i in range(4):
        plt.plot(timesteps, history["GPI"]["cog"][:,i], c='k', lw=.5)
    ax.set_xticks([])

    ax = subplot(5,3,13)
    ax.set_ylabel("Thalamus", fontsize=24)
    for i in range(4):
        plt.plot(timesteps, history["THL"]["mot"][:,i], c='k', lw=.5)
    ax.set_xticks([])

    ax = subplot(5,3,14)
    for i in range(4):
        plt.plot(timesteps, history["THL"]["cog"][:,i], c='k', lw=.5)
    ax.set_xticks([])

    if filename is not None:
        plt.savefig(filename)
    plt.show()
