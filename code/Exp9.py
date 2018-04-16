# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Mario Senden mario.senden@maastrichtuniversity.nl
# -----------------------------------------------------------------------------
# References:
#
# Gancarz, G., Grossberg, S. "A Neural Model of the Saccade Generator in the
# Reticular Formation." Neural Networks 11, no. 7-8 (October 1998): 1159-74.
# doi:10.1016/S0893-6080(98)00096-3.
# -----------------------------------------------------------------------------
# File description:
#
# Simulates reimplemented saccade generation (SG) model of
# Gancarz & Grossberg (1998).
# Generates Fig. 9: Interrupted saccade resulting from OPN stimulation
# (figure 12 in original publication).
# Input to left side of saccade generator was set to 0.7 for 100 ms.
# OPN stimulation J was set to 1.8 for 5 ms starting 45 ms after input onset.
# -----------------------------------------------------------------------------

import sys
import pylab as pl
from matplotlib import rcParams

# set up model
from setup_model import *

# define additional variables
cm2inch = .394  # inch/cm

# figure setup
rcParams.update({'figure.autolayout': True})

fig_name = '../article/figures/fig9.eps'
fig_size = np.multiply([11.6, 17.6], cm2inch)
fig_rows = 6
fig_cols = 1
fig_plots = fig_rows * fig_cols
ppi = 1200
face = 'white'
edge = 'white'

ax = [None] * fig_plots
fig = pl.figure(facecolor=face, edgecolor=edge, figsize=fig_size)
for i in range(0, fig_plots):
    ax[i] = fig.add_subplot(fig_rows, fig_cols, i + 1)
    ax[i].get_xaxis().set_visible(False)
    ax[i].get_yaxis().set_ticks([])
    ax[i].set_xlim([0, 200])
    ax[i].set_ylim([-1., 1.5])
    ax[i].set_ylabel([], fontsize=14, rotation=0, labelpad=20)
    ax[i].tick_params(right='off', top='off', labelsize=12)

'''set up experiment

input & external electric stimulation
-------------------------------------
I : input applied to left LLBN
J : external electric stimulation of OPN

time protocol (in ms)
-----------------------
t_relax  : relaxation period (prior to experiment)
preStim  : time interval prior to stimulus presentation
Stim_1   : stimulus is presented; OPN is not stimulated
Stim_2   : stimulus is presented; OPN is stimulated
Stim_3   : stimulus is presented; OPN is not stimulated
postStim : time interval after stimulus presentation

time vector T
-------------
since neuron activity is plotted as a function of time,
time vector T is pre-calculated based on time protocol
t_start : start of experiment
t_end   : end of experiment
t_steps : number of simulated time steps
'''
I = .7
J = [1.8, 0.]

t_relax = 100
preStim = 50
Stim_1 = 45
Stim_2 = 5
Stim_3 = 50
postStim = 50

t_start = 0
t_end = preStim + Stim_1 + Stim_2 + Stim_3 + postStim
t_steps = int((t_end - t_start) / dt) - 1
T = np.linspace(t_start, t_end, t_steps)

'''set up recording devices

neuron activity is recorded using NEST's multimeter object
recording the 'rate' property of a neuron after each time step
'''
MM = [None] * 2
for s in range(0, 2):
    MM[s] = nest.Create('multimeter')
    nest.SetStatus(MM[s], {'interval': dt, 'record_from': ['rate']})

    '''relaxation

    let the system reach equilibrium
    in the absence of input and stimulation.
    '''
    nest.Simulate(t_relax)

    '''connect recording devices

    multimeters are connected to neurons of interest
    '''
    nest.Connect(MM[s], OPN, syn_spec={'delay': dt})
    nest.Connect(MM[s], LLBN[0], syn_spec={'delay': dt})
    nest.Connect(MM[s], EBN[0], syn_spec={'delay': dt})
    nest.Connect(MM[s], TN[0], syn_spec={'delay': dt})

    # simulate pre-stimulus period
    nest.SetStatus(TN[0], {'rate': .5})
    nest.SetStatus(TN[1], {'rate': .5})
    nest.Simulate(preStim)

    # simulate stimulus period 1
    nest.SetStatus(LLBN[0], {'mean': I})
    nest.Simulate(Stim_1)

    # simulate stimulus period 2
    # stimulate OPN
    nest.SetStatus(EXT, {'mean': J[s]})
    nest.Simulate(Stim_2)

    # simulate stimulus period 3
    # cease OPN stimulation
    nest.SetStatus(EXT, {'mean': 0.})
    nest.Simulate(Stim_3)

    # simulate post-stimulus period
    nest.SetStatus(LLBN[0], {'mean': 0.})
    nest.Simulate(postStim)

    # gather data from recording device
    data = nest.GetStatus(MM[s])
    senders = data[0]['events']['senders']
    voltages = data[0]['events']['rate']

    # plot
    if (s < 1):
        Input = np.zeros(t_end - t_start)
        Input[preStim + 1:t_end - postStim] = I
        Ext = np.zeros(t_end - t_start)
        Ext[preStim + Stim_1 + 1:preStim + Stim_1 + Stim_2] = J[s]

        ax[0].plot(range(t_start, t_end), Input, 'k')
        ax[0].set_ylabel('Input')

        ax[1].plot(T, voltages[np.where(senders == LLBN[0])], 'k')
        ax[1].set_ylabel('LLBN')

        ax[2].plot(T, voltages[np.where(senders == EBN[0])], 'k')
        ax[2].set_ylabel('EBN')

        ax[3].plot(T, voltages[np.where(senders == OPN)], 'k')
        ax[3].set_ylabel('OPN')

        ax[4].plot(T, voltages[np.where(senders == TN[0])], 'k')
        ax[4].set_ylim([.49, .575])
        ax[4].set_ylabel('TN')

        ax[5].plot(range(t_start, t_end), Ext, 'k')
        ax[5].get_xaxis().set_visible(True)
        ax[5].set_xlabel('time (ms)', fontsize=14)
        ax[5].set_ylabel('OPNSTIM', labelpad=40)
        ax[5].set_ylim([-.1, 2.5])
    else:
        ax[4].plot(T, voltages[np.where(senders == TN[0])], 'k--')

# save and display figure
pl.savefig(fig_name, format='eps', dpi=ppi, bbox_inches='tight')
pl.show()
