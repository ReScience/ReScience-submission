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
# Generates Fig. 7a: Trading saccade velocity for duration
# (figure 10 in original publication).
#
# The weight W was set equal to 2. For the high velocity trial (dashed line),
# stimulation frequency F was 3 and stimulation duration was 62 ms.
# For the low velocity trial (solid line), F was 1.3 and stimulation duration
# was 117 ms.
#
# The input curve is computed as described in the original publication.
# However, its shape does not match the one shown in the original publication.
# This leads to undesired second saccade.
# -----------------------------------------------------------------------------

import sys
import pylab as pl
from matplotlib import rcParams
from itertools import cycle

# set up model
from setup_model import *

# define additional variables
cm2inch = .394  # inch/cm

# figure setup
rcParams.update({'figure.autolayout': True})

fig_name = '../article/figures/fig7a.eps'
fig_size = np.multiply([11.6, 17.6], cm2inch)
fig_rows = 5
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
    ax[i].set_ylim([-.5, 2.])
    ax[i].set_ylabel([], fontsize=14, rotation=0, labelpad=20)
    ax[i].tick_params(right='off', top='off', labelsize=12)

lines = ['k--', 'k-']
linecycler = cycle(lines)

'''set up experiment

input & external electric stimulation
-------------------------------------
F : stimulation intensity (applied to SC)
W : connection weight from SC to right LLBN

time protocol (in ms)
-----------------------
t_relax  : relaxation period (prior to experiment)
preStim  : time interval prior to stimulus presentation
Stim     : time interval for which stimulus is presented
postStim : time interval after stimulus presentation

time vector T
-------------
since neuron activity is plotted as a function of time,
time vector T is pre-calculated based on time protocol
t_start : start of experiment
t_end   : end of experiment
t_steps : number of simulated time steps
'''
F = [3., 1.3]
W = 2.
nest.Connect(gS, LLBN[0], 'all_to_all', {
             'model': 'rate_connection_instantaneous', 'weight': W})

t_relax = 100
preStim = 50
Stim = [68, 117]
postStim = [132, 83]

t_start = 0
t_end = preStim + Stim[0] + postStim[0]
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
    nest.Connect(MM[s], gS, syn_spec={'delay': dt})
    nest.Connect(MM[s], OPN, syn_spec={'delay': dt})
    nest.Connect(MM[s], LLBN[0], syn_spec={'delay': dt})
    nest.Connect(MM[s], EBN[0], syn_spec={'delay': dt})
    nest.Connect(MM[s], TN[0], syn_spec={'delay': dt})

    # simulate pre-stimulus period
    nest.SetStatus(SC, {'mean': 0.})
    nest.Simulate(preStim)

    # simulate stimulus period
    nest.SetStatus(SC, {'mean': F[s]})
    nest.Simulate(Stim[s])

    # simulate post-stimulus period
    nest.SetStatus(SC, {'mean': 0.})
    nest.Simulate(postStim[s])

    # reset rates for next simulation
    nest.SetStatus(SC, {'rate': 0.})
    nest.SetStatus(TN[0], {'rate': .5})
    nest.SetStatus(TN[1], {'rate': .5})

    # gather data from recording device
    data = nest.GetStatus(MM[s])
    senders = data[0]['events']['senders']
    voltages = data[0]['events']['rate']

    line = next(linecycler)

    # plot
    ax[0].plot(T, voltages[np.where(senders == gS)], line)
    ax[0].set_ylabel('Input')
    ax[1].plot(T, voltages[np.where(senders == LLBN[0])], line)
    ax[1].set_ylabel('LLBN')

    ax[2].plot(T, voltages[np.where(senders == EBN[0])], line)
    ax[2].set_ylabel('EBN')

    ax[3].plot(T, voltages[np.where(senders == OPN)], line)
    ax[3].set_ylabel('OPN')

    ax[4].plot(T, voltages[np.where(senders == TN[0])], line)
    ax[4].get_xaxis().set_visible(True)
    ax[4].set_xlabel('time (ms)', fontsize=14)
    ax[4].set_ylabel('TN')
    ax[4].set_ylim([.3, .7])

ax[0].text(-0.075, 1.1, 'A', transform=ax[0].transAxes, size=16, weight='bold')

# save and display figure
pl.savefig(fig_name, format='eps', dpi=ppi, bbox_inches='tight')
pl.show()
