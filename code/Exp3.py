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
# Generates Fig. 3: Oblique saccades
# (figure 6 in original publication).
# Input to the horizontal (right) and vertical (up) circuits were:
# (0.67, 0.08), (0.7, 0.22), (0.74,0.4), (0.75, 0.6), (0.7, 0.9).
# These inputs were left on for 75 ms.
# -----------------------------------------------------------------------------

import sys
import pylab as pl
from matplotlib import rcParams


# set up model
from setup_model import *

# define additional variables
cm2inch = .394  # inch/cm
g_pos = 260.    # gain of eye position

# figure setup
rcParams.update({'figure.autolayout': True})

fig_name = '../article/figures/fig3.eps'
fig_size = np.multiply([11.6, 11.6], cm2inch)
ppi = 1200
face = 'white'
edge = 'white'
fig = pl.figure(facecolor=face, edgecolor=edge, figsize=fig_size)
ax = fig.add_subplot(1, 1, 1)

ax.axhline(color='k', linestyle='--')
ax.axvline(color='k', linestyle='--')
ax.set_xlim([-20., 20.])
ax.set_ylim([-20., 20.])
ax.tick_params(labelsize=12)
ax.locator_params(axis='x', nbins=5)
ax.locator_params(axis='y', nbins=5)
ax.set_xlabel('horizontal eye position (deg)', fontsize=14)
ax.set_ylabel('vertical eye position (deg)', fontsize=14)

'''set up experiment

input & external electric stimulation
-------------------------------------
I_h : input applied to right LLBN
I_v : input applied to up LLBN

time protocol (in ms)
-----------------------
t_relax  : relaxation period (prior to experiment)
preStim  : time interval prior to stimulus presentation
Stim     : time interval for which stimulus is presented
postStim : time interval after stimulus presentation
'''
I_h = [.67, .70, .74, .75, .70]
I_v = [.08, .22, .40, .60, .90]

t_relax = 100
preStim = 0
Stim = 75
postStim = 0

'''set up recording devices

neuron activity is recorded using NEST's multimeter object
recording the 'rate' property of a neuron after each time step
'''
MM = [None] * 5
for s in range(0, 5):
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
    nest.Connect(MM[s], TN[1], syn_spec={'delay': dt})
    nest.Connect(MM[s], TN[3], syn_spec={'delay': dt})

    # simulate pre-stimulus period
    nest.SetStatus(TN[1], {'rate': .5})
    nest.SetStatus(TN[3], {'rate': .5})
    nest.Simulate(preStim)

    # simulate stimulus period
    nest.SetStatus(LLBN[1], {'mean': I_h[s]})
    nest.SetStatus(LLBN[3], {'mean': I_v[s]})
    nest.Simulate(Stim)

    # simulate post-stimulus period
    nest.SetStatus(LLBN[1], {'mean': 0.})
    nest.SetStatus(LLBN[3], {'mean': 0.})
    nest.Simulate(postStim)

    # gather data from recording device
    data = nest.GetStatus(MM[s])
    senders = data[0]['events']['senders']
    voltage = data[0]['events']['rate']

    # compute output variables (horizontal and vertical eye position)
    theta_h = g_pos * (voltage[np.where(senders == TN[1])] - .5)
    theta_v = g_pos * (voltage[np.where(senders == TN[3])] - .5)

    # plot
    ax.plot(theta_h, theta_v, linewidth=2)

# save and display figure
pl.savefig(fig_name, format='eps', dpi=ppi, bbox_inches='tight')
pl.show()
