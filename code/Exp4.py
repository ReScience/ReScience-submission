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
# Generates Fig. 4: Saccadic staircase
# (figure 7 in original publication).
# Inputs to the hizontal and vical circuits were held at (0.2, 0.33) for 250 ms
# -----------------------------------------------------------------------------

import sys
import pylab as pl
from matplotlib import rcParams
from scipy.signal import argrelextrema


# set up model
from setup_model import *

# define additional variables
cm2inch = .394  # inch/cm
g_pos = 260.    # gain eye position
sr = 40   		# sampling rate

# figure setup
rcParams.update({'figure.autolayout': True})

fig_name = '../article/figures/fig4.eps'
fig_size = np.multiply([8.5, 11.6], cm2inch)
ppi = 1200
face = 'white'
edge = 'white'
fig = pl.figure(facecolor=face, edgecolor=edge, figsize=fig_size)
ax = fig.add_subplot(1, 1, 1)

ax.set_xlim([0., 10.])
ax.set_ylim([0., 15.5])
ax.tick_params(labelsize=12)
ax.set_xlabel('hizontal eye position (deg)', fontsize=14)
ax.set_ylabel('vical eye position (deg)', fontsize=14)

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
I_h = .20
I_v = .33

t_relax = 100
preStim = 0
Stim = 250
postStim = 50

'''set up recording device

neuron activity is recorded using NEST's multimeter object
recording the 'rate' property of a neuron after each time step
'''
multimeter = nest.Create('multimeter')
nest.SetStatus(multimeter, {'interval': dt, 'record_from': ['rate']})

'''relaxation

let the system reach equilibrium
in the absence of input and stimulation.
'''
nest.Simulate(t_relax)

'''connect recording device

multimeter is connected to neurons of interest
'''
nest.Connect(multimeter, TN[1], syn_spec={'delay': dt})
nest.Connect(multimeter, TN[3], syn_spec={'delay': dt})

# simulate pre-stimulus period
nest.Simulate(preStim)

# simulate stimulus period
nest.SetStatus(LLBN[1], {'mean': I_h})
nest.SetStatus(LLBN[3], {'mean': I_v})
nest.Simulate(Stim)

# simulate post-stimulus period
nest.SetStatus(LLBN[1], {'mean': 0.})
nest.SetStatus(LLBN[3], {'mean': 0.})
nest.Simulate(postStim)

# gather data from recording device
data = nest.GetStatus(multimeter)
senders = data[0]['events']['senders']
voltage = data[0]['events']['rate']

# compute output variables (hizontal and vical eye position)
theta_h = g_pos * (voltage[np.where(senders == TN[1])] - .5)
theta_v = g_pos * (voltage[np.where(senders == TN[3])] - .5)

# plot
ax.plot(theta_h[::sr], theta_v[::sr], 'k.')

# save and display figure
pl.savefig(fig_name, format='eps', dpi=ppi, bbox_inches='tight')
pl.show()
