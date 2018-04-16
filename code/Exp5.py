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
# Generates Fig. 5: EBN tuning curve (figure 8 in original publication).
# Inputs to the horizontal and vertical circuits:
# (.2, 0., .63, .0),(0., 0., .7, 0.),(0., .2, .63, 0.),(0., .45, .45, 0.),
# (0., .7, 0., 0.),(0., .45, 0., .45),(0., .2, 0., .63),(0., 0., 0., .7),
# (.2, 0., 0., .63),(.45, 0., 0., .45),(0.63, 0., 0., .2),(0.7, 0., 0., 0.),
# (0.63, 0., .2, 0.),(0.45, 0., .45, 0.).
# The input was applied for 50 ms.
# -----------------------------------------------------------------------------

import sys
import pylab as pl
from matplotlib import rcParams

# set up model
from setup_model import *

# define additional variables
cm2inch = .394 			      # inch/cm
r = np.zeros(15)		      # EBN activity
a = np.zeros(15)		      # saccade angle

# figure setup
rcParams.update({'figure.autolayout': True})

fig_name = '../article/figures/fig5.eps'
fig_size = np.multiply([8.5, 8.5], cm2inch)
ppi = 1200
face = 'white'
edge = 'white'
fig = pl.figure(facecolor=face, edgecolor=edge, figsize=fig_size)
ax = fig.add_subplot(1, 1, 1, projection='polar')
ax.tick_params(labelsize=12)
ax.set_rticks([.02, .1, 1.8])
ax.set_yticklabels([])

'''set up experiment

input & external electric stimulation
-------------------------------------
I : input applied to the LLBN

time protocol (in ms)
-----------------------
t_relax  : relaxation period (prior to experiment)
preStim  : time interval prior to stimulus presentation
Stim     : time interval for which stimulus is presented
postStim : time interval after stimulus presentation
'''
I = [[.2, 0., .63, .0], [0., 0., .7, 0.], [0., .2, .63, 0.],
     [0., .45, .45, 0.], [0., .7, 0., 0.], [0., .45, 0., .45],
     [0., .2, 0., .63], [0., 0., 0., .7], [.2, 0., 0., .63],
     [.45, 0., 0., .45], [0.63, 0., 0., .2], [0.7, 0., 0., 0.],
     [0.63, 0., .2, 0.], [0.45, 0., .45, 0.]]

t_relax = 100
preStim = 0
Stim = 50
postStim = 50

'''set up recording devices

neuron activity is recorded using NEST's multimeter object
recording the 'rate' property of a neuron after each time step
'''
MM = [None] * 14
for s in range(0, 14):
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
    nest.Connect(MM[s], EBN[0], syn_spec={'delay': dt})

    # simulate pre-stimulus period
    nest.Simulate(preStim)

    # simulate stimulus period
    nest.SetStatus(LLBN[0], {'mean': I[s][0]})
    nest.SetStatus(LLBN[1], {'mean': I[s][1]})
    nest.SetStatus(LLBN[2], {'mean': I[s][2]})
    nest.SetStatus(LLBN[3], {'mean': I[s][3]})
    nest.Simulate(Stim)

    # simulate post-stimulus period
    nest.SetStatus(LLBN[0], {'mean': 0.})
    nest.SetStatus(LLBN[1], {'mean': 0.})
    nest.SetStatus(LLBN[2], {'mean': 0.})
    nest.SetStatus(LLBN[3], {'mean': 0.})
    nest.Simulate(postStim)

    # gather data from recording device
    data = nest.GetStatus(MM[s])
    senders = data[0]['events']['senders']
    voltage = data[0]['events']['rate']

    # compute output variables (angle given by stimulus, radius given by EBN
    # response)
    a[s] = np.math.atan2(I[s][3] - I[s][2], I[s][1] - I[s][0])
    r[s] = np.mean(np.maximum(voltage[np.where(senders == EBN[0])], 0.))
r[-1] = r[0]
a[-1] = a[0]

# plot
ax.plot(a, r, '-ko', linewidth=2)

# save and display figure
pl.savefig(fig_name, format='eps', dpi=ppi, bbox_inches='tight')
pl.show()
