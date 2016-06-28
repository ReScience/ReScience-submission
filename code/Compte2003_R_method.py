# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Andrei Maksimov (maksimov.andrei7@gmail.com)
# -----------------------------------------------------------------------------
# References:
#
# *Cellular and network mechanisms of slow oscillatory activity (<1 Hz)
# and wave propagations in a cortical network model*, A. Compte,
# M.V. Sanchez-Vives, D.A. McCormick, X.-J. Wang,
# Journal of Neurophysiology, 2707--2725, 2003"
# -----------------------------------------------------------------------------
# File description:
#
# generate a schematic representation of resistance measurement method (Fig. 3)
# -----------------------------------------------------------------------------

import numpy as np
import awesomizer as aw
import matplotlib.pyplot as plt

# numpy seed for reproducibility
np.random.seed(11)

# apply default matplotlib style
style = aw.style_Plos()
aw.style_apply(style)


fig = plt.figure()


# create trace of V_m data

# create noisy data
time = np.arange(0., 1000., 1.)
signal = np.ones(len(time))
signal[time < 200.] = 0.
signal[time > 800.] = 0.
signal += np.random.normal(1., 0.3, len(time)) * 2 - 2

# smooth data
bin_width = 50
signal = [np.mean(signal[i:i + bin_width])
          for i in xrange(len(time) - bin_width)]
signal = np.array(signal)
time = time[0:-bin_width]


# plot V_m trace
ax = fig.add_subplot(111)
ax.plot(time, signal, '0.5')
ax.plot(time, signal, '0.5')

# add point to V_m trace where measurement is performed
t_meas = 505
ax.plot([t_meas], signal[t_meas], 'ko', markersize=15)

# add V0 label
ax.text(t_meas, signal[t_meas] + 0.05, '$V_0$', ha='right', va='bottom', 
        color='k', size='xx-large')


# add scheme of virtual injection as an inset

dc_dur = 200.    # the inset spans this time window. Time window includes phases
                 # of initial steady period, hyperpolarization, and recovery.  
t0 = t_meas - 100.    # x position of the inset drawing left top corner
V0 = signal[t_meas] * 0.6  # y position of the inset drawing left top corner

# add trace for initial steady period
t_start = t0
t_stop = t_start + int(dc_dur * 0.2)
time_temp = time[t_start:t_stop]
trace = V0 * np.ones(len(time_temp))

ax.plot(time_temp, trace, 'k', linewidth=2)

# add trace for hyperpolarizing behavior
t_start = t_stop
t_stop = t0 + int(dc_dur * 0.7)
time_temp = time[t_start:t_stop]
V1 = signal[t_meas] * 0.3
tau = 10.
trace = V1 + (V0 - V1) * np.exp(-(time_temp - t_start) / tau)

ax.plot(time_temp, trace, 'k', linewidth=2)

# add trace for recovery from hyperpolarization
t_start = t_stop
t_stop = t0 + dc_dur
time_temp = time[t_start:t_stop]
trace = V0 - (V0 - V1) * np.exp(-(time_temp - t_start) / tau)

ax.plot(time_temp, trace, 'k', linewidth=2)


# add horizontal lines depicting V0 and V1
ax.hlines([V0, V1], t0 - 20., t_stop + 20., linestyle='dashed', color='0.5')

# place V0 and V1 labels
ax.text(t0, V0, r'$V_0 \rightarrow I_0$', ha='center', va='bottom', color='k',
        size='xx-large')
ax.text(t_start, V1 - 0.05, r'$V_1 \rightarrow I_1$', ha='right', va='top', 
        color='k', size='xx-large')

# add annotations
t_ann = t0 + dc_dur * 0.7 - 10.
ax.annotate('', (t_ann, V1), (t_ann, V0), arrowprops={'arrowstyle': '<->',
                                                      'shrinkA': 0., 
                                                      'shrinkB': 0., 
                                                      'linestyle': 
                                                      'dashed'}, color='0.5')

ax.text(t_ann - 15., (V0 + V1) / 2, '$\Delta V$\n$\Delta I$', ha='right', 
        va='center', size='xx-large')

# add connection lines between two figures
ax.annotate('', (t0, V0 + 0.1), (t_meas, signal[t_meas]), 
            arrowprops={'arrowstyle': '-',
                        'shrinkA': 0.,
                        'shrinkB': 0., 
                        'linestyle': 'dashed'}, color='0.5')
ax.annotate('', (t0 + dc_dur, V0 + 0.1), (t_meas, signal[t_meas]), 
            arrowprops={'arrowstyle': '-',
                        'shrinkA': 0.,
                        'shrinkB': 0.,
                        'linestyle': 'dashed'}, color='0.5')


# finalize figure
ax.set_xlabel('time')
ax.set_ylabel('$V_{\mathrm{m}}$')
ax.set_xticks([0., t_meas, time[-1]])
ax.set_xticklabels(['$0$', '$t_0$', '$t_1$'])
ax.set_yticks([])


axes = [{'ax': ax, 'cb': None}]
aw.tight_layout(axes, style, './fig3.png', fig_width='medium',
                label_order=[''])


plt.show()
