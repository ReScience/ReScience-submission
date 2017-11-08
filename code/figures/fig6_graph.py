# ----------------------------------------------------------------------------
# Contributors: Renan O. Shimoura
#               Nilton L. Kamiji
#               Rodrigo F. O. Pena
#               Vinicius L. Cordeiro
#               Cesar C. Ceballos
#               Cecilia Romaro
#               Antonio C. Roque
# ----------------------------------------------------------------------------
# References:
#
# *The Cell-Type Specific Cortical Microcircuit: Relating Structure and Activity
# in a Full-Scale Spiking Network Model*,
# Tobias C. Potjans and Markus Diesmann,
# Cerebral Cortex, 24(3):785-806, 2014.
# ----------------------------------------------------------------------------
# File description:
#
# Code to construct figure 6: A) firing rate vs background rate; B) Colorgraph
# with AIness%; C) firing rate vs relative inhibitory synaptic strength
# ----------------------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
import numpy as np

# loading data
data = pd.read_csv('ainess.dat', header = None, sep=' ')
freq_g4 = pd.read_csv('freq_g4.csv', index_col=0); freq_g4 = freq_g4.T
freq_bg8 = pd.read_csv('freq_bg8.csv', index_col=0); freq_bg8 = freq_bg8.T

g = np.arange(2,11)
bg = np.arange(3.0,15.5,0.5)

# figure for the graphs
fig = plt.figure(figsize=(8,8))
gs = gridspec.GridSpec(2, 3, width_ratios=[1, 3, 0.2], height_ratios=[3, 1])

# fig. 6A
ax1 = plt.subplot(gs[0,0])
plt.plot(freq_g4.L23e,bg)
plt.plot(freq_g4.L4e,bg)
plt.plot(freq_g4.L5e,bg)
plt.plot(freq_g4.L6e,bg)
plt.ylim([3,15])
plt.gca().invert_xaxis()

# fig. 6C
ax2 = plt.subplot(gs[1,1])
freq_bg8.plot(y=['L23e','L4e','L5e','L6e'], ax=ax2)
plt.gca().invert_yaxis()

# fig. 6B
ax3 = plt.subplot(gs[0,1])
levels1 = np.linspace(0, 100, num=5, endpoint='True')
levels2 = np.linspace(0, 100, num=25, endpoint='True')
CS1 = plt.contour(g, bg, data, levels=levels1, colors='k', linestyles='dashed')
plt.clabel(CS1, inline=1, fontsize=10, fmt='%1.0f')
CS2 = plt.contourf(g, bg, data, levels=levels2, cmap=cm.jet)
plt.title('AIness')
plt.xlabel('relative inh. synaptic strength')
plt.ylabel('background rate [Hz]')

# colorbar
ax4 = plt.subplot(gs[0,2])
plt.colorbar(CS2, cax=ax4)

plt.tight_layout()
plt.savefig('ainess.jpg', dpi=600)
