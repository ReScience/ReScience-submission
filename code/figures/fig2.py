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
# Code to construct figure 2: A) Raster plot; B) box plot of firing rates;
# C) barplot of CV's; D) barplot of synchrony index.
# ----------------------------------------------------------------------------

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

###############################################################################
# Loading data and defining parameters
###############################################################################
data = pd.read_csv('data_raster_g4.0_w87.8.dat', sep=" ",
                    header=None, names=['i','t'])

tsim = 10.0     # time of simulation in seconds

# cortical layer labels: e for excitatory; i for inhibitory
lname = ['L23e', 'L23i', 'L4e', 'L4i', 'L5e', 'L5i','L6e', 'L6i']

# number of neurons by layer
n_layer = [0, 20683, 5834, 21915, 5479, 4850, 1065, 14395, 2948];
l_bins = np.cumsum(n_layer) # cumulative number of neurons by layer

# graphs color codes: different colors for different layers
dotsize = 2.5
dotcolor = np.array([[0.0, 0.0, 255.0],
                    [102.0, 178.0, 255.0],
                    [255.0, 128.0, 0.0],
                    [255.0, 178.0, 102.0],
                    [0.0,   128.0, 0.0],
                    [153.0, 255.0, 153.0],
                    [255.0, 0.0,   0.0],
                    [255.0, 153.0, 153.0]])/255.0

# grouping spiking times for each neuron
keys,values = data.sort_values(['i','t']).values.T
ukeys,index=np.unique(keys,True)
arrays=np.split(values,index[1:])
spk_neuron = pd.DataFrame({'i':ukeys,'t':[list(a) for a in arrays]})

# creating a flag to identify cortical layers
spk_neuron['layer'] = pd.cut(spk_neuron['i'], l_bins, labels=lname, right=False)
data['layer'] = pd.cut(data['i'], l_bins, labels=lname, right=False)

# sampling data:
psample = 0.025 # percentage of neurons by layer for the raster plot
n_sample = 1000 # number of neurons by layer for sampled measures
spk_neuron = spk_neuron.groupby(['layer']).apply(lambda x: x.sample(n=n_sample))

# measures DataFrame:
measures_layer = pd.DataFrame(index=lname)

# cleaning variables
del keys, values, ukeys, index, arrays

# figure size for the graphs
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(7.5,10))

###############################################################################
# Raster plot
###############################################################################
plt.subplot2grid((3,2),(0,0), rowspan=3)
plt.gca().set_yticklabels([])
acum_index = 0

for i in range(len(lname)):
    index_start = l_bins[i]
    index_end = l_bins[i]+int(psample*n_layer[i+1])

    x = data.t[data.i.isin(range(index_start,index_end))]
    y = data.i[data.i.isin(range(index_start,index_end))] + acum_index - index_start

    plt.plot(x/1000.0,y,'.',markersize=dotsize,color=dotcolor[i])

    # layers labels
    xpos = tsim-440/1000.0
    ypos = acum_index + (index_end-index_start)/2.0
    plt.text(xpos,ypos,lname[i],horizontalalignment='center', fontweight='bold')

    acum_index = acum_index + (index_end-index_start)

plt.xlim(tsim-400/1000.0,tsim)
plt.ylim(0,acum_index)
plt.xlabel('time [s]')
plt.ylabel(' ')
plt.gca().invert_yaxis()

###############################################################################
# Firing rates
###############################################################################
freq = []
freq = [float(len(spk_neuron.t[i]))/tsim for i in range(len(spk_neuron))]
spk_neuron['f'] = freq
measures_layer['f'] = spk_neuron.groupby(['layer'])['f'].mean()

# boxplot of firing rates by layer
bplot = spk_neuron.boxplot(column = 'f', by = 'layer', showmeans=True,
                    vert = False, rot = 30, ax = axes[0,1],
                    patch_artist=True, sym='+', return_type='dict', grid=False)

[bplot[0]['boxes'][i].set_facecolor(dotcolor[i]) for i in range(0,len(bplot[0]['boxes']))]
[bplot[0]['means'][i].set_markerfacecolor('white') for i in range(0,len(bplot[0]['boxes']))]
[bplot[0]['means'][i].set_markeredgecolor('k') for i in range(0,len(bplot[0]['boxes']))]

axes[0,1].set_title("")
axes[0,1].set_ylabel("")
axes[0,1].set_xlabel('firing rates[Hz]')
axes[0,1].invert_yaxis()
fig.suptitle("")

###############################################################################
# Interspike intervals + coefficient of variation
###############################################################################
# interspike intervals
isi = []
isi = [np.diff(spk_neuron.t[i]) for i in range(len(spk_neuron))]

# coefficient of variation
cv = []
cv = [np.std(isi[i])/np.mean(isi[i]) if len(isi[i])>1 else np.nan\
        for i in range(len(spk_neuron))]
spk_neuron['cv'] = cv

measures_layer['cv'] = spk_neuron.groupby(['layer'])['cv'].mean()

# barplot of mean CV
plt.subplot2grid((3,2),(1,1))
measures_layer['cv'].plot.barh(edgecolor='k' ,color=dotcolor, rot=30, width=0.8)
plt.ylabel("")
plt.xlabel('irregularity')
plt.gca().invert_yaxis()

###############################################################################
# Synchrony index
###############################################################################
sync = []
bins = np.arange(0,tsim*1000.0+3.0,3)

for i in range(len(lname)):
    index_sample = spk_neuron.i[spk_neuron.layer.isin([lname[i]])]
    count, division = np.histogram(data.t[data.i.isin(index_sample)],bins=bins)
    sync.append(np.var(count[166:])/np.mean(count[166:]))

measures_layer['sync'] = sync

# barplot of synchrony index
y_pos = np.arange(len(lname))
plt.subplot2grid((3,2),(2,1))
measures_layer['sync'].plot.barh(edgecolor='k' ,color=dotcolor, rot=30, width=0.8)
plt.ylabel("")
plt.xlabel('synchrony')
plt.gca().invert_yaxis()

plt.tight_layout()
plt.subplots_adjust(left=0.07)
plt.savefig('fig6.jpg', dpi=600)
