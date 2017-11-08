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
# Code to construct figure 5: A1 and B1 Raster plots of network activity;
# A2 and B2 barplots of firing rates
# ----------------------------------------------------------------------------

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

###############################################################################
# Filenames
###############################################################################
#             DC current                      Layer-independent
filename = (['dc_data_raster_g4.0_w87.8.dat','bg_data_raster_g4.0_w87.8.dat'])

###############################################################################
# Parameters
###############################################################################
tsim = 10.0     # time of simulation in seconds
n = 77168       # total number of neurons

# cortical layer labels: e for excitatory; i for inhibitory
lname = ['L23e', 'L23i', 'L4e', 'L4i', 'L5e', 'L5i','L6e', 'L6i']

# number of neurons by layer
n_layer = [0, 20683, 5834, 21915, 5479, 4850, 1065, 14395, 2948];
l_bins = np.cumsum(n_layer) # cumulative number of neurons by layer

# graphs color codes: different colors for different layers
dotsize = 1.5
dotcolor = np.array([[0.0, 0.0, 255.0],
                    [102.0, 178.0, 255.0],
                    [255.0, 128.0, 0.0],
                    [255.0, 178.0, 102.0],
                    [0.0,   128.0, 0.0],
                    [153.0, 255.0, 153.0],
                    [255.0, 0.0,   0.0],
                    [255.0, 153.0, 153.0]])/255.0

# figure size for the graphs
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(10,10))

for k in range(0,2):

    # loading data
    data = pd.read_csv(filename[k], sep=" ", header=None, names=['i','t'])

    # grouping spiking times for each neuron
    keys,values = data.sort_values(['i','t']).values.T
    ukeys,index=np.unique(keys,True)
    arrays=np.split(values,index[1:])
    spk_neuron = pd.DataFrame({'i':ukeys,'t':[list(a) for a in arrays]})

    # creating a flag to identify cortical layers
    spk_neuron['layer'] = pd.cut(spk_neuron['i'], l_bins, labels=lname, right=False)
    data['layer'] = pd.cut(data['i'], l_bins, labels=lname, right=False)

    # cleaning variables
    del keys, values, ukeys, index, arrays

    # sampling data:
    psample = 0.025 # percentage of neurons by layer for the raster plot
    n_sample = 1000 # number of neurons by layer for sampled measures
    spk_neuron = spk_neuron.groupby(['layer']).apply(lambda x: x.sample(n=1000))

    # measures DataFrame:
    measures_layer = pd.DataFrame(index=lname)

    ############################################################################
    # Raster plot
    ############################################################################
    plt.subplot2grid((2,3),(0,k*2), rowspan=2)
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
        plt.text(xpos,ypos,lname[i],horizontalalignment='center',fontweight='bold')

        acum_index = acum_index + (index_end-index_start)

    plt.xlim(tsim-400/1000.0,tsim) # ploting the last 400 ms of simulation time
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
    plt.subplot2grid((2,3),(k,1))
    measures_layer['f'].plot.barh(edgecolor='k' ,color=dotcolor, rot=30, width=0.8)
    plt.ylabel("")
    plt.xlabel('firing rate')
    plt.gca().invert_yaxis()

    # cleaning variables
    del data, measures_layer, spk_neuron

plt.tight_layout()
plt.subplots_adjust(left=0.05)
plt.savefig('fig7.jpg', dpi=600)
