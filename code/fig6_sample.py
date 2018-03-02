import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec

import pandas as pd
import numpy as np

# ----------------------------------------------------------------------------
# File description:
#
# Calculation of the %AIness index used to construct figure 6
# ----------------------------------------------------------------------------
tsim = 10.0

g = np.arange(2.0,4.5,0.5);            # relative inh. synaptic strength values
bg_rate = np.arange(3.0,6.5,0.5);  # background rate values

###############################################################################
# Parameters
###############################################################################

# cortical layer labels: e for excitatory; i for inhibitory
lname = ['L23e', 'L23i', 'L4e', 'L4i', 'L5e', 'L5i','L6e', 'L6i']

# number of neurons by layer
n_layer = [0, 20683, 5834, 21915, 5479, 4850, 1065, 14395, 2948];
l_bins = np.cumsum(n_layer) # cumulative number of neurons by layer
N = np.sum(n_layer)         # total number of neurons

# creating dataframes to store measures
ainess = np.zeros((len(bg_rate),len(g))) # AIness index

for row in range(len(g)):
    for col in range(len(bg_rate)):

        # loading data to a DataFrame structure
        filename = '../data/data_raster_g'+str(g[row])+'_bgrate'+str(bg_rate[col])+'.dat'
        data = pd.read_csv(filename, sep=" ", header=None, names=['i','t'])

        # grouping spiking times for each neuron
        keys,values = data.sort_values(['i','t']).values.T
        ukeys,index=np.unique(keys,True)
        arrays=np.split(values,index[1:])

        spk_neuron = pd.DataFrame({'i':range(0,N),'t':[[]]*N})
        spk_neuron.iloc[ukeys.astype(int),1] = arrays

        # creating a flag to identify cortical layers
        spk_neuron['layer'] = pd.cut(spk_neuron['i'], l_bins, labels=lname, right=False)
        data['layer'] = pd.cut(data['i'], l_bins, labels=lname, right=False)

        # sampling data
        n_sample = 1000 # number of neurons by layer for sampled measures
        spk_neuron = spk_neuron.groupby(['layer']).apply(lambda x: x.sample(n=n_sample))

        # measures DataFrame:
        measures_layer = pd.DataFrame(index=lname)

        # cleaning variables
        del keys, values, ukeys, index, arrays

        ########################################################################
        # Mean firing rates
        ########################################################################
        freq = []
        freq = [float(len(spk_neuron.t[i]))/tsim for i in range(len(spk_neuron))]
        spk_neuron['f'] = freq

        measures_layer['f'] = spk_neuron.groupby(['layer'])['f'].mean()

        ########################################################################
        # Interspike intervals + coefficient of variation
        ########################################################################
        # interspike intervals
        isi = []
        isi = [np.diff(spk_neuron.t[i]) for i in range(len(spk_neuron))]

        # coefficient of variation
        cv = []
        cv = [np.std(isi[i])/np.mean(isi[i]) if len(isi[i])>1 else np.nan \
                for i in range(len(spk_neuron))]
        spk_neuron['cv'] = cv

        measures_layer['cv'] = spk_neuron.groupby(['layer'])['cv'].mean()

        ########################################################################
        # Synchrony index
        ########################################################################
        sync = []
        bins = np.arange(0,tsim*1000.0+3.0,3)

        for i in range(len(lname)):
            index_sample = spk_neuron.i[spk_neuron.layer.isin([lname[i]])]
            count, division = np.histogram(data.t[data.i.isin(index_sample)],bins=bins)

            # removing first 100 ms of simulation
            sync.append(np.var(count[166:])/np.mean(count[166:]))

        measures_layer['sync'] = sync

        ########################################################################
        # AIness measure: f<30Hz & 0.7<=cv<1.2 & sync_index <= 8
        ########################################################################
        measures_layer['AI'] = (measures_layer.f<30)&(measures_layer.sync<=8)&\
                                (measures_layer.cv>=0.7)&(measures_layer.cv<1.2)

        # % of layers in the AIness range
        ainess[col][row] = 100*sum(measures_layer.AI)/8.0

        # cleaning variables
        del measures_layer, spk_neuron, data

# saving files
np.savetxt('../data/ainess_sample.dat', ainess)

# ----------------------------------------------------------------------------
# File description:
#
# Code to construct figure 6: A) firing rate vs background rate; B) Colorgraph
# with AIness%; C) firing rate vs relative inhibitory synaptic strength
# ----------------------------------------------------------------------------
# loading data
data = pd.read_csv('../data/ainess_sample.dat', header = None, sep=' ')

plotstyle = ['-.','-^','-*','-s']

# figure size for the graphs
fig = plt.figure(figsize=(8,8))
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 0.2])

# fig. 8B: contour plot for %AIness
data = data.sort_index(ascending=False)

plt.subplot(gs[0,0])
levels1 = np.linspace(0, 100, num=5, endpoint='True')
CS1 = plt.contour(g, bg_rate[::-1], data, levels=levels1, colors='k', linestyles='dashed')
plt.clabel(CS1, inline=1, fontsize=10, fmt='%1.0f')
CS2 = plt.imshow(data, extent=[min(g),max(g),min(bg_rate),max(bg_rate)], cmap=cm.jet, interpolation='gaussian', aspect='auto')

plt.xlabel('relative inh. synaptic strength')
plt.ylabel('background rate [Hz]')

# colorbar of contour plot
ax4 = plt.subplot(gs[0,1])
plt.colorbar(CS2, cax=ax4, ticks=np.arange(0,101,50))
plt.title('%AIness')

plt.tight_layout()
plt.savefig('./figures/fig6_sample.pdf', dpi=600)
