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
# Calculation of the %AIness index used to construct figure 6
# ----------------------------------------------------------------------------

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import numpy as np

g = np.arange(2.0,11.0);            # relative inh. synaptic strength values
bg_rate = np.arange(3.0,15.5,0.5);  # background rate values

###############################################################################
# Parameters
###############################################################################
tsim = 10.0     # time of simulation in seconds

# cortical layer labels: e for excitatory; i for inhibitory
lname = ['L23e', 'L23i', 'L4e', 'L4i', 'L5e', 'L5i','L6e', 'L6i']

# number of neurons by layer
n_layer = [0, 20683, 5834, 21915, 5479, 4850, 1065, 14395, 2948];
l_bins = np.cumsum(n_layer) # cumulative number of neurons by layer

# creating dataframes to store measures
ainess = np.zeros((len(bg_rate),len(g))) # AIness index
freq_g4 = pd.DataFrame(index=lname)      # firing rates for g = 4.0
freq_bg8 = pd.DataFrame(index=lname)     # firing rates for bg = 8.0 Hz

for row in range(len(g)):
    for col in range(len(bg_rate)):

        # loading data to a DataFrame structure
        filename = 'data_raster_g'+str(g[row])+'_bgfreq'+str(bg_rate[col])+'.dat'
        data = pd.read_csv(filename, sep=" ", header=None, names=['i','t'])

        # grouping spiking times for each neuron
        keys,values = data.sort_values(['i','t']).values.T
        ukeys,index=np.unique(keys,True)
        arrays=np.split(values,index[1:])
        spk_neuron = pd.DataFrame({'i':ukeys,'t':[list(a) for a in arrays]})

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

        if g[row]==4.0: freq_g4[bg_rate[col]] = measures_layer.f
        if bg_rate[col] == 8.0: freq_bg8[g[row]] = measures_layer.f

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
            sync.append(np.var(count[33:])/np.mean(count[33:]))

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
np.savetxt('ainess.dat', ainess)
freq_g4.to_csv('freq_g4.csv')
freq_bg8.to_csv('freq_bg8.csv')
