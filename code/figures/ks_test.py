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
# P-value calculated with Kolmogorov-Smirnov statistic.
# Figures with cummulated histograms of firing rate and CV's distribution.
# ----------------------------------------------------------------------------

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as sc

###############################################################################
# Filenames
###############################################################################
filename = (['raster_brian.dat','raster_nest.dat'])

###############################################################################
# Parameters
###############################################################################
tsim = 60.5     # time of simulation in seconds
n = 77169       # total number of neurons

# cortical layer labels: e for excitatory; i for inhibitory
lname = ['L23e', 'L23i', 'L4e', 'L4i', 'L5e', 'L5i','L6e', 'L6i']

# number of neurons by layer
n_layer = [0, 20683, 5834, 21915, 5479, 4850, 1065, 14395, 2948];
l_bins = np.cumsum(n_layer) # cumulative number of neurons by layer

# dataframe to store CV's and firing rates
measures = pd.DataFrame()

for k in range(0,2):

    # loading data
    if k == 0: data = pd.read_csv(filename[k], sep=" ", header=None, names=['i','t'])
    else:
        data = pd.read_csv(filename[k], sep="\t", header=None, names=['i','t','lixo'])
        data = data.drop('lixo',1)
        data.i = data.i - 1

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
    # n_sample = 1000 # number of neurons by layer for sampled measures
    # spk_neuron = spk_neuron.groupby(['layer']).apply(lambda x: x.sample(n=1000))

    ###############################################################################
    # Interspike intervals + coefficient of variation
    ###############################################################################
    # interspike intervals
    isi = []
    isi = [np.diff(spk_neuron.t[i]) for i in range(len(spk_neuron))]

    # coefficient of variation
    aux = []
    aux = [np.std(isi[i])/np.mean(isi[i]) if len(isi[i])>1 else np.nan\
            for i in range(len(spk_neuron))]

    cv = np.zeros(n)*np.nan
    cv[spk_neuron.i.astype(int)] = aux

    if k == 1:
        measures['cv_nest'] = cv
    else:
        measures['cv_brian'] = cv

    ###############################################################################
    # Firing rates
    ###############################################################################
    aux = []
    aux = [float(len(spk_neuron.t[i]))/tsim for i in range(len(spk_neuron))]

    freq = np.zeros(n)
    freq[spk_neuron.i.astype(int)] = aux

    if k == 1:
        measures['f_nest'] = freq
    else:
        measures['f_brian'] = freq

    # cleaning variables
    del data, spk_neuron

measures['layer'] = pd.cut(measures.index, l_bins, labels=lname, right=False)

# arrays to store the p-value from Kolmogorov-Smirnov test
ks_cv = np.zeros(8)
ks_f = np.zeros(8)

# create figures to plot cummulative histograms
fig1, ax1 = plt.subplots(nrows=2, ncols=4, figsize=(10,6))
fig2, ax2 = plt.subplots(nrows=2, ncols=4, figsize=(10,6))

for i in range(0,8):

    # grid positions for the graph
    xgrid = int(i%2)
    ygrid = int(i/2)

    ################################################################
    # CV's distribution comparison
    ################################################################
    plt.figure(1)
    plt.subplot2grid((2,4),(xgrid,ygrid))
    plt.title(lname[i])

    # excluding NAN values
    aux1 = measures.groupby('layer')['cv_brian'].get_group(lname[i])
    aux1 = aux1[~np.isnan(aux1)]

    aux2 = measures.groupby('layer')['cv_nest'].get_group(lname[i])
    aux2 = aux2[~np.isnan(aux2)]

    # CV's histograms
    cv1,edges1 = np.histogram(aux1,bins='doane',range=(0,2.0))
    cv2,edges2 = np.histogram(aux2,bins=np.linspace(0,2.0,len(edges1)))

    # cummulated histograms
    cv1 = np.cumsum(cv1)
    cv2 = np.cumsum(cv2)

    # normalizing cummulated histograms
    cv1 = cv1/float(max(cv1))
    cv2 = cv2/float(max(cv2))

    # plot data
    plt.plot(edges1[:-1],cv1,label='brian')
    plt.plot(edges2[:-1],cv2,label='nest')
    if (i==0): plt.legend(loc=4)
    if (i>1): plt.yticks([])
    if (i%2==0): plt.xticks([])
    plt.autoscale(enable=True, axis='x', tight=True)

    _,ks_cv[i] = sc.ks_2samp(cv1,cv2)
    del aux1, aux2, cv1, cv2

    ################################################################
    # firing rate distribution comparison
    ################################################################
    plt.figure(2)
    plt.subplot2grid((2,4),(xgrid,ygrid))
    plt.title(lname[i])

    # excluding NAN values
    aux1 = measures.groupby('layer')['f_brian'].get_group(lname[i])
    aux1 = aux1[~np.isnan(aux1)]

    aux2 = measures.groupby('layer')['f_nest'].get_group(lname[i])
    aux2 = aux2[~np.isnan(aux2)]

    # firing rate histograms
    f1,edges1 = np.histogram(aux1,bins='doane',range=(0,max(aux1)))
    f2,edges2 = np.histogram(aux2,bins=np.linspace(0,max(aux1),len(edges1)))

    # cummulated histograms
    f1 = np.cumsum(f1)
    f2 = np.cumsum(f2)

    # normalizing cummulated histograms
    f1 = f1/float(max(f1))
    f2 = f2/float(max(f2))

    # plot data
    plt.plot(edges1[:-1],f1,label='brian')
    plt.plot(edges1[:-1],f2,label='nest')
    if (i==0): plt.legend(loc=4)
    if (i>1): plt.yticks([])
    plt.autoscale(enable=True, axis='x', tight=True)

    _,ks_f[i] = sc.ks_2samp(f1,f2)
    del aux1, aux2, f1, f2

# formating figures:

plt.figure(1)
fig1.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.xlabel('CV', fontsize=12, labelpad=10)
plt.ylabel('cumulative distribution',fontsize=12, labelpad=12)
plt.tight_layout()
plt.savefig('cv_hists.jpg',dpi=600)

plt.figure(2)
fig2.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.xlabel('firing rate [Hz]', fontsize=12, labelpad=10)
plt.ylabel('cumulative distribution',fontsize=12, labelpad=12)
plt.tight_layout()
plt.savefig('f_hists.jpg',dpi=600)

# save Kolmogorov-Smirnov test results
np.savetxt('ks_cv.dat',ks_cv)
np.savetxt('ks_freq.dat',ks_f)
