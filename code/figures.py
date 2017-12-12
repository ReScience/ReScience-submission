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
# Code to construct figures 2, 5 and 6 from the article
# ----------------------------------------------------------------------------

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec

import pandas as pd
import numpy as np

# ----------------------------------------------------------------------------
# Function description:
#
# Code to construct figure 2: A) Raster plot; B) box plot of firing rates;
# C) barplot of CV's; D) barplot of synchrony index.
# ----------------------------------------------------------------------------
def createfig2(tsim, filename):

    ###############################################################################
    # Loading data and defining parameters
    ###############################################################################
    data = pd.read_csv(filename, sep=" ",
                        header=None, names=['i','t'])

    # cortical layer labels: e for excitatory; i for inhibitory
    lname = ['L23e', 'L23i', 'L4e', 'L4i', 'L5e', 'L5i','L6e', 'L6i']

    # number of neurons by layer
    n_layer = [0, 20683, 5834, 21915, 5479, 4850, 1065, 14395, 2948];
    l_bins = np.cumsum(n_layer) # cumulative number of neurons by layer
    N = np.sum(n_layer)         # total number of neurons

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

    spk_neuron = pd.DataFrame({'i':range(0,N),'t':[[]]*N})
    spk_neuron.iloc[ukeys.astype(int),1] = arrays

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
    plt.savefig('./figures/fig2.pdf', dpi=600)

# ----------------------------------------------------------------------------
# Function description:
#
# Code to construct figure 5: A1 and B1 Raster plots of network activity;
# A2 and B2 barplots of firing rates
# ----------------------------------------------------------------------------
def createfig5(tsim, filename):

    ###############################################################################
    # Filenames
    ###############################################################################
    #             DC current                      Layer-independent
    filename = (filename)

    ###############################################################################
    # Parameters
    ###############################################################################
    # cortical layer labels: e for excitatory; i for inhibitory
    lname = ['L23e', 'L23i', 'L4e', 'L4i', 'L5e', 'L5i','L6e', 'L6i']

    # number of neurons by layer
    n_layer = [0, 20683, 5834, 21915, 5479, 4850, 1065, 14395, 2948];
    l_bins = np.cumsum(n_layer) # cumulative number of neurons by layer
    N = np.sum(n_layer)         # total number of neurons

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

        spk_neuron = pd.DataFrame({'i':range(0,N),'t':[[]]*N})
        spk_neuron.iloc[ukeys.astype(int),1] = arrays

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
        plt.xlabel('firing rate [Hz]')
        plt.gca().invert_yaxis()

        # cleaning variables
        del data, measures_layer, spk_neuron

    plt.tight_layout()
    plt.subplots_adjust(left=0.05)
    plt.savefig('./figures/fig5.pdf', dpi=600)

def createfig5_hist(tsim):

    s = np.arange(0,100);            # seed of pseudo-random number

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
    freqs = pd.DataFrame(index=lname)      # mean firing rates by layer

    # measures DataFrame:
    measures_layer = pd.DataFrame(index=lname)

    for row in range(len(s)):

        # loading data to a DataFrame structure
        filename = '../data/data_raster_g4.0_bgrate8.0_bg_random'+str(s[row])+'.dat'
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

        # cleaning variables
        del keys, values, ukeys, index, arrays

        ########################################################################
        # Mean firing rates
        ########################################################################
        freq = []
        freq = [float(len(spk_neuron.t[i]))/tsim for i in range(len(spk_neuron))]
        spk_neuron['f'] = freq

        measures_layer['f_'+str(s[row])] = spk_neuron.groupby(['layer'])['f'].mean()

        # cleaning variables
        del spk_neuron, data

    freqs = measures_layer.T

    fig, ax = plt.subplots(nrows=4, ncols=1, figsize=(4,10))
    hist_lim = np.array([10, 10, 30, 14])

    for i in range(0,len(lname),2):
        plt.subplot(4,1,1+i/2)
        freqs[lname[i]].plot.hist(bins=np.arange(0,hist_lim[i/2]), histtype='stepfilled', alpha=0.6)
        freqs[lname[i+1]].plot.hist(bins=np.arange(0,hist_lim[i/2]), histtype='stepfilled', alpha=0.6)
        plt.ylabel('#trials')
        plt.xlim(0,hist_lim[i/2])
        plt.legend()

    plt.xlabel('population firing rate [Hz]')
    plt.tight_layout()
    plt.savefig('./figures/fig5_hist.pdf', dpi=600)

# ----------------------------------------------------------------------------
# File description:
#
# Calculation of the %AIness index used to construct figure 6
# ----------------------------------------------------------------------------
def datafig6(tsim):

    g = np.arange(2.0,10.5,0.5);            # relative inh. synaptic strength values
    bg_rate = np.arange(3.0,15.5,0.5);  # background rate values

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
    freq_g4 = pd.DataFrame(index=lname)      # firing rates for g = 4.0
    freq_bg8 = pd.DataFrame(index=lname)     # firing rates for bg = 8.0 Hz

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
    np.savetxt('../data/ainess.dat', ainess)
    freq_g4.to_csv('../data/freq_g4.csv')
    freq_bg8.to_csv('../data/freq_bg8.csv')

# ----------------------------------------------------------------------------
# File description:
#
# Code to construct figure 6: A) firing rate vs background rate; B) Colorgraph
# with AIness%; C) firing rate vs relative inhibitory synaptic strength
# ----------------------------------------------------------------------------
def createfig6():

    # loading data
    data = pd.read_csv('../data/ainess.dat', header = None, sep=' ')
    freq_g4 = pd.read_csv('../data/freq_g4.csv', index_col=0); freq_g4 = freq_g4.T
    freq_bg8 = pd.read_csv('../data/freq_bg8.csv', index_col=0); freq_bg8 = freq_bg8.T

    bg = np.arange(3.0,15.5,0.5)
    g = np.arange(2,10.5,0.5)

    plotstyle = ['-.','-^','-*','-s']

    # figure size for the graphs
    fig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(2, 3, width_ratios=[1, 3, 0.2], height_ratios=[3, 1])

    # fig. 8A: firing rates for fixed parameter g = 4.0
    ax1 = plt.subplot(gs[0,0])
    plt.plot(freq_g4.L23e,bg, plotstyle[0])
    plt.plot(freq_g4.L4e,bg, plotstyle[1])
    plt.plot(freq_g4.L5e,bg, plotstyle[2])
    plt.plot(freq_g4.L6e,bg, plotstyle[3])
    plt.ylim([3,15])
    plt.gca().invert_xaxis()
    plt.xlabel('firing rate [Hz]')
    plt.xlim(16, 0)

    # fig. 8C: firing rates for fixed parameter background rate = 8.0Hz
    ax2 = plt.subplot(gs[1,1])
    freq_bg8.plot(y=['L23e','L4e','L5e','L6e'], ax=ax2, style=plotstyle)
    plt.gca().invert_yaxis()
    plt.ylabel('firing rate [Hz]')

    # fig. 8B: contour plot for %AIness
    data = data.sort_index(ascending=False)
    bg = np.arange(15.5,3.0,-0.5)

    ax3 = plt.subplot(gs[0,1])
    levels1 = np.linspace(0, 100, num=5, endpoint='True')
    CS1 = plt.contour(g, bg, data, levels=levels1, colors='k', linestyles='dashed')
    plt.clabel(CS1, inline=1, fontsize=10, fmt='%1.0f')
    CS2 = plt.imshow(data, extent=[2,10,3,15], cmap=cm.jet, interpolation='gaussian', aspect='auto')

    plt.xlabel('relative inh. synaptic strength')
    plt.ylabel('background rate [Hz]')

    # colorbar of contour plot
    ax4 = plt.subplot(gs[0,2])
    plt.colorbar(CS2, cax=ax4, ticks=np.arange(0,101,50))
    plt.title('%AIness')

    plt.tight_layout()
    plt.savefig('./figures/fig6.pdf', dpi=600)
