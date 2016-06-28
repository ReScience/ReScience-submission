# -----------------------------------------------------------------------------
#
# Distributed under the GNU General Public License.
#
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
# Contains functions used for analysis of simulated data
# -----------------------------------------------------------------------------

import numpy as np


def f_Sd_sort(Spike_data, num_neuron):
    '''
    Description
    -----------
    This function sorts population spike-times into spike trains of 
    individual neurons

    Parameters
    ----------
    Spike_data: dict{'times':array[|spike time|], 'senders': array[|neuron gid|]} 
                spiking data mixed from multiple neurons. 

    -num_neuron: float - full number of neurons recorded. 

    Returns
    -------
    Spike_sorted: dict('ids': list of neuronal ids, 
                       'times': list of array[|spike times|]}
                  Spiking data per individual neuron. Includes empty fields for 
                  non-spiking neurons. 
    '''

    Spike_sorted = {'ids': [], 'times': []}

    # set of IDs of neurons that emit at least 1 spike
    Id_list = list(set(Spike_data['senders']))

    for neuron_id in range(num_neuron):

        # for neurons that emitted at least 1 spike
        if neuron_id < len(Id_list):

            # mask to select elements related to particular neuron ID
            mask = Spike_data['senders'] == Id_list[neuron_id]

            Spike_sorted['ids'] += [Id_list[neuron_id]]
            Spike_sorted['times'] += [Spike_data['times'][mask]]

        # other neurons are considered silent
        else:
            Spike_sorted['ids'] += [None]
            Spike_sorted['times'] += [np.array([])]

    return Spike_sorted


def f_psth(Data, bin_width=1., t_min=0, t_max=1):
    '''
    Description
    -----------
    calculate population time-histogram: for every time bin [ms] calculate 
    average spiking rate for each spike train in Data.   

    Parameters
    ----------
    Data: list(1D array) - [ms] array[|spike times|] - list of spike trains
    t_min and t_max: float  - [ms] boundaries of time window where PSTH is built 
    bin_width : float    - [ms] bin width

    Return
    ------
    PSTH: dict -{'times':|bin times [ms]|,'rates':|neuron_id x psth rates [Hz]|}
    '''

    assert t_max > t_min + \
        bin_width, "(t_max - t_min) should include at least 1 bin_width"

    num_bins = int((t_max - t_min) / bin_width)
    t_max = t_min + num_bins * bin_width

    num_neurons = len(Data)
    PSTH = {'times': np.arange(t_min + bin_width / 2., t_max + bin_width / 2., bin_width),
            'rates': np.zeros([num_neurons, num_bins])}

    train_id = 0

    for Spike_train in Data:
        # element-wise approach: correspond each time value to time bin

        # select data inside time window
        mask = (Spike_train > t_min) * (Spike_train < t_max)
        Spike_train = Spike_train[mask]

        # bin data
        binned_train = map(int, np.floor(Spike_train - t_min) / bin_width)
        PSTH['rates'][train_id] = np.histogram(binned_train,
                                               bins=np.arange(num_bins + 1))[0]
        train_id += 1

    PSTH['rates'] *= 1. / bin_width * 1000.
    return PSTH
