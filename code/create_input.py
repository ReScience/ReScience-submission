# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Pamela Hathway p.hathway16@imperial.ac.uk
# -----------------------------------------------------------------------------
# These functions are called by run_simulation.py to create the spike trains.
# They were adapted to Python from the Matlab code from
# Masquelier T, Guyonneau R, Thorpe SJ (2008). Spike Timing Dependent Plasticity Finds the Start of
# Repeating Patterns in Continuous Spike Trains. PLoS ONE 3(1): e1377 doi:10.1371/journal.pone.0001377

# This script can be run independently.
# -----------------------------------------------------------------------------

from brian2 import *
from numba import jit
import numpy
from numpy import array


def make_input(min_rate, max_rate, max_time_wo_spike, max_change_speed, runduration, number_neurons, dt_createpattern, random_seed):
    # this code was adapted from the original Matlab code
    # this function creates input spikes one neuron at a time
    spiketimes = []
    afferents = []
    runduration1 = min(runduration, 150)  # input will be tripled later, only need 150s
    for n in range(number_neurons):
        # call function to make a single spike train with a random seed for numba function
        # (seed for this number is defined in run_simulation.py, so this is deterministic)
        st = numpy.array(
            make_single_train(min_rate, max_rate, max_time_wo_spike, max_change_speed, runduration, number_neurons,
                              dt_createpattern, numpy.random.randint(2**30)))
        spiketimes.append(st)
        afferents.append(n * ones(len(st)))
    # convert spike trains into array and sort by time
    spiketimes = hstack(spiketimes)
    afferents = hstack(afferents)
    sortarray = argsort(spiketimes)
    spiketimes = spiketimes[sortarray]
    afferents = afferents[sortarray]
    return afferents, spiketimes


@jit(nopython=True)  # numba decorator compiles this function into low level code to run faster
def make_single_train(min_rate, max_rate, max_time_wo_spike, max_change_speed, runduration, number_neurons,
                      dt_createpattern, random_seed):
    # this code was adapted from the original Matlab code
    # get random seed from make_input() (changes for every neuron)
    numpy.random.seed(int(random_seed))
    runduration1 = min(runduration, 150)  # input will be tripled later, only need 150s
    st = []
    # add a virtual spike at beginning to avoid all neurons spiking at 0
    virtual_pre_sim_spike = - numpy.random.rand() * max_time_wo_spike
    # firing rates change at each time step
    firing_rate = min_rate + numpy.random.rand() * (max_rate - min_rate)
    rate_change = 2 * (numpy.random.rand() - .5) * max_change_speed
    mtws = max_time_wo_spike
    # create spike trains from changing firing rates: go through 150s 1ms at a time and
    # add spikes according to current firing rate (changes every 1ms)
    for t in numpy.arange(dt_createpattern, runduration1, dt_createpattern):
        if numpy.random.rand() < dt_createpattern * firing_rate or \
                (len(st) < 1 and t - virtual_pre_sim_spike > mtws) or \
                (len(st) > 0 and t - st[-1] > mtws):
            # add a spike
            tmp = t - dt_createpattern * numpy.random.rand()
            tmp = max(0, tmp)
            tmp = min(runduration1, tmp)
            st.append(tmp)
            mtws = max_time_wo_spike
        # calculate firing rate change
        firing_rate = firing_rate + rate_change * dt_createpattern
        rate_change = rate_change + 1 / 5 * 2 * (numpy.random.rand() - .5) * max_change_speed
        rate_change = max(min(rate_change, max_change_speed), -max_change_speed)
        firing_rate = max(min(firing_rate, max_rate), min_rate)
    return st


def make_pattern_presentation_array(runduration, patternlength, pattern_freq, random_seed):
    # make array with a 1 (pattern) or 0 (not pattern) for each time window (window size = length of pattern (0.05s))
    runduration1 = min(runduration, 150)
    if pattern_freq == 0.5:
        # if pattern frequency is 0.5, then this array is just alternating 0's and 1's
        position_copypaste = array([0, 1] * runduration1 * int(pattern_freq / patternlength))
    else:
        # add a number of 1's (depends on frequency of pattern presentation)
        # into a zeros array, but avoid consecutive 1's
        position_copypaste = zeros(int(runduration1 / patternlength))
        while sum(position_copypaste) < floor(int(runduration1 / patternlength) * pattern_freq):
            random_index = numpy.random.randint(0, len(position_copypaste))
            if position_copypaste[random_index] == 0:
                if random_index > 0 and random_index < int(runduration1 / patternlength) - 1 and position_copypaste[
                    random_index - 1] == 0 and position_copypaste[random_index + 1] == 0:
                    position_copypaste[random_index] = 1
                elif random_index == 0 and position_copypaste[random_index + 1] == 0:
                    position_copypaste[random_index] = 1
                elif random_index == int(runduration1 / patternlength) - 1 and position_copypaste[
                    random_index - 1] == 0:
                    position_copypaste[random_index] = 1
    return position_copypaste


def copy_and_paste_jittered_pattern(times, indices, position_copypaste, patternlength, jitter_sd, spike_del,
                                    number_pat, random_seed):
    # choose first segment as pattern to be the pattern spikes that are copy-pasted at all other times
    startCPindex = where(position_copypaste == 1)[0][0]
    # get times and indices in pattern window
    tim = times[
          searchsorted(times, startCPindex * patternlength):searchsorted(times, (startCPindex + 1) * patternlength)]
    ind = indices[
          searchsorted(times, startCPindex * patternlength):searchsorted(times, (startCPindex + 1) * patternlength)]
    tim = tim[ind < number_pat]
    ind = ind[ind < number_pat]
    tim -= startCPindex * patternlength
    # replace indices and times for 50ms when position_copypaste == 1 with ind, tim
    indices1 = []
    times1 = []
    counti2 = 0
    for in2, i2 in enumerate(position_copypaste):
        ind1 = copy(ind)
        tim1 = copy(tim)
        if i2 == 1:  # pattern is present
            counti2 += 1
            if spike_del > 0:
                # randomly delete percentage of spike in pattern
                keep_array = numpy.random.rand(len(ind))
                keep_array = keep_array > spike_del
                ind1 = ind[keep_array]
                tim1 = tim[keep_array]
                # add random spikes to the same neurons to keep the spike density constant
                ind1_add = ind[invert(keep_array)]
                tim1_add = numpy.random.rand(sum(invert(keep_array))) * patternlength
                ind1 = concatenate((ind1, ind1_add))
                tim1 = concatenate((tim1, tim1_add))
            indices1.append(ind1)
            if jitter_sd > 0:
                jitter = numpy.random.normal(0, jitter_sd, len(tim1))
            else:
                jitter = zeros(len(tim1))
            tim_jit = tim1 + jitter / 1000
            # if we are at the start of simulation, avoid negative spike times
            if in2 == 0:
                tim_jit[tim_jit < 0] = 0

            times1.append(tim_jit + in2 * patternlength)
            start_pattern1 = searchsorted(times, in2 * patternlength)
            end_pattern1 = searchsorted(times, (in2 + 1) * patternlength)
            tim_npat = times[start_pattern1:end_pattern1]
            ind_npat = indices[start_pattern1:end_pattern1]
            indices1.append(ind_npat[ind_npat >= number_pat])
            times1.append(tim_npat[ind_npat >= number_pat])
        else:
            # find index where pattern window starts
            start_pattern1 = searchsorted(times, in2 * patternlength)
            # find index where pattern window ends
            end_pattern1 = searchsorted(times, (in2 + 1) * patternlength)
            indices1.append(indices[start_pattern1:end_pattern1])
            times1.append(times[start_pattern1:end_pattern1])
    indices1 = hstack(indices1)
    times1 = hstack(times1)
    # sort input according to time
    sortarray = times1.argsort()
    indices1 = indices1[sortarray]
    times1 = times1[sortarray]
    return times1, indices1


def add_noise(times, indices, times_add, indices_add):
    # combine the basic activity and the 10Hz additional noise to one input
    times = concatenate((times, times_add))
    indices = concatenate((indices, indices_add))
    sortarray = argsort(times)
    times = times[sortarray]
    indices = indices[sortarray]
    return times, indices


def triple_input_runtime(times, indices):
    # To shorten time spent on creating input, 150s input is tripled to give 450s
    times = concatenate((times, times + 150, times + 300))
    indices = concatenate((indices, indices, indices))
    return times, indices


@jit(nopython=True)
def remove_close_spikes(times, indices, dt):
    # remove spikes that are too close in time, depends on time resolution chosen for simulation
    last_spike = -2 * numpy.ones(int(numpy.amax(indices) + 1))
    keep_flag = numpy.ones(len(times), dtype=numpy.bool_)
    # calculations of spike distance
    for j, st in enumerate(times):
        if st - last_spike[int(indices[j])] < dt:
            keep_flag[j] = False
        else:
            last_spike[int(indices[j])] = st
    # print('    Number of spikes to be deleted: ', len(indices) - sum(keep_flag), 'or', round(100*(len(indices) - sum(keep_flag))/len(indices), 2), '%')
    times = times[keep_flag]
    indices = indices[keep_flag]
    return times, indices


if __name__ == '__main__':
    random_seed = 2
    defaultclock.dt = 1e-4 * second  # length of time step for simulation

    jitter_sd = 1  # amount of jitter
    number_pat = 1000  # number of neurons that take part in the pattern presentation
    pattern_freq = 0.25  # frequency at which the pattern is presented
    spike_del = 0  # percentage of spikes within the pattern to be deleted
    T = 0.5 * (1 - spike_del) * number_pat  # threshold
    win = 1.9 * T / 2000  # initial weight for 2000 input neurons - this number was not changed

    # parameters for making input - same as in Masquelier et al., 2008
    runduration = 450  # length of simulation [s], full length = 450
    tripling = True  # instead of creating 450 s of unique input, make 150 s and concatenate those spikes
    number_neurons = 2000
    dt_createpattern = 0.001  # [s]
    number_npat = number_neurons - number_pat  # number of neurons not in pattern
    patternlength = 0.05  # [s]
    # parameters for basic spikes (2000 neurons, firing rate changes over time 0-90 Hz, avg 54Hz)
    max_rate_pat = 90
    min_rate_pat = 0
    max_time_wo_spike_pat = 0.05
    max_change_speed_pat = max_rate_pat / max_time_wo_spike_pat
    # parameters for additional noise spikes (2000 neurons, firing rate does not change, avg 10 Hz)
    max_rate_add = 10
    min_rate_add = 10
    max_time_wo_spike_add = 1000
    max_change_speed_add = 0

    print('Start making input with seed', random_seed, time.strftime('%H:%M'))
    # seed(int(random_seed))
    indices, times = make_input(min_rate_pat, max_rate_pat, max_time_wo_spike_pat, max_change_speed_pat, runduration,
                                number_neurons, dt_createpattern, random_seed)
    print('Start making noise with seed', random_seed, time.strftime('%H:%M'))
    indices_add, times_add = make_input(min_rate_add, max_rate_add, max_time_wo_spike_add, max_change_speed_add,
                                        runduration, number_neurons, dt_createpattern, random_seed)
    position_copypaste = make_pattern_presentation_array(runduration, patternlength, pattern_freq, random_seed)
    print('Start copy paste pattern', time.strftime('%H:%M'))
    times, indices = copy_and_paste_jittered_pattern(times, indices, position_copypaste, patternlength, jitter_sd,
                                                     spike_del, number_pat, random_seed)
    print('Start adding noise', time.strftime('%H:%M'))
    times, indices = add_noise(times, indices, times_add, indices_add)
    if tripling and runduration > 300:
        print('Starting tripling input', time.strftime('%H:%M'))
        times, indices = triple_input_runtime(times, indices)
        position_copypaste = concatenate((position_copypaste, position_copypaste, position_copypaste))
    timing_pattern = where(position_copypaste > 0)[0] * patternlength
    print('Start removing too-close spikes', time.strftime('%H:%M'))
    times, indices = remove_close_spikes(times, indices, defaultclock.dt / second)
    times = times * second
    patneurons = range(0, number_pat)
    nonpatneurons = range(number_pat, number_neurons)

    print(times[0:30])
    print(indices[0:30])
