# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Pamela Hathway p.hathway16@imperial.ac.uk
# ----------------------------------------------------------------------------- 
# Runs the main pattern finding algorithm and saved the results.
# -----------------------------------------------------------------------------

from brian2 import *
import time as time
import argparse
import pickle
from numba import jit
import numpy

set_device('cpp_standalone', build_on_run=True, directory='STDP_standalone')

argparser = argparse.ArgumentParser()
argparser.add_argument("runid")
args = argparser.parse_args()
runid = args.runid
rest, numbr = runid.split('=')
seedyy = int(numbr)
savedir = '../../data'  # this should be the cluster directory where the saved files should go
para = load('../../data/para.npy')  # where the para.npy file is located

seedy = mod(seedyy, 100)  # get random seed for simulation, keep it between 1 and 100 as for other simulations
pararow = seedyy // 100

win = para[pararow, 1]  # 1.9*T/2000
jitter_sd = para[pararow, 2]
number_pat = int(para[pararow, 3])
pattern_freq = para[pararow, 4]
spike_del = para[pararow, 5]
T = para[pararow, 6]  #  0.5*(1-spike_del)*number_pat
K2 = 3  # 4
A = - K2 * T  # 4 * T  # 1500
refract = 1

name_run = ''
defaultclock.dt = 1e-4 * second  # length of time step for simulation
runduration = 450  # length of simulation [second], full length = 450
tripling = True
now = datetime.datetime.now()

result_onerun = {}
result_onerun['name'] = name_run
result_onerun['runid'] = runid
result_onerun['seedyy'] = seedyy
result_onerun['win'] = win
result_onerun['jitter'] = jitter_sd
result_onerun['n_pat'] = number_pat
result_onerun['freq_pat'] = pattern_freq
result_onerun['del'] = spike_del
result_onerun['threshold'] = T
result_onerun['runduration'] = runduration
result_onerun['date'] = now.strftime('%y%m%d')
result_onerun['t_res'] = defaultclock.dt / second
result_onerun['K2'] = K2
result_onerun['A'] = A
result_onerun['refract_N2'] = refract

with open('%s/result_seed%s_%s.pickle' % (savedir, runid, name_run), 'wb') as handle:
    pickle.dump(result_onerun, handle, protocol=pickle.HIGHEST_PROTOCOL)


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


print('####################################')
print('Simulation parameters')
print('####################################')
print('Name =                   ', name_run)
print('dt =                     ', defaultclock.dt)
print('initial weight =         ', win)
print('jitter (SD) =            ', jitter_sd)
print('% of neurons in pattern =', 100 * number_pat / 2000)
print('pattern freq =           ', pattern_freq)
print('% spikes deleted =       ', spike_del * 100)
print('random seed =            ', seedy)
print('runduration =            ', runduration)
print('threshold =              ', T)
print('refractory =             ', refract)
print('K2 =                     ', K2)
print('A =                      ', A)
print('####################################')
print('Create input')
print('####################################')

# parameters for making input - exactly same as in Masquelier code
number_neurons = 2000
dt_createpattern = 0.001
number_npat = number_neurons - number_pat
patternlength = 0.05
# parameters for basic spikes (2000 neurons, firing rate changes over time, avg 60Hz)
max_rate_pat = 90
min_rate_pat = 0
max_time_wo_spike_pat = 0.05
max_change_speed_pat = max_rate_pat / max_time_wo_spike_pat
# parameters for noise spikes (2000 neurons, firing rate does not change, avg 10 Hz)
max_rate_add = 10
min_rate_add = 10
max_time_wo_spike_add = 1000
max_change_speed_add = 0

# LIF output neuron model parameters
taum = 10 * ms
taus = 2.5 * ms
tausyn = 2.5 * ms
u_rest = 0
X = (taus / taum) ** (taum / (taus - taum))
# A = -1500
deltax = 1
deltaa = 1
deltau = 1.0

# STDP synapse parameters
tauplus = 16.8 * ms
tauminus = 33.7 * ms
aplus = 2.0 ** -5
aminus = 0.85 * aplus
wmin = 0
wmax = 1

# equations
eqs = '''du/dt = (A*a)/taus + (X*x-u)/taum : 1
        dx/dt = -x/tausyn : 1
        da/dt = -a/taus : 1'''

eqs_reset = '''x = 0
                u = 2*T
                a = deltaa'''

stdp_model = ''' wi : 1
                dLTPtrace/dt = -LTPtrace / tauplus  : 1 (event-driven)
                dLTDtrace/dt = -LTDtrace / tauminus : 1 (event-driven)'''

stdp_on_pre = (''' x_post += deltax*wi
                    LTPtrace = aplus
                    wi = clip(wi + LTDtrace, wmin, wmax)
                    LTDtrace = 0''')

stdp_on_post = ('''LTDtrace = -aminus
					wi = clip(wi + LTPtrace, wmin, wmax)
					LTPtrace = 0''')

# stdp_on_post = ('''LTDtrace = -aminus
#                wi = clip(wi + LTPtrace + 0.5*(aminus-aplus)*(lastspike_pre==lastspike_post), wmin, wmax)
#                LTPtrace = 0''')

# Make input
numpy.random.seed(int(seedyy))
indices, times = make_input(min_rate_pat, max_rate_pat, max_time_wo_spike_pat,
                            max_change_speed_pat, runduration, number_neurons, dt_createpattern, seedyy)
indices_add, times_add = make_input(min_rate_add, max_rate_add, max_time_wo_spike_add,
                                    max_change_speed_add, runduration, number_neurons, dt_createpattern, seedyy)
position_copypaste = make_pattern_presentation_array(runduration, patternlength, pattern_freq, seedyy)
times, indices = copy_and_paste_jittered_pattern(times, indices, position_copypaste,
                                                 patternlength, jitter_sd, spike_del, number_pat, seedyy)
times, indices = add_noise(times, indices, times_add, indices_add)
if tripling and runduration > 300:
    times, indices = triple_input_runtime(times, indices)
    position_copypaste = concatenate((position_copypaste, position_copypaste, position_copypaste))
timing_pattern = where(position_copypaste > 0)[0] * patternlength
times, indices = remove_close_spikes(times, indices, defaultclock.dt / second)
times = times * second

# Make neuron layers N0(2000 input spikes) N1(2000 presynpatic neurons) N2(1 postsynaptic neuron)
N0 = SpikeGeneratorGroup(number_neurons, indices, times)
N1 = NeuronGroup(number_neurons, '''du/dt  = -u/taum : 1''', threshold='u > T', reset='u = u_rest', refractory=0 * ms,
                 method='linear')
N1spm = SpikeMonitor(N1)
N2 = NeuronGroup(1, eqs, threshold='u > T', reset=eqs_reset, refractory=refract * ms, method='linear')
N2spm = SpikeMonitor(N2)
# N2stm = StateMonitor(N2, ['u'], record=True)

# Make synapses
syn01 = Synapses(N0, N1, on_pre='''u_post += T + 100''', method='linear')
syn01.connect(j='i')
syn12 = Synapses(N1, N2, model=stdp_model, on_pre=stdp_on_pre, on_post=stdp_on_post, method='linear')
syn12.connect(i=range(0, number_neurons), j=0)
syn12.wi = win
# syn12stm = StateMonitor(syn12, ['wi'], record=range(0, 2000), dt=100*ms)

print('####################################')
print('Simulation')
print('####################################')

net = Network(collect())
net.add(N0, N1, N2, N2spm, syn01, syn12)
net.run(runduration * second, report='text')

print('####################################')
print('Results')
print('####################################')

print('N2 spikes                 = ', len(N2spm.t))
result_onerun['N2_spikes'] = len(N2spm.t)

# # Calculate latency
eval_last_sec = min(runduration / 3, 150)
latency = empty(0)
# if pattern is not present at start, add zeros (for each N2 spike) at beginning of latency array
latency = concatenate((latency, zeros(len(where(N2spm.t / second < timing_pattern[0])[0]))))
for i2 in range(len(timing_pattern)):
    if i2 < len(timing_pattern) - 1:
        in_window = [p for indx, p in enumerate(N2spm.t / second) if
                     p >= timing_pattern[i2] and p < timing_pattern[i2 + 1]]
    else:
        in_window = [p for indx, p in enumerate(N2spm.t / second) if p >= timing_pattern[i2]]
    latency = concatenate((latency, in_window - timing_pattern[i2]))
latency[where(latency > patternlength)] = 0

# calculate average latency of last 150s
lat_end = latency[where(N2spm.t / second > runduration - eval_last_sec)[0]]
if len(lat_end) > sum(position_copypaste) * 0.8 / 3:
    avg_lat = 1000 * mean(lat_end[lat_end > 0])
    result_onerun['avg_lat'] = avg_lat
else:
    print('only %s spikes in last 3rd of simulation' % len(lat_end))
    avg_lat = -1
    result_onerun['avg_lat'] = -1
print('Avg latency               = ', round(result_onerun['avg_lat'], 3))

fa = sum(lat_end == 0)
spikes50msbins, binns = histogram(N2spm.t / second, bins=arange(0, runduration + patternlength, patternlength))
true_hits_end = spikes50msbins[int(eval_last_sec / patternlength):][
    position_copypaste[int(eval_last_sec / patternlength):] > 0]
hits = sum(true_hits_end > 0) / len(true_hits_end)
result_onerun['hits'] = round(hits, 3)
result_onerun['fa'] = fa
if hits > 0.98 and fa == 0 and avg_lat < 10 and avg_lat > 0:
    result_onerun['success'] = 1
else:
    result_onerun['success'] = 0

print('Hit rate (>98)            = ', result_onerun['hits'])
print('Number false alarms (!=0) = ', result_onerun['fa'])
print('Success                   = ', result_onerun['success'])

# see when pattern was found
if hits > 0.9:
    lat_begin = latency[where(N2spm.t / second < 50)[0]]
    find_spike = where(lat_begin == 0)[0][-1] + 1
    find_t = N2spm.t[find_spike] / second
    result_onerun['find_t'] = round(find_t, 2)
    result_onerun['find_spike'] = find_spike
else:
    result_onerun['find_t'] = -1
    result_onerun['find_spike'] = -1
print('find_t                    = ', result_onerun['find_t'])
print('find_spike                = ', result_onerun['find_spike'])

with open('%s/result_seed%s_%s.pickle' % (savedir, runid, name_run), 'wb') as handle:
    pickle.dump(result_onerun, handle, protocol=pickle.HIGHEST_PROTOCOL)

print('####################################')
print('End of simulation')
print('####################################')
