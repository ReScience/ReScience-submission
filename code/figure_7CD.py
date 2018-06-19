# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Pamela Hathway p.hathway16@imperial.ac.uk
# ----------------------------------------------------------------------------- 
# Creates figure 7 CD
# This script can be run independently.
# -----------------------------------------------------------------------------

from brian2 import *
import time as time


def fig_7CD(run_idx):
    from create_input import make_input, make_pattern_presentation_array, \
        copy_and_paste_jittered_pattern, add_noise, triple_input_runtime, remove_close_spikes

    device.reinit()
    device.activate()
    set_device('cpp_standalone', build_on_run=True, directory='STDP_standalone_ata3')
    print('%s Preparing Figure 7 C and D: success with ATA rule' % time.strftime('%H:%M'))

    # parameters that were varied for the publication
    random_seed = run_idx
    defaultclock.dt = 1e-4 * second  # length of time step for simulation

    jitter_sd = 0  # amount of jitter
    number_pat = 1000  # number of neurons that take part in the pattern presentation
    number_neurons = 2000
    pattern_freq = 0.25  # frequency at which the pattern is presented
    spike_del = 0  # percentage of spikes within the pattern to be deleted
    T = 500  # threshold
    win = 0.2375  # initial weight for 2000 input neurons - this number was not changed

    # parameters for making input - same as in Masquelier et al., 2008
    runduration = 150  # length of simulation [s], full length = 450
    tripling = False  # instead of creating 450 s of unique input, make 150 s and concatenate those spikes
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

    # LIF output neuron model parameters
    taum = 10 * ms
    taus = 2.5 * ms
    tausyn = 2.5 * ms
    u_rest = 0
    X = (taus / taum) ** (taum / (taus - taum))
    deltax = 1
    deltaa = 1
    deltau = 1.754
    K2 = 3
    A = - K2 * T
    refract = 1

    # STDP synapse parameters
    tauplus = 16 * ms
    tauminus = 33 * ms
    aplus = 0.006
    aminus = 0.5 * aplus
    wmin = 0
    wmax = 0.5

    # equations
    eqs = '''du/dt = (-(u-u_rest))/taum : 1'''

    eqs_reset = '''u = -50'''

    stdp_model = ''' wi : 1
					dLTPtrace/dt = -LTPtrace / tauplus  : 1 (event-driven)
					dLTDtrace/dt = -LTDtrace / tauminus : 1 (event-driven)'''

    stdp_on_pre = (''' u_post += deltau*wi
						LTPtrace += aplus
						wi = clip(wi + LTDtrace, wmin, wmax)''')

    stdp_on_post = ('''LTDtrace -= aminus
						wi = clip(wi + LTPtrace, wmin, wmax)''')

    numpy.random.seed(int(random_seed))
    indices, times = make_input(min_rate_pat, max_rate_pat, max_time_wo_spike_pat,
                                max_change_speed_pat, runduration, number_neurons, dt_createpattern, random_seed)
    indices_add, times_add = make_input(min_rate_add, max_rate_add, max_time_wo_spike_add,
                                        max_change_speed_add, runduration, number_neurons, dt_createpattern,
                                        random_seed)
    position_copypaste = make_pattern_presentation_array(runduration, patternlength, pattern_freq, random_seed)
    times, indices = copy_and_paste_jittered_pattern(times, indices, position_copypaste,
                                                     patternlength, jitter_sd, spike_del, number_pat, random_seed)
    times, indices = add_noise(times, indices, times_add, indices_add)
    if tripling and runduration > 300:
        times, indices = triple_input_runtime(times, indices)
        position_copypaste = concatenate((position_copypaste, position_copypaste, position_copypaste))
    timing_pattern = where(position_copypaste > 0)[0] * patternlength
    times, indices = remove_close_spikes(times, indices, defaultclock.dt / second)
    times = times * second

    # Make neuron layers N0(input spikes) N1(2000 input neurons) N2(1 output neuron)
    # Carry out the sort by hand because it's more efficient than the brian version
    I = lexsort((indices, times))
    indices = indices[I]
    times = times[I]
    N0 = SpikeGeneratorGroup(number_neurons, indices, times, sorted=True)
    N1 = NeuronGroup(number_neurons, '''du/dt  = -u/taum : 1''', threshold='u > T', reset='u = u_rest',
                     refractory=0 * ms, method='linear')
    N2 = NeuronGroup(1, eqs, threshold='u > T', reset=eqs_reset, refractory=refract * ms, method='linear')
    N2spm = SpikeMonitor(N2)

    # Make synapses
    syn01 = Synapses(N0, N1, on_pre='''u_post += 600''', method='linear')
    syn01.connect(j='i')
    syn12 = Synapses(N1, N2, model=stdp_model, on_pre=stdp_on_pre, on_post=stdp_on_post, method='linear')
    syn12.connect(i=range(0, number_neurons), j=0)
    syn12.wi = win

    # Some speed hacks
    N0._previous_dt = N0.dt_[:]
    N0._spikes_changed = False

    net = Network(collect())
    net.add(N0, N1, N2, N2spm, syn01, syn12)
    net.run(runduration * second)

    # # Calculate latency for RNN
    latency = []
    # if pattern is not present at start, add zeros (for each N2 spike) at beginning of latency array
    latency.append(zeros(len(where(N2spm.t / second < timing_pattern[0])[0])))
    for i2 in range(len(timing_pattern)):
        if i2 < len(timing_pattern) - 1:
            in_window = [p for indx, p in enumerate(N2spm.t / second) if
                         p >= timing_pattern[i2] and p < timing_pattern[i2 + 1]]
        else:
            in_window = [p for indx, p in enumerate(N2spm.t / second) if p >= timing_pattern[i2]]
        latency.append(in_window - timing_pattern[i2])
    latency = hstack(latency)
    latency[where(latency > patternlength)] = 0

    # calculate what spikes fall into what bin
    spikes50msbins, binns = histogram(N2spm.t / second, bins=arange(0, 150.05, 0.05))
    # calc if those spikes were within pattern or not
    true = position_copypaste * spikes50msbins
    # get sum of spikes per second
    correctsecond = empty(runduration)
    for i in range(runduration):
        numspikesinsec = sum(spikes50msbins[i * 20:i * 20 + 20])
        correctsecond[i] = sum(true[i * 20:i * 20 + 20]) / numspikesinsec

    fig7 = figure(figsize=(12, 5))
    suptitle('Figure 7 C and D')
    subplot(121)
    plot(N2spm.t / second, latency * 1000, 'b.', markersize=3)
    ylabel('Postsynaptic spike latency [ms]')
    xlabel('time [s]')

    subplot(122)
    plot(range(runduration), 100 * correctsecond, 'g-', markersize=3)
    ylabel('% correct spikes')
    xlabel('time [s]')
    xlim([0, runduration])
    ylim([-5, 105])
    yticks([0, 50, 100])
    xticks([0, 50, 100, 150])

    # tight_layout()

    savefig('../article/figures/figure_7CD_seed%s_created_%s.pdf' % (
    random_seed, datetime.datetime.now().strftime('%y%m%d')))

    return fig7


if __name__ == '__main__':
    # example seeds that show a similar behaviour to the one shown in fig 7CD in the article:
    # 13, 20, 28, 49
    run_idx = 28
    fig7CD = fig_7CD(run_idx)
    show()
