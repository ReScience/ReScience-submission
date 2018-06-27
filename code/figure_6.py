# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Pamela Hathway p.hathway16@imperial.ac.uk
# ----------------------------------------------------------------------------- 
# Creates figure 6 either from scratch (default) or using saved data (to switch,
# comment out lines 266-267 and comment in line 265). 
# This script can be run independently.
# -----------------------------------------------------------------------------

from brian2 import *
import pickle as pickle


def fig_6_saved():
    print('%s Preparing Figure 6: Time of finding pattern' % time.strftime('%H:%M'))
    plt.rcParams['savefig.dpi'] = 75
    plt.rcParams['figure.autolayout'] = False
    plt.rcParams['figure.figsize'] = 6, 5
    plt.rcParams['axes.labelsize'] = 17
    plt.rcParams['axes.titlesize'] = 17
    plt.rcParams['font.size'] = 17
    plt.rcParams['lines.linewidth'] = 2.0
    plt.rcParams['lines.markersize'] = 5
    plt.rcParams['legend.fontsize'] = 16

    resolutions = [10**(-4), 10**(-5), 10**(-6)]

    # ### load arrays from saved runs
    with open('../data/results_saved_0_e4.pickle', 'rb') as handle:
            run = pickle.load(handle)
    # only load runs that were successful in finding
    find_time_e4 = run['find_t'][run['success'] > 0]

    with open('../data/results_saved_0_e5.pickle', 'rb') as handle:
        run = pickle.load(handle)
    # only load runs that were successful in finding
    find_time_e5 = run['find_t'][run['success'] > 0]

    with open('../data/results_saved_0_e6.pickle', 'rb') as handle:
            run = pickle.load(handle)
    # only load runs that were successful in finding
    find_time_e6 = run['find_t'][run['success'] > 0]

    # ### plot time finding
    fig = figure(figsize=(6, 5))
    title('Figure 6 saved results')
    errorbar(resolutions, [mean(find_time_e4), mean(find_time_e5), mean(find_time_e6)], [std(find_time_e4), std(find_time_e5), std(find_time_e6)], color='k')
    plot(resolutions, [mean(find_time_e4), mean(find_time_e5), mean(find_time_e6)], 'ko')
    ylabel('pattern finding time [s]')
    xlabel('simulation time step [s]')
    semilogx(10**(-4), 10**(-6), 'k-')
    xlim([0.0000006, 0.00015])
    ylim([0, 30])
    # yticks(array([12, 14, 16, 18, 20, 22]))
    axhline(y=13.5, linestyle='--', color='g')
    tight_layout()

    savefig('../article/figures/figure_6_from_saved_created_%s.pdf' % datetime.datetime.now().strftime('%y%m%d'))

    return fig


def fig_6_new(reps_per_resolution):
    print('Preparing Figure 6: Time of finding pattern')
    plt.rcParams['savefig.dpi'] = 75
    plt.rcParams['figure.autolayout'] = False
    plt.rcParams['figure.figsize'] = 6, 5
    plt.rcParams['axes.labelsize'] = 17
    plt.rcParams['axes.titlesize'] = 17
    plt.rcParams['font.size'] = 17
    plt.rcParams['lines.linewidth'] = 2.0
    plt.rcParams['lines.markersize'] = 3
    plt.rcParams['legend.fontsize'] = 16

    # call run_find 5 times with different resolutions and set runtime to 30seconds
    # plot results
    resolutions = [1e-4, 1e-5, 1e-6]
    find_t_all = zeros((len(resolutions), reps_per_resolution))
    find_spike_all = zeros((len(resolutions), reps_per_resolution))

    for res_idx, res in enumerate(resolutions):
        # res = resolutions[0]
        # rep = 0
        for rep in range(0, reps_per_resolution):
            find_t, find_spike = run_find(res, rep)
            find_t_all[res_idx, rep] = find_t
            find_spike_all[res_idx, rep] = find_spike

    mean_find_t = sum(find_t_all, axis=1)/sum(find_t_all > 0, axis=1)
    std_find_t = array([std(find_t_all[0, find_t_all[0, :] > 0]), std(find_t_all[1, find_t_all[1, :] > 0]), std(find_t_all[2, find_t_all[2, :] > 0])])

    fig = figure(figsize=(6, 5))
    title('Figure 6')
    errorbar(resolutions, mean_find_t, std_find_t, color='k')
    plot(resolutions, mean_find_t, 'ko')
    ylabel('pattern finding time [s]')
    xlabel('simulation time step [s]')
    semilogx(10**(-4), 10**(-6), 'k-')
    xlim([0.0000006, 0.00015])
    ylim([0, 30])
    # yticks(array([12, 14, 16, 18, 20, 22]))
    axhline(y=13.5, linestyle='--', color='g')
    tight_layout()

    savefig('../article/figures/figure_6_created_%s.pdf' % datetime.datetime.now().strftime('%y%m%d'))

    return fig


def run_find(res, rep):
    print('    Start rep %s at resolution %s' % (rep, res))
    from create_input import make_input, make_pattern_presentation_array, \
        copy_and_paste_jittered_pattern, add_noise, triple_input_runtime, remove_close_spikes

    device.reinit()
    device.activate()
    set_device('cpp_standalone', build_on_run=True, directory='STDP_standalone0')

    # parameters that were varied for the publication
    random_seed = rep
    defaultclock.dt = res * second  # length of time step for simulation

    jitter_sd = 1  # amount of jitter
    number_pat = 1000  # number of neurons that take part in the pattern presentation
    number_neurons = 2000
    pattern_freq = 0.25  # frequency at which the pattern is presented
    spike_del = 0  # percentage of spikes within the pattern to be deleted
    T = 0.5 * (1 - spike_del) * number_pat  # threshold
    win = 1.9 * T / number_neurons  # initial weight for 2000 input neurons - this number was not changed

    # parameters for making input - same as in Masquelier et al., 2008
    runduration = 30  # length of simulation [s], full length = 450
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
    deltau = 1.0
    K2 = 3
    A = - K2 * T
    refract = 1

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

    numpy.random.seed(int(random_seed))
    indices, times = make_input(min_rate_pat, max_rate_pat, max_time_wo_spike_pat,
                                max_change_speed_pat, runduration, number_neurons, dt_createpattern, random_seed)
    indices_add, times_add = make_input(min_rate_add, max_rate_add, max_time_wo_spike_add,
                                        max_change_speed_add, runduration, number_neurons, dt_createpattern, random_seed)
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

    # Output neuron using RNN STDP rule
    N2 = NeuronGroup(1, eqs, threshold='u > T', reset=eqs_reset, refractory=refract * ms, method='linear')
    N2spm = SpikeMonitor(N2)

    # Make synapses
    syn01 = Synapses(N0, N1, on_pre='''u_post += T + 100''', method='linear')
    syn01.connect(j='i')
    # RNN
    syn12 = Synapses(N1, N2, model=stdp_model, on_pre=stdp_on_pre, on_post=stdp_on_post, method='linear')
    syn12.connect(i=range(0, number_neurons), j=0)
    syn12.wi = win
    syn12stm = StateMonitor(syn12, ['wi'], record=range(0, 2000), dt=1 * second)

    # Some speed hacks
    N0._previous_dt = N0.dt_[:]
    N0._spikes_changed = False

    net = Network(collect())
    net.add(N0, N1, syn01, N2, N2spm, syn12, syn12stm)
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

    # calculate whether run was successful or not
    spikes50msbins, binns = histogram(N2spm.t / second, bins=arange(0, runduration + patternlength, patternlength))

    # see when pattern was found
    if sum(latency[-10:] > 0) == 10:
        find_spike = where(latency == 0)[0][-1] + 1
        find_t = N2spm.t[find_spike] / second
    else:
        find_t = -1
        find_spike = -1

    return find_t, find_spike


if __name__=='__main__':
    # fig_6_saved()
    reps_per_resolution = 10
    fig_6_new(reps_per_resolution)

    show()

