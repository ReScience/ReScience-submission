# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Pamela Hathway p.hathway16@imperial.ac.uk
# ----------------------------------------------------------------------------- 
# Runs the main pattern finding algorithm and creates figures 3, 4, 7AB.
# This script can be run independently.
# -----------------------------------------------------------------------------

from brian2 import *
import time as time


def run_sim(run_idx, pararow=0, only_success=False):
	from figure_3_4_7AB import fig_3, fig_4, fig_7AB
	from create_input import make_input, make_pattern_presentation_array, \
		copy_and_paste_jittered_pattern, add_noise, triple_input_runtime, remove_close_spikes

	device.reinit()
	device.activate()
	set_device('cpp_standalone', build_on_run=True, directory='STDP_standalone')

	if only_success == True:
		run_idx = pararow * 100 + run_idx

	# parameters that were varied for the publication
	random_seed = run_idx
	defaultclock.dt = 1e-4 * second  # length of time step for simulation

	if only_success == False:
		jitter_sd = 1  # amount of jitter
		number_pat = 1000  # number of neurons that take part in the pattern presentation
		number_neurons = 2000
		pattern_freq = 0.25  # frequency at which the pattern is presented
		spike_del = 0  # percentage of spikes within the pattern to be deleted
		T = 0.5*(1-spike_del)*number_pat  # threshold
		win = 1.9*T/number_neurons  # initial weight for 2000 input neurons - this number was not changed
	else:
		para = load('../data/para.npy')
		win = para[pararow, 1]  # 1.9*T/2000
		jitter_sd = para[pararow, 2]
		number_pat = int(para[pararow, 3])
		pattern_freq = para[pararow, 4]
		spike_del = para[pararow, 5]
		T = para[pararow, 6]  # 0.5*(1-spike_del)*number_pat
		number_neurons = 2000


	if only_success == False:
		print('%s Preparing Figure 3 and 4 and 7AB: latency and convergence and weights' % time.strftime('%H:%M'))
		print('    #### Simulation parameters:')
		print('    random seed =            ', random_seed)
		print('    dt =                     ', defaultclock.dt)
		print('    initial weight =         ', win)
		print('    jitter (SD) =            ', jitter_sd)
		print('    % of neurons in pattern =', 100*number_pat/number_neurons)
		print('    pattern freq =           ', pattern_freq)
		print('    % spikes deleted =       ', spike_del*100)


	# parameters for making input - same as in Masquelier et al., 2008
	runduration = 450  # length of simulation [s], full length = 450
	tripling = True  # instead of creating 450 s of unique input, make 150 s and concatenate those spikes
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
	X = (taus/taum)**(taum/(taus - taum))
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

	stdp_on_pre_nn = (''' x_post += deltax*wi
						LTPtrace = aplus
						wi = clip(wi + LTDtrace, wmin, wmax)''')

	stdp_on_post_nn = ('''LTDtrace = -aminus
						wi = clip(wi + LTPtrace, wmin, wmax)''')

	stdp_on_pre_ata = (''' x_post += deltax*wi
						LTPtrace += aplus
						wi = clip(wi + LTDtrace, wmin, wmax)''')

	stdp_on_post_ata = ('''LTDtrace -= aminus
						wi = clip(wi + LTPtrace, wmin, wmax)''')

	if only_success == False:
		print('    #### Creating input')
	start_time_input = time.time()
	seed(int(random_seed))
	indices, times = make_input(min_rate_pat, max_rate_pat, max_time_wo_spike_pat,
								max_change_speed_pat, runduration, number_neurons, dt_createpattern)
	indices_add, times_add = make_input(min_rate_add, max_rate_add, max_time_wo_spike_add,
										max_change_speed_add, runduration, number_neurons, dt_createpattern)
	position_copypaste = make_pattern_presentation_array(runduration, patternlength, pattern_freq)
	times, indices = copy_and_paste_jittered_pattern(times, indices, position_copypaste,
													 patternlength, jitter_sd, spike_del, number_pat)
	times, indices = add_noise(times, indices, times_add, indices_add)
	if tripling and runduration > 300:
		times, indices = triple_input_runtime(times, indices)
		position_copypaste = concatenate((position_copypaste, position_copypaste, position_copypaste))
	timing_pattern = where(position_copypaste > 0)[0] * patternlength
	times, indices = remove_close_spikes(times, indices, defaultclock.dt / second)
	times = times * second

	start_time_simulation = time.time()
	# Make neuron layers N0(input spikes) N1(2000 input neurons) N2(1 output neuron)
	# Carry out the sort by hand because it's more efficient than the brian version
	I = lexsort((indices, times))
	indices = indices[I]
	times = times[I]
	N0 = SpikeGeneratorGroup(number_neurons, indices, times, sorted=True)
	N1 = NeuronGroup(number_neurons, '''du/dt  = -u/taum : 1''', threshold='u > T', reset='u = u_rest', refractory=0*ms, method='linear')
	N1spm = SpikeMonitor(N1)

	# Output neuron using RNN STDP rule
	N2 = NeuronGroup(1, eqs, threshold='u > T', reset=eqs_reset, refractory=refract * ms, method='linear')
	N2spm = SpikeMonitor(N2)
	# N2stm = StateMonitor(N2, ['u'], record=True)
	if only_success == False:
		# Output neuron using NN STDP rule
		N2nn = NeuronGroup(1, eqs, threshold='u > T', reset=eqs_reset, refractory=refract * ms, method='linear')
		N2nnspm = SpikeMonitor(N2nn)
		# Output neuron using ATA STDP rule
		N2ata = NeuronGroup(1, eqs, threshold='u > T', reset=eqs_reset, refractory=refract * ms, method='linear')
		N2ataspm = SpikeMonitor(N2ata)


	# Make synapses
	syn01 = Synapses(N0, N1, on_pre='''u_post += T + 100''', method='linear')
	syn01.connect(j='i')
	# RNN
	syn12 = Synapses(N1, N2, model=stdp_model, on_pre=stdp_on_pre, on_post=stdp_on_post, method='linear')
	syn12.connect(i=range(0, number_neurons), j=0)
	syn12.wi = win
	syn12stm = StateMonitor(syn12, ['wi'], record=range(0, 2000), dt=1*second)
	if only_success == False:
		# NN
		syn12nn = Synapses(N1, N2nn, model=stdp_model, on_pre=stdp_on_pre_nn, on_post=stdp_on_post_nn, method='linear')
		syn12nn.connect(i=range(0, number_neurons), j=0)
		syn12nn.wi = win
		syn12nnstm = StateMonitor(syn12nn, ['wi'], record=range(0, 2000), dt=1*second)
		# ATA
		syn12ata = Synapses(N1, N2ata, model=stdp_model, on_pre=stdp_on_pre_ata, on_post=stdp_on_post_ata, method='linear')
		syn12ata.connect(i=range(0, number_neurons), j=0)
		syn12ata.wi = win
		syn12atastm = StateMonitor(syn12ata, ['wi'], record=range(0, 2000), dt=1*second)

	# Some speed hacks
	N0._previous_dt = N0.dt_[:]
	N0._spikes_changed = False


	if only_success == False:
		print('    #### Simulation (ca. 70s)')
		net = Network(collect())
		net.add(N0, N1, syn01, N2, N2spm, syn12, syn12stm,
				N2nn, N2nnspm, syn12nn, syn12nnstm,
				N2ata, N2ataspm, syn12ata, syn12atastm)
	else:
		net = Network(collect())
		net.add(N0, N1, syn01, N2, N2spm, syn12, syn12stm)

	net.run(runduration * second)


	start_time_results = time.time()
	if only_success == False:
		print('    #### Results')

	# # Calculate latency for RNN
	eval_last_sec = min(runduration/3, 150)
	latency = []
	# if pattern is not present at start, add zeros (for each N2 spike) at beginning of latency array
	latency.append(zeros(len(where(N2spm.t / second < timing_pattern[0])[0])))
	for i2 in range(len(timing_pattern)):
		if i2 < len(timing_pattern)-1:
			in_window = [p for indx, p in enumerate(N2spm.t / second) if p >= timing_pattern[i2] and p < timing_pattern[i2+1]]
		else:
			in_window = [p for indx, p in enumerate(N2spm.t / second) if p >= timing_pattern[i2]]
		latency.append(in_window - timing_pattern[i2])
	latency = hstack(latency)
	latency[where(latency > patternlength)] = 0


	# ### Analyse results
	# calculate average latency of last 150s and print results
	lat_end = latency[where(N2spm.t / second > runduration - eval_last_sec)[0]]
	if len(lat_end) > 500:
		avg_lat = 1000*mean(lat_end[lat_end > 0])
		if only_success == False:
			print('    Avg latency               = ', round(avg_lat, 3))
	else:
		if only_success == False:
			print('    only %s spikes in last 3rd of simulation - this run was probably unsuccessful' % len(lat_end))
		avg_lat = -1

	# calculate whether run was successful or not
	fa = sum(lat_end == 0)
	spikes50msbins, binns = histogram(N2spm.t / second, bins=arange(0, runduration + patternlength, patternlength))
	true_hits_end = spikes50msbins[int(eval_last_sec/patternlength):][position_copypaste[int(eval_last_sec/patternlength):] > 0]
	hits = sum(true_hits_end > 0) / len(true_hits_end)
	if hits > 0.98 and fa == 0 and avg_lat < 10 and avg_lat > 0:
		success = 1
	else:
		success = 0

	if only_success == False:
		print('    Hit rate (>98)            = ', hits)
		print('    Number false alarms (!=0) = ', fa)
		print('    Success                   = ', success)

	# ###see when pattern was found
	if hits > 0.9 and fa < 50:
		lat_begin = latency[where(N2spm.t / second < 50)[0]]
		find_spike = where(lat_begin == 0)[0][-1] + 1
		find_t = N2spm.t[find_spike] / second
		if only_success == False:
			print('    find_t                    = ', round(find_t, 2))
			print('    find_spike                = ', find_spike)
	else:
		find_t = -1
		find_spike = -1
		if only_success == False:
			print('    it does not look like the pattern was found successfully')


	if only_success == False:
		print('    #### Make figures')
		fig3 = fig_3(N2spm, latency)  # plot the latency of all output neuron spikes

		fig4 = fig_4(N1spm, syn12, timing_pattern)  # plot the weights from the second to last pattern

		fig7AB = fig_7AB(syn12stm, syn12nnstm, syn12atastm)

		return fig3, fig4, fig7AB

	else:
		return success


if __name__=='__main__':
	run_sim(0)
	show()
