from brian2 import *
from numba import jit
import numpy

@jit(nopython=True) # numba decorator compiles this function into low level code to run faster
def make_single_train(min_rate, max_rate, max_time_wo_spike, max_change_speed, runduration, number_neurons, dt_createpattern):
	runduration1 = min(runduration, 150)  # input will be tripled later, only need 150s
	st = []
	virtual_pre_sim_spike = - numpy.random.rand() * max_time_wo_spike
	firing_rate = min_rate + numpy.random.rand() * (max_rate - min_rate)
	rate_change = 2 * (numpy.random.rand() - .5) * max_change_speed
	mtws = max_time_wo_spike
	for t in numpy.arange(dt_createpattern, runduration1, dt_createpattern):
		if numpy.random.rand() < dt_createpattern * firing_rate or \
			(len(st) < 1 and t - virtual_pre_sim_spike > mtws) or \
			(len(st) > 0 and t - st[-1] > mtws):
			tmp = t - dt_createpattern * numpy.random.rand()
			tmp = max(0, tmp)
			tmp = min(runduration1, tmp)
			st.append(tmp)
			mtws = max_time_wo_spike
		firing_rate = firing_rate + rate_change * dt_createpattern
		rate_change = rate_change + 1 / 5 * 2 * (numpy.random.rand() - .5) * max_change_speed
		rate_change = max(min(rate_change, max_change_speed), -max_change_speed)
		firing_rate = max(min(firing_rate, max_rate), min_rate)
	return array(st)


def make_input(min_rate, max_rate, max_time_wo_spike, max_change_speed, runduration, number_neurons, dt_createpattern):
	# make input spikes
	spiketimes = []
	afferents = []
	runduration1 = min(runduration, 150)  # input will be tripled later, only need 150s
	for n in range(number_neurons):
		st = make_single_train(min_rate, max_rate, max_time_wo_spike, max_change_speed, runduration, number_neurons, dt_createpattern)
		spiketimes.append(st)
		afferents.append(n*ones(len(st)))
	spiketimes = hstack(spiketimes)
	afferents = hstack(afferents)
	sortarray = argsort(spiketimes)
	spiketimes = spiketimes[sortarray]
	afferents = afferents[sortarray]
	return afferents, spiketimes


def make_pattern_presentation_array(runduration, patternlength, pattern_freq):
	# make array with a 1 (pattern) or 0 (not pattern) for each time window (window size = length of pattern (0.05s))
	runduration1 = min(runduration, 150)
	if pattern_freq == 0.5:
		# if pattern frequency is 0.5, then this array s just alternating 0's and 1's
		position_copypaste = array([0, 1] * runduration1 * int(pattern_freq / patternlength))
	else:
		position_copypaste = zeros(int(runduration1 / patternlength))
		while sum(position_copypaste) < floor(int(runduration1 / patternlength) * pattern_freq):
			random_index = numpy.random.randint(0, len(position_copypaste))
			if position_copypaste[random_index] == 0:
				if random_index > 0 and random_index < int(runduration1 / patternlength) - 1 and position_copypaste[
					random_index - 1] == 0 and position_copypaste[random_index + 1] == 0:
					position_copypaste[random_index] = 1
				elif random_index == 0 and position_copypaste[random_index + 1] == 0:
					position_copypaste[random_index] = 1
				elif random_index == int(runduration1 / patternlength) - 1 and position_copypaste[random_index - 1] == 0:
					position_copypaste[random_index] = 1
	return position_copypaste


def copy_and_paste_jittered_pattern(times, indices, position_copypaste, patternlength, jitter_sd, spike_del, number_pat):
	# choose first segment as pattern to be copy-pasted
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
				tim1_add = numpy.random.rand(sum(invert(keep_array)))*patternlength
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


def remove_close_spikes(times, indices, dt):
	# remove spikes that are too close in time, depends on time resolution chosen for simulation
	last_spike = -2 * ones(len(set(indices)))
	keep_flag = ones(len(times), dtype=bool)
	# calculations of spike distance
	for j, st in enumerate(times):
		if st - last_spike[int(indices[j])] < dt:
			keep_flag[j] = False
		else:
			last_spike[int(indices[j])] = st
	print('Number of spikes to be deleted: ', len(indices) - sum(keep_flag))
	times = times[keep_flag]
	indices = indices[keep_flag]
	return times, indices
