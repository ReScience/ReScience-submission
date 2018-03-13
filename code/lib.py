from __future__ import print_function, division

from collections import Counter

import numpy as np
from scipy.special import comb
from joblib.memory import Memory
from joblib import Parallel, delayed
import brian2 as br2
from brian2 import mV, ms

mem = Memory(cachedir='./__joblib_cache__', verbose=0, compress=True)


def group_size_evolutions(group_size, n_rep=1, linear=True, dt=0.1*ms,
                          bias_correction=False, seeds=None):
    """Return the evolution of the group size after one cycle. Returns a
    tuple (g0, g1), containing the size of the group size at stimulation time
    (g0), and the group size in the next iteration (g1). If ``bias_correction``
    is set to True, both of these values are corrected for the bias introduced
    by the finite time step size."""
    group_ev = []
    for rep in range(n_rep):
        if seeds is None:
            seed = None
        else:
            seed = seeds[rep]
        if bias_correction:
            # Run without stimulation, but otherwise identical (same seed)
            results_no_stim = run_with_cache(56 * ms, group_size=0,
                                             stim_time=None, linear=linear,
                                             repetition=rep, dt=dt, seed=seed)
        results = run_with_cache(56 * ms, group_size=group_size,
                                 stim_time=50 * ms, linear=linear,
                                 repetition=rep, dt=dt, seed=seed)
        next_iteration = np.sum(np.abs(results['t'] - 55 * ms) < dt / 2)
        if bias_correction:
            # Use the actual number of synchronous spike at the time of the
            # stimulation -- due to the background activity, this is usually
            # higher than the size of the stimulation group
            g0 = np.sum(np.abs(results['t'] - 50 * ms) < dt / 2)
            # Only count the synchronous spikes that were added due to the
            # stimulation
            next_iteration_no_stim = np.sum(np.abs(results_no_stim['t'] - 55 * ms) < dt / 2)
            g1 = next_iteration - next_iteration_no_stim
        else:
            g0 = group_size
            g1 = next_iteration
        group_ev.append((g0, g1))
    return group_ev


def F(bins, v_hist, epsilon):
    """Probability that an input of strength epsilon elicits a spike."""
    theta = bins[-1]
    return np.sum(v_hist[bins[:-1] > theta-epsilon/mV])


def epsilon(n_ex, n_in, w_ex, w_in, linear=True):
    """Total input for a neuron with given number of connections and given
       connection weights."""
    in_input = -n_in*w_in
    ex_input = n_ex*w_ex
    if not linear:
        ex_input = np.clip(ex_input, 0*mV, 2*mV) + np.clip(2*(ex_input-2*mV),
                                                           0*mV, 4*mV)
    return in_input + ex_input


def P_conn(group_size, n_ex, n_in, p_0, p_ex, p_in):
    """Probability that a neuron receives exactly n_ex excitatory and n_in
       inhibitory connections."""
    n_remain = group_size - n_ex - n_in
    return (comb(group_size, n_ex+n_in, exact=True) * comb(n_ex+n_in, n_ex, exact=True) *
            (p_0*p_ex)**n_ex * (p_0*p_in)**n_in * (1 - p_0)**n_remain)


def calc_group_evolution(bins, v_hist, group_size, w_ex=0.2*mV, w_in=0.2*mV,
                         N=1000, p_0=0.3, p_in=0.5, p_ex=0.5, linear=True):
    """Calculate the expected group size for the next iteration, based on a
       distribution of voltages and the group size of the previous iteration.
    """
    P_spike = 0.0
    for n_ex in range(group_size):
        for n_in in range(group_size - n_ex):
            P_triggered = F(bins, v_hist, epsilon(n_ex, n_in, w_ex,
                                                  w_in, linear=linear))
            P_spike += P_triggered * P_conn(group_size, n_ex, n_in, p_0, p_ex, p_in)
    return (N - group_size)*P_spike

# Note that this is the same as adding the @mem.cache decorator to the function,
# but the decorator leaves the name of the function unchanged which leads to
# problems with joblib's Parallel/delayed mechanism
# (see https://github.com/joblib/joblib/issues/226)
calc_group_evolution_cached = mem.cache(calc_group_evolution)

def run(runtime=250*ms, stim_time=150*ms, N=1000,
        linear=True, w_ex=0.2*mV, w_in=0.2*mV, repetition=0,
        store_spikes=True, store_v=False, dt=0.1*ms, seed=None, group_size=100):
    """Run a simulation and return the resulting spiketrain.

    Note: the `repetition` argument is only used to distinguish multiple runs
    of the identical network -- otherwise joblib's cache would return an
    identical result for each of them."""
    # Change Brian's schedule to apply synaptic events before checking the
    # threshold, in line with the original publication
    br2.magic_network.schedule = ['start', 'groups', 'synapses', 'thresholds',
                                  'resets', 'end']
    # Set a seed different for every simulation
    br2.seed(seed)

    # Basic constants used in the model
    br2.defaultclock.dt = dt
    Theta = 16 * mV
    tau_m = 8 * ms
    gamma = 1 / tau_m
    V0 = 17.6 * mV / tau_m

    if linear:
        # Build the group of neurons to use
        eq = "dV/dt = -gamma*V + V0 : volt"
        G = br2.NeuronGroup(N, model=eq, threshold="V>Theta",
                            reset="V=0*mV", method='euler')
    else:
        eq = """dV/dt = -gamma*V + V0 : volt
                ve : volt
             """
        G = br2.NeuronGroup(N, model=eq, threshold="V>Theta",
                            reset="V=0*mV", method='euler')

    if stim_time is not None:
        # Build stimulation
        stim = br2.SpikeGeneratorGroup(1, [0], [stim_time - dt])
        stim_syn = br2.Synapses(stim, G, on_pre="V = 2*Theta")
        stim_syn.connect(i=0, j=np.arange(group_size))

    # Number of spikes initially "in transit" (between 1 and 50 (inclusive))
    n_transit_spikes = np.random.random_integers(50)
    # "In transit" means that the spikes arrive in the first 5ms, i.e. as if
    # triggered by events before the start of the simulation (the synaptic delay
    # is 5ms and therefore all spikes triggered by neurons within the network
    # will arrive later than that).
    transit_spikes = br2.SpikeGeneratorGroup(n_transit_spikes,
                                             np.arange(n_transit_spikes),
                                             np.random.rand(n_transit_spikes)*5*ms)

    # We'll assume the initial spikes are 50% excitatory and 50% inhibitory
    n_inh_transit_spikes = n_transit_spikes//2
    n_exc_transit_spikes = n_transit_spikes - n_inh_transit_spikes
    if n_inh_transit_spikes > 0:
        inh_transit_syn = br2.Synapses(G, G, 'w : volt (constant, shared)',
                                       on_pre='V -= w', delay=5 * ms - dt)
        inh_transit_syn.connect(i=np.arange(n_inh_transit_spikes),
                                j=np.random.randint(N, size=n_inh_transit_spikes))
        inh_transit_syn.w = w_in

    if linear:  # See explanations about linear vs. non-linear coupling below
        exc_transit_syn = br2.Synapses(G, G, 'w : volt (constant, shared)',
                                       on_pre='V += w', delay=5 * ms - dt)
    else:
        exc_transit_syn = br2.Synapses(G, G, 'w : volt (constant, shared)',
                               on_pre='ve += w', delay=5 * ms - dt)
    exc_transit_syn.connect(i=np.arange(n_exc_transit_spikes),
                            j=np.random.randint(N, size=n_exc_transit_spikes))
    exc_transit_syn.w = w_ex

    # Build recurrent synaptic connections
    connections = np.random.rand(N, N) < 0.3
    exc_or_inh = np.random.rand(N, N) < 0.5
    exc_i, exc_j = (connections & exc_or_inh).nonzero()
    inh_i, inh_j = (connections & ~exc_or_inh).nonzero()

    if linear:
        # For linearly summing delta synapses, we can directly update the
        # neuron's membrane potential
        exc_syn = br2.Synapses(G, G, 'w : volt (constant, shared)',
                               on_pre='V += w', delay=5 * ms - dt)
        exc_syn.connect(i=exc_i, j=exc_j)
        exc_syn.w = w_ex
    else:
        # For nonlinearly summing synapses, we do not target the membrane
        # potential directly, but store the sum of all synaptic activity in
        # a separate variable.
        exc_syn = br2.Synapses(G, G, 'w : volt (constant, shared)',
                               on_pre='ve += w', delay=5 * ms - dt)
        exc_syn.connect(i=exc_i, j=exc_j)
        exc_syn.w = w_ex
        # We have to calculate the change in v across all synapses first before
        # updating v, and then setting the impact of the excitatory synapses
        # back to 0 before the next time step.
        G.run_regularly('''
                        V += clip(ve, 0*mV, 2*mV) + clip(2*(ve-2*mV), 0*mV, 4*mV)
                        ve = 0*mV
                        ''', when='after_synapses')

    # Inhibitory synapses always sum linearly
    inh_syn = br2.Synapses(G, G, 'w : volt (constant, shared)',
                           on_pre='V -= w', delay=5 * ms - dt)
    inh_syn.connect(i=inh_i, j=inh_j)
    inh_syn.w = w_in

    # Set random initial conditiong
    G.V = np.random.rand(N) * Theta

    if store_spikes:
        spikes = br2.SpikeMonitor(G)
    else:
        rates = br2.PopulationRateMonitor(G)

    if store_v:
        v_mon = br2.StateMonitor(G, 'V', record=True, when='before_synapses')

    br2.run(runtime)

    results = {}
    if store_spikes:
        results['i'] = spikes.i[:]
        results['t'] = spikes.t[:]
    else:
        # Revert the normalization performed by PopulationRateMonitor, i.e.
        # return the absolute number of spikes
        results['rate'] = np.array(rates.rate[:]*dt*N, dtype=int)
    if store_v:
        bins = np.linspace(-Theta/8/mV, Theta/mV, 100, endpoint=True)
        results['v'], results['bins'] = np.histogram(v_mon.V/mV, bins=bins)
    return results

# (See note for calc_group_evolution_cached above)
run_with_cache = mem.cache(run)


def run_and_bin(ext, inh, linear, repetition, stim_time=150*ms, dt=0.1*ms,
                seed=None, group_size=100):
    """Run the network and bin the resulting spike rates with the `bin_spikes`
    function. Return the extracted information needed to classify runs into the
    four categories (see `classify_run`). This information is returned as a
    tuple ``(max_f_before, min_g_size, max_f_after)``, where:
    ``max_f_before``: Maximum firing rate before stimulation
    ``min_g_size``: Minimal "group size" (number of spikes in the pulse
                    triggered by the stimulation), for 10 pulses after the
                    stimulation.
    ``max_f_after``: Maximum firing rate after stimulation (excluding the
                     pulses triggered by the stimulation)"""
    # Run simulation
    results = run_with_cache(w_ex=ext, w_in=inh, linear=linear,
                             stim_time=stim_time, dt=dt, repetition=repetition,
                             store_spikes=False, seed=seed, group_size=group_size)
    rates = results['rate']
    delay_steps = int(round(5*ms/dt))
    stim_step = int(round(stim_time/dt))
    check_until = int(round(stim_time/dt)) + 11*delay_steps
    group_size = rates[stim_step + delay_steps:check_until:delay_steps]
    before_stim = rates[:stim_step]
    after_stim = rates[check_until:]
    # Exclude the bins with synchronous activity
    after_stim = after_stim[np.arange(len(after_stim)) % delay_steps != 0]

    result = (np.max(before_stim),
              np.min(group_size),
              np.max(after_stim))
    return result

# (See note for calc_group_evolution_cached above)
run_and_bin_cached = mem.cache(run_and_bin)


@mem.cache
def do_repetitions(ext, inh, linear, n_rep, dt=0.1*ms, seeds=None,
                   group_size=100):
    """Perform ``n_rep`` repeated simulation runs for the given parameters and
    return a list of the results (see `run_and_bin` for the result format).
    The simulation are run in parallel with multiprocessing."""
    if seeds is None:
        seeds = [None]*n_rep
    # Note that "n_jobs=-2" means: use all available processor except 1
    return Parallel(n_jobs=-2)(delayed(run_and_bin_cached)(ext, inh, linear, i,
                                                           dt=dt, seed=seeds[i],
                                                           group_size=group_size)
                               for i in range(n_rep))


@mem.cache
def grid_search(weights_exc=np.arange(0.16, 0.4, 0.375/150.)*mV,
                weights_inh=None,
                n_rep=1, linear=True, quiet=False, dt=0.1*ms, seeds=None,
                group_size=100):
    """Perform a grid search over the given parameter range for excitatory and
    inhibitory synaptic connection strengths. If only the parameter range for
    excitatory connections is specified, then the same range is used for
    inhibitory connections.
    Returns a dictionary mapping each combination of excitatory and inhibitory
    weights (provided as a tuple ``(exc, inh)`` to a list of results (as
    returned by `run_and_bin`)."""
    if weights_inh is None:
        weights_inh = weights_exc
    results = {}
    for exc in weights_exc:
        if not quiet:
            print('exc=%.3fmV (%d values for inh, '
                  '%d repetitions each): ' % (exc/mV, len(weights_inh), n_rep),
                  end='')
        for inh in weights_inh:
            print('.', end='')
            results[(exc/mV, inh/mV)] = do_repetitions(exc, inh, linear, n_rep,
                                                       dt=dt, seeds=seeds,
                                                       group_size=group_size)
        print()
    return results


def classify_run(max_firing_before, min_group_size, max_firing_after):
    """
    Classify the results of a single run (see `run_and_bin` for the meaning
    of the three arguments) as either:
    ``'UNSTABLE_BEFORE'``: network is already unstable before stimulation
    ``'UNSTABLE_AFTER'``: network is unstable after stimulation
    ``'PERSISTENT'``: stimulation is persistently propagated
    ``'SHORT-LIVED'``: network is stable, but stimulation is not propagated
    In this context, "unstable" means that more than 10% of the population fire
    within the same time step, and "persistently propagated" means that the
    stimulus-triggered pulses are distinguishable from the maximum background
    activity for at least 10 pulses after stimulation.
    """
    if max_firing_before >= 100:
        # unstable "spontaneously"
        return 'UNSTABLE_BEFORE'
    if max_firing_after >= 100:
        # unstable after initiation of synchronous activity
        return 'UNSTABLE_AFTER'
    if min_group_size > max_firing_after:
        # group distinguishable from background --> persistent propagation
        return 'PERSISTENT'
    # Stable background activity, no persistent propagation
    return 'SHORT-LIVED'


def decide_color(results, method='max'):
    """
    Decide the color to use for a point in parameter space (Fig 3 of the
    original publication). Uses `classify_run` to classify each point, and
    returns a color based on the most common classification for the repeated
    runs that share the same parameters.
    """
    classifications = [classify_run(max_f_before,
                                    min_g_size,
                                    max_f_after)
                       for max_f_before, min_g_size, max_f_after in results]
    if method == 'max':
        # Use the classification that came up the most
        classification = Counter(classifications).most_common(1)[0][0]

        if classification == 'UNSTABLE_BEFORE':
            return 1, 0, 0  # red
        elif classification == 'UNSTABLE_AFTER':
            return 1, 1, 0  # yellow
        elif classification == 'PERSISTENT':
            return 0, 0, 1  # blue
        elif classification == 'SHORT-LIVED':
            return 0, 1, 0  # green
        else:
            raise AssertionError('Unknown classification: %s' % classification)
    elif method == 'mean':
        total = len(results)
        # Fractions for the different outcomes
        unstable_before = sum(c == 'UNSTABLE_BEFORE' for c in classifications)
        unstable_after = sum(c == 'UNSTABLE_AFTER' for c in classifications)
        short_lived = sum(c == 'SHORT-LIVED' for c in classifications)
        persistent = sum(c == 'PERSISTENT' for c in classifications)
        assert unstable_before + unstable_after + short_lived + persistent == total
        return np.array([unstable_before + unstable_after,
                        short_lived + unstable_after,
                        persistent], dtype=float)/total
    else:
        raise AssertionError('Unkown method: %s' % method)
