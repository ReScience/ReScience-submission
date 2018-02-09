from __future__ import print_function

from itertools import product
import os

import numpy as np
import brian2 as br2
from brian2 import mV, us, ms
from joblib import Parallel, delayed
import matplotlib.pyplot as plt

from lib import run_with_cache, group_size_evolutions, calc_group_evolution_cached, decide_color, grid_search
from plot import plot_population_activity, plot_markov, plot_grid

# Note that 'weave', Brian's standard target on Python 2 is not compatible with
# multiprocessing and can therefore not be used (we want to do parameter
# explorations in parallel)
br2.prefs.codegen.target = 'cython'

FIG_FOLDER = "../article/Figs/"

main_seeds = [3982293089, 1861487226, 2570878100, 3472600360, 2780156760]

if not os.path.exists(FIG_FOLDER):
    os.makedirs(FIG_FOLDER)

# Set figure options from configuration file
plt.style.use('figures.conf')

def generate_seeds(main_seed, n_seeds):
    np.random.seed(main_seed)
    return np.random.randint(0, 2**32-1, size=n_seeds, dtype=np.uint32)


def fig2(quiet=False, dt=0.1*ms, seed=None):
    """Generate the raster plot and population statistics for Figure 2
    of the original paper (https://doi.org/10.1371/journal.pcbi.1002384.g002)"""
    if not quiet:
        print('Running simulations for Figure 2')
    seeds = generate_seeds(seed, 2)
    results = run_with_cache(w_ex=0.2*mV, w_in=0.2*mV, linear=True, dt=dt, seed=seeds[0])
    results_nl = run_with_cache(w_ex=0.2*mV, w_in=0.2*mV, linear=False, dt=dt, seed=seeds[1])

    figname = "spikes_l_dt_%.0fus.pdf" % (dt/us)
    plot_population_activity(results, save=FIG_FOLDER + figname, dt=dt)
    figname = "spikes_nl_dt_%.0fus.pdf" % (dt/us)
    plot_population_activity(results_nl, save=FIG_FOLDER + figname, dt=dt)


def fig3(linear=True, short=True, quiet=False, dt=0.1*ms, seed=None):
    """Perform a grid search over the excitatory and inhibitory coupling
    strengths and generate Figures 3a and 3b of the original paper
    (https://doi.org/10.1371/journal.pcbi.1002384.g003)"""
    if not quiet:
        print('Running simulations for Figure 3 ({} '
              'coupling)'.format('linear' if linear else 'nonlinear'))

    if short:
        par_range = np.linspace(0.16, 0.4, 10, endpoint=True)*mV
        repetitions = 5
        name = "grid_short"
    else:
        par_range = np.linspace(0.16, 0.4, 100, endpoint=True)*mV
        repetitions = 10
        name = "grid"

    seeds = generate_seeds(seed, len(par_range)**2*repetitions)

    # Generate the data
    grid_results = grid_search(par_range, linear=linear, n_rep=repetitions,
                               dt=dt, seeds=seeds)
    # Plot
    for method, method_string in zip(['max', 'mean'], ['', '_mean']):
        if linear:
            figname = name + "_l_dt_%.0fus%s.pdf" % (dt/us, method_string)
        else:
            figname = name + "_nl_dt_%.0fus%s.pdf" % (dt/us, method_string)
        colors = np.zeros((len(par_range), len(par_range), 3))
        for i, j in product(range(len(par_range)), range(len(par_range))):
            c = decide_color(grid_results[par_range[i]/mV, par_range[j]/mV],
                             method=method)
            colors[j][i] = np.array(c)
        plot_grid(par_range * 150, colors, save=FIG_FOLDER + figname)


def fig4(linear=True, short=True, quiet=False, dt=0.1*ms, seed=None):
    """Simulate linearly and nonlinearly coupled networks for pulses of
     different sizes. Generates Figure 4 of the original paper
     (https://doi.org/10.1371/journal.pcbi.1002384.g004)"""
    if not quiet:
        print('Running simulations for Figure 4 ({} '
              'coupling)'.format('linear' if linear else 'nonlinear'))

    # Generate the data
    if short:
        n_rep = 5
        step_size = 5
    else:
        n_rep = 50
        step_size = 1

    gsize_x = np.arange(1, 180, 5)
    analytic_x = np.arange(1, 180, step_size)
    seeds = generate_seeds(seed, len(gsize_x)*n_rep + n_rep)
    # Numerical results
    gsize = Parallel(n_jobs=-2)(delayed(group_size_evolutions)(i, n_rep,
                                                               linear=linear,
                                                               dt=dt,
                                                               seeds=seeds[idx*n_rep:(idx+1)*n_rep])
                                for idx, i in enumerate(gsize_x))
    # Semi-analytical results
    seed_start = len(gsize_x)*n_rep
    result_list = Parallel(n_jobs=-2)(delayed(run_with_cache)(stim_time=None,
                                                              linear=linear,
                                                              dt=dt,
                                                              store_v=True,
                                                              repetition=i,
                                                              seed=seeds[seed_start+i])
                                     for i in range(1, n_rep))
    bins = result_list[0]['bins']
    v_hist = np.sum(np.array([r['v'] for r in result_list],dtype=float), axis=0)
    v_hist /= np.sum(v_hist)  # Normalize the sum to 1
    if not quiet:
        print('\tCalculating semi-analytical solution')
    semi_analytic = Parallel(n_jobs=-2)(delayed(calc_group_evolution_cached)(bins, v_hist, i, linear=linear)
                                        for i in analytic_x)
    # Plot the data
    if linear:
        fname = "markov_l_dt_%.0fus.pdf" % (dt/us)
    else:
        fname = "markov_nl_dt_%.0fus.pdf" % (dt/us)

    plot_markov(gsize_x, gsize, analytic_x, semi_analytic, linear=linear, save=FIG_FOLDER + fname)


def all_figs(short=True, quiet=False, dt=0.1*ms):
    """Generate all figures. With ``short=True`` (the default), coarser
    parameter explorations and less repetitions are used compared to the
    original paper, to reduce the total simulation time."""
    fig2(quiet=quiet, dt=dt, seed=main_seeds[0])

    fig3(linear=False, short=short, quiet=quiet, dt=dt, seed=main_seeds[1])
    fig3(linear=True, short=short, quiet=quiet, dt=dt, seed=main_seeds[2])

    fig4(linear=False, short=short, quiet=quiet, dt=dt, seed=main_seeds[3])
    fig4(linear=True, short=short, quiet=quiet, dt=dt, seed=main_seeds[4])


if __name__ == "__main__":
    import time
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--long',
                        help='Run the full set of parameters/repetitions',
                        action='store_true')
    parser.add_argument('-q', '--quiet',
                        help='Do not print progress information',
                        action='store_true')
    parser.add_argument('--dt', type=float,
                        help='Simulation time step to use (in ms). '
                             'Defaults to 0.1ms',
                        default=0.1)
    args = parser.parse_args()
    short = not args.long
    quiet = args.quiet
    dt = args.dt*ms
    if not quiet:
        if short:
            print('Running a reduced number of parameters/repetitions.')
        else:
            print('Running the full number of parameters/repetitions.')
            print('These simulations will take several hours to complete. Storing '
                  'the results will take about ~6GB of disk space.')
        print('Note that you can interrupt the simulation at any time, intermediate '
              'results will be stored and reused the next time.')
        print('Using a simulation time step of %s.' % str(dt))
    start = time.time()
    all_figs(short=short, quiet=quiet, dt=dt)
    took = time.time() - start
    hours = int(took / 60. / 60.)
    minutes = int(60*(took / 60. / 60. - hours))
    if not quiet:
        print('Finished after %dh%dm' % (hours, minutes))
