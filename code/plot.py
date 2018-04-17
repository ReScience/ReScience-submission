import matplotlib.pyplot as plt
import numpy as np

from brian2 import ms, mV


def save_fig(save=None, fig=None):
    """Code snippet to use at the end of making a figure"""
    if save:
        if fig is not None:
            fig.savefig(save)
        else:
            plt.savefig(save)
        plt.close(fig)
    else:
        plt.show()


def plot_population_activity(results, ymax1=250, save=None, dt=0.1*ms):
    """Draw and/or save the figure"""
    spikes_i, spikes_t = results['i'], results['t']
    fig, (ax_groups, ax_binned, ax_raster) = plt.subplots(3, 1,
                                                          gridspec_kw={'height_ratios': (1, 1, 3),
                                                                       'left': 0.2,
                                                                       'bottom': 0.15,
                                                                       'top': 0.95,
                                                                       'right': 0.95},
                                                          sharex='col')
    plt.figure()
    # Draw the raster plot
    ax_raster.vlines(np.arange(150, 250, 5), 0, 200, color="gray", linestyle='dashed')
    ax_raster.plot(spikes_t/ms, spikes_i + 1, '|', color="black")  # Plot a 1-based neuron index
    red_spikes = ((spikes_t >= 150*ms) &
                  ((spikes_t % (5*ms) < dt / 2) |
                   (spikes_t % (5*ms) > (5*ms - dt / 2))))

    ax_raster.plot(spikes_t[red_spikes]/ms, spikes_i[red_spikes] + 1, 'r|')

    ax_raster.set(xlim=(110, 220), ylim=(0, 200), xlabel='Time (ms)',
                  ylabel='Neuron', xticks=[120, 150, 180, 210],
                  yticks=[1, 100, 200])

    # Note that we use the simulation dt as the bin size here, so that we don't
    # lose any time resolution when checking for spikes triggered by the
    # synchronous stimulation. When plotting the population rate, we however
    # bin into bins of 1ms
    bins = np.arange(0, 250 + dt/ms, dt/ms)
    binned_activity, _ = np.histogram(spikes_t / ms, bins)
    # Select the spikes that (most likely) arise from the synchronous
    # activation of cells at time stim_time: Cells spiking *exactly* a multiple
    # of synaptic_delay after the initial activity will be assumed to having
    # been triggered by it.
    delay_steps = int(round(5*ms/dt))
    group_size = binned_activity[int(round(150*ms/dt))::delay_steps]

    # Draw the population rate (with 1ms bins)
    steps_per_bin = int(round(1*ms/dt))
    binned_activity_1ms = binned_activity.reshape(-1, steps_per_bin).sum(axis=1)
    ax_binned.bar(bins[:-1:steps_per_bin], binned_activity_1ms,
                  color="black")
    ax_binned.set(ylim=(0, ymax1), ylabel='Rate', yticks=[0, 100, 200])

    # Draw the group size
    ax_groups.bar(np.arange(150, 250, 5), group_size, color='red')
    ax_groups.set(ylim=(0, ymax1), ylabel="g'", yticks=[0, 100, 200])

    save_fig(save, fig)


def plot_grid(par_range_exc, par_range_inh, colors, show_ylabel=True, save=None):
    """Draw the result of a parameter search"""
    fig, ax = plt.subplots(figsize=(1.72, 1.72),  # three panels
                           gridspec_kw={'top': 0.95,
                                        'right': 0.95,
                                        'left': 0.25,
                                        'bottom': 0.23})
    ax.imshow(colors, extent=[min(par_range_exc/mV),
                              max(par_range_exc/mV),
                              max(par_range_inh/mV),
                              min(par_range_inh/mV)])
    if par_range_exc[-1] - par_range_exc[0] > 30*mV:  # full exploration
        xticks = [30, 45, 60]
        yticks = [30, 45, 60]
    else:  # exploration of small zone
        xticks = [35, 40, 45]
        yticks = [25, 30]
    ax.set(xlabel="Total excitatory weight", xticks=xticks, yticks=yticks)
    if show_ylabel:
        ax.set_ylabel("Total inhibitory weight")
    save_fig(save, fig=fig)


def plot_markov(group_x, group_ev, analytic_x, semi_analytic, save=None,
                linear=True, bias_correction=False):
    """Plot the most probable evolution of a group size at time t+1 given the
    size at time t"""
    fig, ax = plt.subplots(gridspec_kw={'top': 0.95,
                                        'bottom': 0.15,
                                        'left': 0.2,
                                        'right': 0.95})
    # Numeric solution
    g_symbol = '\hat{g}' if bias_correction else 'g'
    ax.set(aspect='equal', xticks=[1] + list(range(25, np.max(group_x)+1, 25)),
           yticks=list(range(0, np.max(group_x)+1, 25)),
           xlim=(-5, np.max(group_x)+3.5/2), ylim=(-5, np.max(group_x)+3.5/2),
           xlabel="# of synchronized neurons ${}'_0$".format(g_symbol),
           ylabel="# of synchronized neurons ${}'_1$".format(g_symbol))
    ax.plot(group_x, np.array(group_ev), marker='_',
            ms=3.5,
            mew=1,
            color="black",
            linestyle='', alpha=0.2)
    ax.plot(np.mean(np.array(group_x), axis=1),
            np.mean(np.array(group_ev), axis=1),
            marker='s', color="lime",
            linestyle='')
    ax.plot(group_x, group_x, color="black")
    ax.tick_params(direction='in', pad=5)

    # (Semi-)analytic solution
    ax.plot(analytic_x, semi_analytic, 'r:', lw=2)

    # Making an inset
    inset = fig.add_axes([0.25, 0.7, .2, .2])
    if linear:
        inset.plot(range(11), range(11), color='black')
    else:
        inset.plot(range(11), range(11), color='gray')
        inset.plot(range(11),
                   list(range(3)) + [4, 6] + [6 for i in range(6)],
                   color='black')
    inset.set(xlim=(0, 10), ylim=(0, 10), xticks=[], yticks=[])

    save_fig(save, fig=fig)
