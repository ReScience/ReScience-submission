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
    ax_raster.plot(spikes_t/ms, spikes_i, '|', color="black")
    red_spikes = ((spikes_t >= 150*ms) &
                  ((spikes_t % (5*ms) < dt / 2) |
                   (spikes_t % (5*ms) > (5*ms - dt / 2))))

    ax_raster.plot(spikes_t[red_spikes]/ms, spikes_i[red_spikes], 'r|')

    ax_raster.set(xlim=(110, 220), ylim=(0, 200), xlabel='Time (ms)',
                  ylabel='Neuron')

    bins = np.arange(0, 250 + dt/ms, dt/ms)
    binned_activity, _ = np.histogram(spikes_t / ms, bins)
    # Select the spikes that (most likely) arise from the synchronous
    # activation of cells at time stim_time: Cells spiking *exactly* a multiple
    # of synaptic_delay after the initial activity will be assumed to having
    # been triggered by it.
    delay_steps = int(round(5*ms/dt))
    group_size = binned_activity[int(round(150*ms/dt))::delay_steps]

    # Draw the population rate
    ax_binned.bar(bins[:-1], binned_activity, color="black")
    ax_binned.set(ylim=(0, ymax1), ylabel='Rate')

    # Draw the group size
    ax_groups.bar(np.arange(150, 250, 5), group_size, color='red')
    ax_groups.set(ylim=(0, ymax1), ylabel="g'")

    save_fig(save, fig)


def plot_grid(par_range, colors, save=None):
    """Draw the result of a parameter search"""
    fig, ax = plt.subplots(gridspec_kw={'top': 0.95,
                                        'right': 0.95,
                                        'left': 0.175,
                                        'bottom': 0.15})
    ax.imshow(colors, extent=[min(par_range/mV),
                              max(par_range/mV),
                              max(par_range/mV),
                              min(par_range/mV)])
    ax.set_xlabel("Total ext weight")
    ax.set_ylabel("Total Inh weight")
    save_fig(save, fig=fig)


def plot_markov(group_x, group_ev, analytic_x, semi_analytic, save=None, linear=True):
    """Plot the most probable evolution of a group size at time t+1 given the
    size at time t"""
    fig, ax = plt.subplots(gridspec_kw={'top': 0.95,
                                        'bottom': 0.15,
                                        'left': 0.2,
                                        'right': 0.95})
    # Numeric solution
    ax.set(aspect='equal', xticks=[1] + list(range(25, group_x[-1]+1, 25)),
       yticks=[1] + list(range(25, group_x[-1]+1, 25)),
       xlim=(0, group_x[-1]), ylim=(0, group_x[-1]),
       xlabel="# of synchronized neurons g'0",
       ylabel="# of synchronized neurons g'1")
    ax.plot(group_x, np.array(group_ev), marker='_',
            ms=3.5,
            mew=1,
            color="black",
            linestyle='', alpha=0.2)
    ax.plot(group_x,
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
