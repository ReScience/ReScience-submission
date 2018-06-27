# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Pamela Hathway p.hathway16@imperial.ac.uk
# ----------------------------------------------------------------------------- 
# Functions to create figures 3, 4, 7AB
# This script cannot be run independently but is called by run_simulation.py
# -----------------------------------------------------------------------------

from brian2 import *


def fig_3(N2spm, latency, random_seed):
    plt.rcParams['savefig.dpi'] = 75
    plt.rcParams['figure.autolayout'] = False
    plt.rcParams['figure.figsize'] = 6, 5
    plt.rcParams['axes.labelsize'] = 17
    plt.rcParams['axes.titlesize'] = 17
    plt.rcParams['font.size'] = 17
    plt.rcParams['lines.linewidth'] = 2.0
    plt.rcParams['lines.markersize'] = 8
    plt.rcParams['legend.fontsize'] = 16
    fig3 = figure()
    suptitle('Figure 3')
    title('latency vs output neuron spikes')
    plot(range(len(latency)), latency * 1000, 'g.', markersize=3)
    ylabel('Postsynaptic spike latency [ms]')
    xlabel('# discharge')
    yticks([0, 10, 20, 30, 40, 50])
    ylim([-2, 52])

    tight_layout(rect=[0, 0.03, 1, 0.95])

    savefig('../article/figures/figure_3_seed%s_created_%s.pdf' % (random_seed, datetime.datetime.now().strftime('%y%m%d')))

    return fig3


def fig_4(N1spm, syn12, timing_pattern, random_seed):
    plt.rcParams['savefig.dpi'] = 75
    plt.rcParams['figure.autolayout'] = False
    plt.rcParams['figure.figsize'] = 12, 5
    plt.rcParams['axes.labelsize'] = 17
    plt.rcParams['axes.titlesize'] = 17
    plt.rcParams['font.size'] = 17
    plt.rcParams['lines.linewidth'] = 2.0
    plt.rcParams['lines.markersize'] = 8
    plt.rcParams['legend.fontsize'] = 16
    # get spikes from end of simulation
    figsta = searchsorted(N1spm.t / second, timing_pattern[-2] - 0.025)
    figsto = searchsorted(N1spm.t / second, timing_pattern[-2] + 0.075)
    greyness = [str((255 - syn12.wi[int(i)]*255) / 255.) for i in N1spm.i[figsta:figsto]]  # high w == white

    ax = figure(figsize=(8, 6)).add_subplot(1, 1, 1)
    ax.set_facecolor("k")
    sc = scatter(N1spm.t[figsta:figsto] / second, N1spm.i[figsta:figsto], c=greyness, s=3, lw=0, cmap='Greys')
    cbar = colorbar(sc)
    cbar.ax.invert_yaxis()
    cbar.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    cbar.set_ticklabels([1, 0.8, 0.6, 0.4, 0.2, 0])
    cbar.set_label('synaptic weight')

    title('Figure 4')
    xlabel('time [s]')
    ylabel('# afferent')
    xlim([timing_pattern[-2] - 0.025, timing_pattern[-2] + 0.075])
    ylim([0, 2000])
    xticks(array([timing_pattern[-2], timing_pattern[-2] + 0.05]))
    a = ax.get_xticks().tolist()
    a[0] = timing_pattern[-2]
    a[1] = timing_pattern[-2] + 0.05
    ax.set_xticklabels(a)
    ax.add_patch(
        matplotlib.patches.Rectangle((timing_pattern[-2] - 0.001, 0), 0.05, 1000, linewidth=3, color='b', fill=False))

    tight_layout(rect=[0, 0.03, 1, 0.95])

    savefig('../article/figures/figure_4_seed%s_created_%s.pdf' % (random_seed, datetime.datetime.now().strftime('%y%m%d')))

    return ax


def fig_7AB(syn12stm, syn12nnstm, syn12atastm, record_wi_dt, random_seed):
    fig = figure(figsize=(12, 5))
    suptitle('Figure 7 A and B')
    subplot(121)
    title('first second')

    plot(arange(0, len(syn12stm.wi[0, :]) * record_wi_dt, record_wi_dt), mean(syn12stm.wi, axis=0), 'k', label='RNN')
    plot(arange(0, len(syn12nnstm.wi[0, :]) * record_wi_dt, record_wi_dt), mean(syn12nnstm.wi, axis=0), 'k--',
         label='NN')
    plot(arange(0, len(syn12atastm.wi[0, :]) * record_wi_dt, record_wi_dt), mean(syn12atastm.wi, axis=0), 'k-.',
         label='ATA')

    ylabel('average weight per synapse')
    xlabel('time [s]')
    ylim([0.15, 0.55])
    xticks([0, 0.5, 1])
    yticks(array([0.2, 0.3, 0.4, 0.5]))
    xlim([0, 1])

    subplot(122)
    title('whole simulation')
    plot(arange(0, len(syn12stm.wi[0, :]) * record_wi_dt, record_wi_dt), mean(syn12stm.wi, axis=0), 'k', label='RNN')
    plot(arange(0, len(syn12nnstm.wi[0, :]) * record_wi_dt, record_wi_dt), mean(syn12nnstm.wi, axis=0), 'k--',
         label='NN')
    plot(arange(0, len(syn12atastm.wi[0, :]) * record_wi_dt, record_wi_dt), mean(syn12atastm.wi, axis=0), 'k-.',
         label='ATA')

    xlabel('time [s]')
    ylim([0.15, 0.55])
    xticks(array([1, 150, 300, 450]))
    yticks(array([0.2, 0.3, 0.4, 0.5]))
    xlim([1, 450])
    legend(loc=1)

    tight_layout(rect=[0, 0.03, 1, 0.95])

    savefig('../article/figures/figure_7AB_seed%s_created_%s.pdf' % (random_seed, datetime.datetime.now().strftime('%y%m%d')))

    return fig


def fig_9(N2stm, timing_pattern, find_t, record_u_dt, random_seed):
    plt.rcParams['lines.linewidth'] = 1.0
    xlim_starts = [0, find_t - 0.25, 449]

    fig = figure(figsize=(12, 7))
    suptitle('Figure 9 (supplemental only, Figure 4 original)')

    for subplot_idx in range(3):
        subplot(3, 1, subplot_idx + 1)
        plot(N2stm.t[int(xlim_starts[subplot_idx]/record_u_dt):int((xlim_starts[subplot_idx]+1)/record_u_dt)] / second,
             N2stm.u[0, int(xlim_starts[subplot_idx]/record_u_dt):int((xlim_starts[subplot_idx]+1)/record_u_dt)], 'b', label='potential')
        axhline(y=500, color='r', linestyle='--', label='threshold')
        axhline(y=0, color='k', linestyle='--', label='resting pot.')

        start_in_timingpattern = where(timing_pattern >= xlim_starts[subplot_idx])[0][0]
        for pattime in timing_pattern[start_in_timingpattern:start_in_timingpattern+10]:
            plt.axvline(pattime, ls='--', c='r', linewidth=2)
            plt.axvspan(pattime, (pattime + 0.05), facecolor='g', alpha=0.3)
        xlim([xlim_starts[subplot_idx], xlim_starts[subplot_idx] + 1])
        yticks([0, 500, 1000])
        ylim([-400, 1100])
        ylabel('Potential [a.u.]')
        xlabel('t [s]')

    legend(loc=1)

    savefig('../article/figures/figure_9sup_seed%s_created_%s.pdf' % (random_seed, datetime.datetime.now().strftime('%y%m%d')))


    return fig


def fig_10(N1spm, timing_pattern, patternlength, random_seed):
    import matplotlib.gridspec as gridspec
    plt.rcParams['lines.markersize'] = 3

    timeshown = 0.5
    starttime = timing_pattern[0] - patternlength
    stoptime = starttime + timeshown

    n_inds2plot = 50

    plot_from = searchsorted(N1spm.t / second, starttime)
    plot_until = searchsorted(N1spm.t / second, stoptime)
    times = N1spm.t[plot_from:plot_until] / second
    inds = N1spm.i[plot_from:plot_until]

    indices_pat50 = inds[where(inds < n_inds2plot)[0]]
    times_pat50 = times[where(inds < n_inds2plot)[0]]
    # indices_pat50 = indices_pat50 * 2 * (indices_pat50 > 0)

    indices_npat50 = inds[where((inds >= 1000) & (inds < 1000 + n_inds2plot))[0]]
    times_npat50 = times[where((inds >= 1000) & (inds < 1000 + n_inds2plot))[0]]
    indices_npat50 = (indices_npat50 - 1000) + n_inds2plot


    fig = figure(figsize=(8, 6))
    gridspec.GridSpec(4, 4)
    suptitle('Figure 10 (supplemental only, original Figure 1)')

    subplot2grid((4, 4), (0, 0), colspan=3, rowspan=3)
    plot(times_pat50, indices_pat50, 'b.')
    plot(times_npat50, indices_npat50, 'b.')
    for pattime in timing_pattern[0:where(timing_pattern < stoptime)[0][-1] + 1]:
        start = where(times_pat50 >= pattime)[0][0]
        stop = where(times_pat50 < pattime + patternlength)[0][-1] + 1
        plot(times_pat50[start:stop], indices_pat50[start:stop], 'r.')
    xticks([])
    ylabel('# afferent')
    xlim([starttime, stoptime])
    ylim([-1, 100])

    subplot2grid((4, 4), (0, 3), colspan=1, rowspan=3)
    uniq, counts = unique(indices_pat50, return_counts=True)
    barh(uniq + 0.5, counts*2, color='r')
    uniq, counts = unique(indices_npat50, return_counts=True)
    barh(uniq + 0.5, counts*2, color='b')
    xlabel('Firing rate [Hz]')
    ylim([0, 100])
    xticks([0, 50, 100])
    yticks([])

    subplot2grid((4, 4), (3, 0), colspan=3, rowspan=1)
    times_100 = concatenate((times_pat50, times_npat50))
    counts_hist, bins = histogram(times_100, 50)
    bar(arange(0.5, 50.5, 1), counts_hist, color='b')
    xlabel('t [ms]')
    ylabel('Firing rate [Hz]')
    xlim([0, 50])
    ylim([0, 100])

    savefig('../article/figures/figure_10sup_seed%s_created_%s.pdf' % (random_seed, datetime.datetime.now().strftime('%y%m%d')))

    return fig




if __name__ == '__main__':
    print('These figures rely on output from run_simulation.py. '
          'Please run that and these figures will be generated automatically.')
