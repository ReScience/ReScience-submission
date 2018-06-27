# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Pamela Hathway p.hathway16@imperial.ac.uk
# ----------------------------------------------------------------------------- 
# Creates figure 5 either from scratch (default) or using saved data (to switch,
# comment out lines 225-226 and comment in line 224). 
# This script can be run independently.
# -----------------------------------------------------------------------------

from brian2 import *
import matplotlib.pyplot as plt
import pickle as pickle
import time



def fig_5_saved():

    import pickle as pickle
    print('%s Preparing Figure 5: Robustness' % time.strftime('%H:%M'))
    plt.rcParams['savefig.dpi'] = 75
    plt.rcParams['figure.autolayout'] = False
    plt.rcParams['figure.figsize'] = 8, 5
    plt.rcParams['axes.labelsize'] = 17
    plt.rcParams['axes.titlesize'] = 17
    plt.rcParams['font.size'] = 17
    plt.rcParams['lines.linewidth'] = 2.0
    plt.rcParams['lines.markersize'] = 3
    plt.rcParams['legend.fontsize'] = 13

    # ### load parameter arrays and original paper results
    with open('../data/results_masquelier2008.pickle', 'rb') as handle:
            results_masq = pickle.load(handle)

    # ### load results from different parameter conditions
    runs_to_plot = ['e4', 'e6']
    names_to_plot = ['rep 10$^{-4}$', 'rep 10$^{-6}$']
    styles = ['-', '--']
    colors = ['k', 'k']

    fig = figure(figsize=(12, 4))
    suptitle('Figure 5 saved results')

    for idx, run_name in enumerate(runs_to_plot):

        run_ids = arange(0, 22, 1)
        success_rates = empty(22)

        # save all success results in one array
        for run_idx in run_ids:

            with open('../data/results_saved_%s_%s.pickle' % (run_idx, run_name), 'rb') as handle:
                run = pickle.load(handle)

            success_rates[run_idx] = run['success_overall']

        # make separate arrays for each condition and add result from standard condition (0)
        win_res = concatenate((success_rates[1:5], array([success_rates[0]])))
        jit_res = concatenate((array([success_rates[5]]), array([success_rates[0]]), success_rates[6:11]))
        npat_res = concatenate((success_rates[11:14], array([success_rates[0]]), array([success_rates[14]])))
        freq_res = concatenate((success_rates[15:18], array([success_rates[0]]), array([success_rates[18]])))
        del_res = concatenate((array([success_rates[0]]), success_rates[19:22]))

        subplot(1, 5, 1)
        plot(results_masq['freq_x'], freq_res, marker='o', color=colors[idx], linestyle=styles[idx])

        subplot(1, 5, 2)
        plot(results_masq['jit_x'], jit_res, marker='o', color=colors[idx], linestyle=styles[idx])

        subplot(1, 5, 3)
        plot(results_masq['npat_x'], npat_res, marker='o', color=colors[idx], linestyle=styles[idx])

        subplot(1, 5, 4)
        plot(results_masq['win_x'], win_res, marker='o', color=colors[idx], linestyle=styles[idx], label=names_to_plot[idx])

        subplot(1, 5, 5)
        plot(results_masq['del_x'], del_res, marker='o', color=colors[idx], linestyle=styles[idx])

    subplot(1, 5, 1)
    plot(results_masq['freq_x'], results_masq['freq_masq'], marker='x', color='g', linestyle='-')
    ylabel('% of success')
    xlim([0, 0.55]), ylim([-5, 105])
    xlabel('Pattern freq')
    xticks(array([0.1, 0.3, 0.5]))
    yticks(array([0, 25, 50, 75, 100]))
    title('A', weight='bold')

    subplot(1, 5, 2)
    plot(results_masq['jit_x'], results_masq['jit_masq'], marker='x', color='g', linestyle='-')
    xlim([-0.5, 6.5]), ylim([-5, 105])
    xlabel('Jitter SD [ms]')
    yticks(array([0, 50, 100]))
    xticks(array([0, 2, 4, 6]))
    yticks(array([]))
    title('B', weight='bold')

    subplot(1, 5, 3)
    plot(results_masq['npat_x'], results_masq['npat_masq'], marker='x', color='g', linestyle='-')
    xlim([0.175, 0.625]), ylim([-5, 105])
    xlabel('Prop. pat. neur.')
    yticks(array([0, 50, 100]))
    xticks(array([0.2, 0.4, 0.6]))
    yticks(array([]))
    title('C', weight='bold')

    subplot(1, 5, 4)
    plot(results_masq['win_x'], results_masq['win_masq'], marker='x', color='g', linestyle='-', label='orig rerun')
    plot(array([-100, -200]), array([-5, -10]), marker='x', color='g', linestyle='--', label='orig paper')
    xlim([0.25, 0.5]), ylim([-5, 105])
    xlabel('Initial weight')
    yticks(array([0, 50, 100]))
    xticks(array([0.3, 0.4, 0.5]))
    yticks(array([]))
    title('D', weight='bold')
    legend(loc=3)

    subplot(1, 5, 5)
    plot(results_masq['del_x'], results_masq['del_masq'], marker='x', color='g', linestyle='--')
    xlim([-0.025, 0.325]), ylim([0 - 5, 105])
    xlabel('Spike deletion')
    yticks(array([0, 50, 100]))
    xticks(array([0, 0.1, 0.2, 0.3]))
    yticks(array([]))
    title('E', weight='bold')

    tight_layout(rect=[0, 0.03, 1, 0.95])

    savefig('../article/figures/figure_5_from_saved_created_%s.pdf' % datetime.datetime.now().strftime('%y%m%d'))

    return fig


def fig_5_new(reps_per_combination):
    # call run_sim 10 times with different resolutions and set runtime to 30seconds
    from run_simulation import run_sim

    param_combinations = 22
    success_rates = zeros(param_combinations)

    for pararow in range(param_combinations):

        for rep in range(0, reps_per_combination):

            print('    Parameter combination %s, repetition %s ' % (pararow, rep))
            success = run_sim(rep, pararow, only_success=True)
            success_rates[pararow] += success

    success_rates = success_rates/(reps_per_combination/100)

    # make separate arrays for each condition and insert result from standard condition (0)
    win_res = concatenate((success_rates[1:5], array([success_rates[0]])))
    jit_res = concatenate((array([success_rates[5]]), array([success_rates[0]]), success_rates[6:11]))
    npat_res = concatenate((success_rates[11:14], array([success_rates[0]]), array([success_rates[14]])))
    freq_res = concatenate((success_rates[15:18], array([success_rates[0]]), array([success_rates[18]])))
    del_res = concatenate((array([success_rates[0]]), success_rates[19:22]))

    with open('../data/results_masquelier2008.pickle', 'rb') as handle:
            results_masq = pickle.load(handle)

    fig = figure(figsize=(12, 4))
    suptitle('Figure 5')

    subplot(1, 5, 1)
    plot(results_masq['freq_x'], freq_res, marker='o', color='k', linestyle='-')
    plot(results_masq['freq_x'], results_masq['freq_masq'], marker='x', color='g', linestyle='-')
    ylabel('% of success')
    xlim([0, 0.55]), ylim([-5, 105])
    xlabel('Pattern freq')
    xticks(array([0.1, 0.3, 0.5]))
    yticks(array([0, 25, 50, 75, 100]))
    title('A', weight='bold')

    subplot(1, 5, 2)
    plot(results_masq['jit_x'], jit_res, marker='o', color='k', linestyle='-')
    plot(results_masq['jit_x'], results_masq['jit_masq'], marker='x', color='g', linestyle='-')
    xlim([-0.5, 6.5]), ylim([-5, 105])
    xlabel('Jitter SD [ms]')
    yticks(array([0, 50, 100]))
    xticks(array([0, 2, 4, 6]))
    yticks(array([]))
    title('B', weight='bold')

    subplot(1, 5, 3)
    plot(results_masq['npat_x'], npat_res, marker='o', color='k', linestyle='-')
    plot(results_masq['npat_x'], results_masq['npat_masq'], marker='x', color='g', linestyle='-')
    xlim([0.175, 0.625]), ylim([-5, 105])
    xlabel('Prop. pat. neur.')
    yticks(array([0, 50, 100]))
    xticks(array([0.2, 0.4, 0.6]))
    yticks(array([]))
    title('C', weight='bold')

    subplot(1, 5, 4)
    plot(results_masq['win_x'], win_res, marker='o', color='k', linestyle='-', label='rep 10$^{-4}$')
    plot(results_masq['win_x'], results_masq['win_masq'], marker='x', color='g', linestyle='-', label='orig rerun')
    plot(array([-100, -200]), array([-5, -10]), marker='x', color='g', linestyle='--', label='orig paper')
    xlim([0.25, 0.5]), ylim([-5, 105])
    xlabel('Initial weight')
    yticks(array([0, 50, 100]))
    xticks(array([0.3, 0.4, 0.5]))
    yticks(array([]))
    title('D', weight='bold')
    legend(loc=3)

    subplot(1, 5, 5)
    plot(results_masq['del_x'], del_res, marker='o', color='k', linestyle='-')
    plot(results_masq['del_x'], results_masq['del_masq'], marker='x', color='g', linestyle='--')
    xlim([-0.025, 0.325]), ylim([0 - 5, 105])
    xlabel('Spike deletion')
    yticks(array([0, 50, 100]))
    xticks(array([0, 0.1, 0.2, 0.3]))
    yticks(array([]))
    title('E', weight='bold')

    tight_layout(rect=[0, 0.03, 1, 0.95])

    savefig('../article/figures/figure_5_created_%s.pdf' % datetime.datetime.now().strftime('%y%m%d'))

    return fig


if __name__=='__main__':
    # fig_5_saved()
    reps_per_resolution = 10
    fig5 = fig_5_new(reps_per_resolution)

    show()