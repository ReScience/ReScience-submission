# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Pamela Hathway p.hathway16@imperial.ac.uk
# ----------------------------------------------------------------------------- 
# Creates figure 2
# This script can be run independently.
# -----------------------------------------------------------------------------

from brian2 import *
import time


def fig_2():
    print('%s Preparing Figure 2: STDP rules' % time.strftime('%H:%M'))
    plt.rcParams['savefig.dpi'] = 75
    plt.rcParams['figure.autolayout'] = False
    plt.rcParams['figure.figsize'] = 12, 5
    plt.rcParams['axes.labelsize'] = 17
    plt.rcParams['axes.titlesize'] = 17
    plt.rcParams['font.size'] = 17
    plt.rcParams['lines.linewidth'] = 2.0
    plt.rcParams['lines.markersize'] = 3
    plt.rcParams['legend.fontsize'] = 16

    device.reinit()
    device.activate()

    # spike times
    times_pre = array([2, 4, 6, 14, 16, 18, 22, 24, 26]) * ms
    indices_pre = zeros(len(times_pre))
    times_post = array([8, 12, 20]) * ms
    indices_post = zeros(len(times_post))
    factor = array([3])  # how long the interval between spikes should be

    # parameters STDP synapses
    tauplus = 16.8 * ms
    tauminus = 33.7 * ms
    aplus = 2.0 ** -5  # 0.03125# 0.0008#0.03125
    aminus = 0.85 * aplus
    defaultclock.dt = 1e-5 * second  # length of time step for simulation

    # equations synapses
    syn_eqs = ''' wi : 1
	dLTPtrace/dt = -LTPtrace / tauplus : 1 (clock-driven)
	dLTDtrace/dt = -LTDtrace / tauminus : 1 (clock-driven)'''

    # equations ATA rule
    synATA_pre = '''LTPtrace += aplus
	wi = clip(wi + LTDtrace, 0, 1)'''

    synATA_post = ''' LTDtrace -= aminus
	wi = clip(wi + LTPtrace, 0, 1)'''

    # equations NN rule
    synNSme_pre = '''LTPtrace = aplus
	wi = clip(wi + LTDtrace, 0, 1)'''

    synNSme_post = ''' LTDtrace = -aminus
	wi = clip(wi + LTPtrace, 0, 1)'''

    # equations RNN
    synRNN_pre = '''LTPtrace = aplus
	wi = clip(wi + LTDtrace, 0, 1)
	LTDtrace = 0'''

    synRNN_post = ''' LTDtrace = -aminus
	wi = clip(wi + LTPtrace, 0, 1)
	LTPtrace = 0'''

    Nin = SpikeGeneratorGroup(1, indices_pre, times_pre * factor)
    Nout = SpikeGeneratorGroup(1, indices_post, times_post * factor)

    synATA = Synapses(Nin, Nout, model=syn_eqs, on_pre=synATA_pre, on_post=synATA_post, method='linear')
    synATA.connect(i=0, j=0)
    synATA.wi = 0.5
    synATAstm = StateMonitor(synATA, ['wi', 'LTPtrace', 'LTDtrace'], record=True, dt=0.01 * ms)

    synNSme = Synapses(Nin, Nout, model=syn_eqs, on_pre=synNSme_pre, on_post=synNSme_post, method='linear')
    synNSme.connect(i=0, j=0)
    synNSme.wi = 0.5
    synNSmestm = StateMonitor(synNSme, ['wi', 'LTPtrace', 'LTDtrace'], record=True, dt=0.01 * ms)

    synRNN = Synapses(Nin, Nout, model=syn_eqs, on_pre=synRNN_pre, on_post=synRNN_post, method='linear')
    synRNN.connect(i=0, j=0)
    synRNN.wi = 0.5
    synRNNstm = StateMonitor(synRNN, ['wi', 'LTPtrace', 'LTDtrace'], record=True, dt=0.01 * ms)

    run((times_pre[-1] + 1 * ms) * factor)

    # parameters for plotting
    height_post = 0.65
    height_pre = 0.7

    fig2 = figure()
    suptitle('Figure 2')
    subplot(1, 3, 1)
    title('ATA')
    plot(times_pre / ms * factor / 1000, height_pre + indices_pre, 'g.', markersize=8, label='pre')
    plot(times_post / ms * factor / 1000, height_post + indices_post, 'b.', markersize=8, label='post')
    plot(synATAstm.t_, synATAstm.wi[0], 'k', label='weight')
    ylim([0.4, height_pre + 0.07])
    y = [0.5, height_post, height_pre]
    labels = ['weight', 'spikes post', 'spikes pre']
    yticks(y, labels)
    xticks([])

    subplot(1, 3, 2)
    title('NN')
    plot(times_pre / ms * factor / 1000, height_pre + indices_pre, 'g.', markersize=8)
    plot(times_post / ms * factor / 1000, height_post + indices_post, 'b.', markersize=8)
    plot(synNSmestm.t_, synNSmestm.wi[0], 'k')
    ylim([0.4, height_pre + 0.07])
    yticks([])
    xticks([])

    subplot(1, 3, 3)
    title('RNN')
    plot(times_pre / ms * factor / 1000, height_pre + indices_pre, 'g.', markersize=8, label='pre')
    plot(times_post / ms * factor / 1000, height_post + indices_post, 'b.', markersize=8, label='post')
    plot(synRNNstm.t_, synRNNstm.wi[0], 'k', label='weight')
    ylim([0.4, height_pre + 0.07])
    yticks([])
    xticks([])

    tight_layout(rect=[0, 0.03, 1, 0.95])

    savefig('../article/figures/figure_2_created_%s.pdf' % datetime.datetime.now().strftime('%y%m%d'))

    return fig2


if __name__ == '__main__':
    fig2 = fig_2()
    show()
