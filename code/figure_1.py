# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Pamela Hathway p.hathway16@imperial.ac.uk
# ----------------------------------------------------------------------------- 
# Creates figure 1
# This script can be run independently.
# -----------------------------------------------------------------------------

from brian2 import *
import time


def fig_1():
    print('%s Preparing Figure 1: Potentials' % time.strftime('%H:%M'))
    plt.rcParams['savefig.dpi'] = 75
    plt.rcParams['figure.autolayout'] = False
    plt.rcParams['figure.figsize'] = 6, 5
    plt.rcParams['axes.labelsize'] = 17
    plt.rcParams['axes.titlesize'] = 17
    plt.rcParams['font.size'] = 17
    plt.rcParams['lines.linewidth'] = 2.0
    plt.rcParams['lines.markersize'] = 3
    plt.rcParams['legend.fontsize'] = 16

    device.reinit()
    device.activate()

    runduration = 0.08
    defaultclock.dt = 1e-4 * second

    # output neuron
    win = 1  # initial weight
    taum = 10 * ms
    taus = 2.5 * ms
    tausyn = 2.5 * ms
    T = 2.8
    u_rest = 0
    X = (taus/taum)**(taum/(taus - taum))
    A = -8.5
    deltax = 1
    deltaa = 1

    # equations
    eqs = '''du/dt = (A*a)/taus + (X*x-u)/taum : 1
                dx/dt = -x/tausyn : 1 
                da/dt = -a/taus : 1'''
    eqs_reset = '''x = 0
                u = 2*T
                a = deltaa'''
    stdp_model = ''' wi : 1'''
    stdp_on_pre = ''' x_post += deltax*wi '''

    # input
    times = array([5, 20, 40, 43, 46, 60]) * ms
    indices = zeros(len(times))

    N0 = SpikeGeneratorGroup(1, indices, times)
    N1 = NeuronGroup(1, '''du/dt  = -u/taum : 1''', threshold='u > T', reset='u = u_rest', method='linear')
    N2 = NeuronGroup(1, eqs, threshold='u > T', reset=eqs_reset, refractory=2 * ms, method='linear')
    N2spm = SpikeMonitor(N2)
    N2stm = StateMonitor(N2, ['u'], record=True, dt=0.1 * ms)

    syn01 = Synapses(N0, N1, on_pre='''u_post += 10''', method='linear')
    syn01.connect(j='i')
    syn12 = Synapses(N1, N2, model=stdp_model, on_pre=stdp_on_pre, method='linear')
    syn12.connect(j='i')
    syn12.wi = win

    net = Network(collect())
    net.add(N0, N1, N2, N2spm, syn01, syn12)
    net.run(runduration * second)

    fig1 = figure()
    title('Figure 1')
    plot(N2stm.t / ms, N2stm.u[0, :], 'k', linewidth=2, label='potential')
    for q in times[:-1] / ms:
        axvline(x=q, color='k', linestyle='-.', linewidth=1)
    axvline(x=(times[-1] / ms), color='k', linestyle='-.', linewidth=1, label='input spikes')
    axhline(y=T, color='r', linestyle='--', linewidth=1, label='threshold')
    axhline(y=0, color='k', linestyle=':', linewidth=1, label='resting potential')
    xlabel('times [ms]')
    ylabel('voltage [a.u.]')
    xlim([0, 80])

    tight_layout()

    savefig('../article/figures/figure_1_created_%s.pdf' % datetime.datetime.now().strftime('%y%m%d'))

    return fig1


if __name__=='__main__':
    fig1 = fig_1()
    show()