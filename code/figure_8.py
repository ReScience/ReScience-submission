# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Pamela Hathway p.hathway16@imperial.ac.uk
# ----------------------------------------------------------------------------- 
# Creates figure 8
# This script can be run independently.
# -----------------------------------------------------------------------------

from brian2 import *


def fig_8():
	print('%s Preparing Figure 8: EPSP shape' % time.strftime('%H:%M'))
	device.reinit()
	device.activate()

	plt.rcParams['savefig.dpi'] = 75
	plt.rcParams['figure.autolayout'] = False
	plt.rcParams['figure.figsize'] = 12, 5
	plt.rcParams['axes.labelsize'] = 17
	plt.rcParams['axes.titlesize'] = 17
	plt.rcParams['font.size'] = 17
	plt.rcParams['lines.linewidth'] = 2.0
	plt.rcParams['lines.markersize'] = 3
	plt.rcParams['legend.fontsize'] = 16

	runduration = 0.08
	defaultclock.dt = 1e-4 * second

	# parameters neuron and STDP
	taum = 10 * ms
	taus = 2.5 * ms
	tausyn = 2.5 * ms
	T = 2.8
	u_rest = 0
	X = (taus/taum)**(taum/(taus - taum))
	A = -8.5
	deltax = 1
	deltaa = 1
	deltau = 1.2
	win = 1

	# equations
	eqs = '''du/dt = (A*a)/taus + (X*x-u)/taum : 1
				dx/dt = -x/tausyn : 1 
				da/dt = -a/taus : 1'''

	eqs_reset = '''x = 0
				u = 2*T
				a = deltaa'''

	stdp_model = ''' wi : 1'''

	stdp_on_pre = ''' x_post += deltax*wi '''

	stdp_on_pre_imm = ''' u_post += deltau*wi '''

	# input
	times = array([5, 20, 40, 43, 46, 60]) * ms
	indices = zeros(len(times))

	N0 = SpikeGeneratorGroup(1, indices, times)
	N2 = NeuronGroup(1, eqs, threshold='u > T', reset=eqs_reset, refractory=2 * ms, method='linear')
	N3 = NeuronGroup(1, eqs, threshold='u > T', reset=eqs_reset, refractory=2 * ms, method='linear')
	N2spm = SpikeMonitor(N2)
	N3spm = SpikeMonitor(N3)
	N2stm = StateMonitor(N2, ['u'], record=True, dt=0.1 * ms)
	N3stm = StateMonitor(N3, ['u'], record=True, dt=0.1 * ms)

	syn02 = Synapses(N0, N2, model=stdp_model, on_pre=stdp_on_pre, method='linear')
	syn02.connect(j='i')
	syn02.wi = win
	syn03 = Synapses(N0, N3, model=stdp_model, on_pre=stdp_on_pre_imm, method='linear')
	syn03.connect(j='i')
	syn03.wi = win

	net = Network(collect())
	net.add(N0, N2, N3, N3spm, N3stm, syn03, N2spm, syn02)
	net.run(runduration * second)

	fig8 = figure(figsize=(6, 5))
	title('Figure 8')
	plot(N2stm.t / ms, N2stm.u[0, :], 'g', linewidth=2, label='kernel')
	plot(N2spm.t / ms, 6 + N2spm.i, 'g.')
	plot(N3stm.t / ms, N3stm.u[0, :], 'k', linewidth=2, label='imm. increase')
	plot(N3spm.t / ms, 6 + N3spm.i, 'k.')
	legend(loc=2)

	axhline(y=T, color='r', linestyle='--', linewidth=2, label='threshold')
	axhline(y=0, color='k', linestyle=':', linewidth=2, label='resting potential')
	xlabel('times [ms]')
	ylabel('Potential [a.u.]')
	xlim([0, 80])

	tight_layout()

	savefig('../article/figures/figure_8_created_%s.pdf' % datetime.datetime.now().strftime('%y%m%d'))

	return fig8


if __name__=='__main__':
	fig8 = fig_8()
	show()
