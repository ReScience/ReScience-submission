# -----------------------------------------------------------------------------
# Distributed under the GNU General Public License.
#
# Contributors: Pamela Hathway p.hathway16@imperial.ac.uk
# ----------------------------------------------------------------------------- 
# Functions to create figures 3, 4, 7AB
# This script cannot be run independently but is called by run_simulation.py
# -----------------------------------------------------------------------------

from brian2 import *


def fig_3(N2spm, latency):
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

	tight_layout()

	savefig('../article/figures/figure_3_created_%s.pdf' % datetime.datetime.now().strftime('%y%m%d'))

	return fig3


def fig_4(N1spm, syn12, timing_pattern):
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
	figsta = searchsorted(N1spm.t /second, timing_pattern[-2] - 0.025)
	figsto = searchsorted(N1spm.t /second, timing_pattern[-2] + 0.075)
	greyness = [str((255 - syn12.wi[int(i)]) / 255.) for i in N1spm.i[figsta:figsto]]  # high w == white

	ax = figure(figsize=(8, 6)).add_subplot(1, 1, 1)
	ax.set_facecolor("k")
	ax.scatter(N1spm.t[figsta:figsto] / second, N1spm.i[figsta:figsto], c=greyness, s=3, lw=0, cmap='Greys')
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

	tight_layout()

	savefig('../article/figures/figure_4_created_%s.pdf' % datetime.datetime.now().strftime('%y%m%d'))

	return ax


def fig_7AB(syn12stm, syn12nnstm, syn12atastm):
	fig = figure(figsize=(12, 5))
	suptitle('Figure 7 A and B')
	subplot(121)
	title('first second')

	plot(arange(0, len(syn12stm.wi[0,:]), 1), mean(syn12stm.wi, axis=0), 'k', label='RNN')
	plot(arange(0, len(syn12nnstm.wi[0,:]), 1), mean(syn12nnstm.wi, axis=0), 'k--', label='NN')
	plot(arange(0, len(syn12atastm.wi[0,:]), 1), mean(syn12atastm.wi, axis=0), 'k-.', label='ATA')

	ylabel('average weight per synapse')
	xlabel('time [s]')
	ylim([0.15, 0.55])
	xticks([0, 0.5, 1])
	yticks(array([0.2, 0.3, 0.4, 0.5]))
	xlim([0, 1])

	subplot(122)
	title('whole simulation')
	plot(arange(0, len(syn12stm.wi[0,:]), 1), mean(syn12stm.wi, axis=0), 'k', label='RNN')
	plot(arange(0, len(syn12nnstm.wi[0,:]), 1), mean(syn12nnstm.wi, axis=0), 'k--', label='NN')
	plot(arange(0, len(syn12atastm.wi[0,:]), 1), mean(syn12atastm.wi, axis=0), 'k-.', label='ATA')

	xlabel('time [s]')
	ylim([0.15, 0.55])
	xticks(array([1, 150, 300, 450]))
	yticks(array([0.2, 0.3, 0.4, 0.5]))
	xlim([1, 450])
	legend(loc=1)

	tight_layout()

	savefig('../article/figures/figure_7AB_created_%s.pdf' % datetime.datetime.now().strftime('%y%m%d'))

	return fig


if __name__=='__main__':
	print('These figures rely on output from run_simulation.py. '
		  'Please run that and these figures will be generated automatically.')


