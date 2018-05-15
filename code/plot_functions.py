from brian2 import *


def make_figure_3(N2spm, latency):
	figure()
	subplot(121)
	title('latency vs time')
	plot(N2spm.t / second, latency * 1000, 'g.', markersize=3)
	xlabel('time [s]')
	yticks([0, 10, 20, 30, 40, 50])
	ylim([-2, 52])

	subplot(122)
	title('latency vs output neuron spikes')
	plot(range(len(latency)), latency * 1000, 'g.', markersize=3)
	ylabel('Postsynaptic spike latency [ms]')
	xlabel('# discharge')
	yticks([0, 10, 20, 30, 40, 50])
	ylim([-2, 52])

	tight_layout()
	show()


def make_figure_4(N1spm, syn12, timing_pattern):
	# get spikes from end of simulation
	figsta = searchsorted(N1spm.t /second, timing_pattern[-2] - 0.025)
	figsto = searchsorted(N1spm.t /second, timing_pattern[-2] + 0.075)
	greyness = [str((255 - syn12.wi[int(i)]) / 255.) for i in N1spm.i[figsta:figsto]]  # high w == white
	# greyness = [str(w[int(i)]/255.) for i in indices[figsta:figsto] ]

	ax = figure(figsize=(8, 6)).add_subplot(1, 1, 1, facecolor='k')
	ax.scatter(N1spm.t[figsta:figsto] / second, N1spm.i[figsta:figsto], c=greyness, s=3, lw=0, cmap='Greys')
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
	title('converged state')

	tight_layout()
	show()

