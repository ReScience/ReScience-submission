from brian2 import *
import os.path
import pickle
import os
import re
import glob

# ####### set options for plotting
plt.rcParams['savefig.dpi'] = 75
plt.rcParams['figure.autolayout'] = False
plt.rcParams['figure.figsize'] = 12, 5
plt.rcParams['axes.labelsize'] = 17
plt.rcParams['axes.titlesize'] = 17
plt.rcParams['font.size'] = 17
plt.rcParams['lines.linewidth'] = 2.0
plt.rcParams['lines.markersize'] = 8
plt.rcParams['legend.fontsize'] = 16

### from seed runs get results and save as dictionary

# para = load('/home/ph416/Documents/pf/para.npy')
# direc= '/home/ph416/Documents/pf/data/'

para = load('/Users/ph/Documents/pf/para.npy')
datadir = '/Users/ph/Documents/pf/data/'
savedir = '/Users/ph/Documents/pf/data/'
figdir = '/Users/ph/Documents/pf/figures/'


loopsta = 0
loopsto = 22
runstr = 'saved'

# loop over all conditions (22)
for run_idx in range(loopsta, loopsto):

	result = {}

	sta = run_idx * 100
	sto = sta + 99
	number_runs = sto - sta + 1

	# check if all files available
	avail_res = []
	missing = []
	for seed_idx in range(sta, sto + 1):
		file_path = ('%sresult_seedrunid=%s_%s.pickle' % (datadir, seed_idx, runstr))
		if os.path.exists(file_path):
			avail_res.append(seed_idx)
		else:
			missing.append(seed_idx)
	if len(missing) > 0:
		print('missing files!', missing)
	# break

	hits_array = empty(len(avail_res))
	fa_array = empty(len(avail_res))
	avglat_array = empty(len(avail_res))
	find_t_array = empty(len(avail_res))
	find_spike_array = empty(len(avail_res))
	N2_spikes_array = empty(len(avail_res))
	A_array = empty(len(avail_res))
	tres_array = empty(len(avail_res))

	no_result = []
	count = 0
	for seed_idx in avail_res:

		with open('%sresult_seedrunid=%s_%s.pickle' % (datadir, seed_idx, runstr), 'rb') as handle:
			run = pickle.load(handle)

		check_params = para[run_idx, :] - array([4, run['win'], run['jitter'], run['n_pat'],
												run['freq_pat'], run['del'], run['threshold']])

		if sum(check_params) != 0:
			print('parameters are not what they should be!')
			print('parameters should have been used')
			print('w', para[run_idx, 1], 'j', para[run_idx, 2], 'n', para[run_idx, 3], 'f', para[run_idx, 4],
				  'd', para[run_idx, 5], 'T', para[run_idx, 6])
			print('parameters that actually have been used')
			print('w', run['win'], 'j', run['jitter'], 'n', run['n_pat'], 'f', run['freq_pat'],
				  'd', run['del'], 'T', run['threshold'])

		if 'hits' in run.keys():
			hits_array[count] = run['hits']
			fa_array[count] = run['fa']
			avglat_array[count] = run['avg_lat']
			find_t_array[count] = run['find_t']
			find_spike_array[count] = run['find_spike']
			N2_spikes_array[count] = run['N2_spikes']
		else:
			# print('runid %s does not have results, writing -99' % seed_idx)
			hits_array[count] = -99
			fa_array[count] = -99
			avglat_array[count] = -99
			find_t_array[count] = -99
			find_spike_array[count] = -99
			N2_spikes_array[count] = -99
			no_result.append(seed_idx)

		A_array[count] = run['A']
		tres_array[count] = run['t_res']
		count += 1


	print('No results for: %s' % no_result)
	success_array_excllat = (hits_array >= 0.98) * (fa_array == 0)
	success_array = (hits_array >= 0.98) * (fa_array == 0) * (avglat_array < 10) * (avglat_array > 0)
	if sum(success_array) != sum(success_array_excllat):
		print('high latency for one run: incl lat ', sum(success_array), 'excl lat', sum(success_array_excllat))
	success = 100 * sum(success_array) / len(avail_res)
	print(run_idx, success)

	# result['name'] = 'result_%s_e%s_s%s-%s' % (runstr, int(-log10(run['t_res'])), avail_res[0], avail_res[-1])
	result['name'] = 'result_%s_s%s-%s' % (runstr, avail_res[0], avail_res[-1])
	result['runduration'] = run['runduration']  # runduration
	result['threshold'] = run['threshold']
	result['win'] = run['win']
	result['jitter'] = run['jitter']
	result['n_pat'] = run['n_pat']
	result['freq_pat'] = run['freq_pat']
	result['del'] = run['del']
	result['date'] = run['date']
	result['t_res'] = run['t_res']
	result['A'] = run['A']
	# result['refract'] = run['refract']
	result['K2'] = run['K2']

	result['N2_spikes'] = N2_spikes_array
	result['n_runs'] = number_runs
	result['n_runs_included'] = len(avail_res)
	result['seeds'] = array(avail_res)
	result['hits'] = hits_array
	result['fa'] = fa_array
	result['avg_lat'] = avglat_array
	result['find_t'] = find_t_array
	result['find_spike'] = find_spike_array
	result['success_overall'] = success
	result['success'] = success_array
	result['no_result'] = no_result
	result['missing'] = missing


	with open('%scluster_%s_e4.pickle' % (savedir, result['name'], result['date']), 'wb') as handle:
		pickle.dump(result, handle, protocol=pickle.HIGHEST_PROTOCOL)


# os.chdir(direc)
# files = glob.glob('dict_result_newsave_e6_s%s-%s_*.pickle' % (run_idx*100, run_idx*100 + 99))

# ### print all results
# run_ids = concatenate((arange(0, 16, 1), arange(19, 22, 1)))

# ### load all results into one array with 22 entries

runs_to_plot = [runstr]
names_to_plot = ['rep 10$^{-4}$']
styles = ['-']
colors = ['k']

win_ = [0.275, 0.325, 0.375, 0.425, 0.475]
jit_ = [0, 1, 2, 3, 4, 5, 6]
npat_ = [0.2, 0.3, 0.4, 0.5, 0.6]
freq_ = [0.05, 0.1, 0.15, 0.25, 0.5]
del_ = [0.0, 0.1, 0.2, 0.3]


figure(figsize=(12, 4))
plt.rcParams['legend.fontsize'] = 13


for idx, run_name in enumerate(runs_to_plot):

	run_ids = arange(0, 22, 1)
	success_rates = empty(22)
	for run_idx in run_ids:

		for filename in os.listdir(savedir):
			if re.match('cluster_result_%s_s%s-%s_\d{6}.pickle' % (run_name, run_idx * 100, run_idx * 100 + 99), filename):
				# print(filename)
				result_file = filename

		with open('%s%s' % (savedir, result_file), 'rb') as handle:
			run = pickle.load(handle)
		# print('run', run_idx, 'n_runs included', run['n_runs_included'], 'result', run['success_overall'])

		success_rates[run_idx] = sum(run['success'])

	win_res = concatenate((success_rates[1:5], array([success_rates[0]])))
	jit_res = concatenate((array([success_rates[5]]), array([success_rates[0]]), success_rates[6:11]))
	npat_res = concatenate((success_rates[11:14], array([success_rates[0]]), array([success_rates[14]])))
	freq_res = concatenate((success_rates[15:18], array([success_rates[0]]), array([success_rates[18]])))
	del_res = concatenate((array([success_rates[0]]), success_rates[19:22]))

	subplot(1, 5, 1)
	plot(freq_, freq_res, marker='o', color=colors[idx], linestyle=styles[idx])

	subplot(1, 5, 2)
	plot(jit_, jit_res, marker='o', color=colors[idx], linestyle=styles[idx])

	subplot(1, 5, 3)
	plot(npat_, npat_res, marker='o', color=colors[idx], linestyle=styles[idx])

	subplot(1, 5, 4)
	plot(win_, win_res, marker='o', color=colors[idx], linestyle=styles[idx], label=names_to_plot[idx])

	subplot(1, 5, 5)
	plot(del_, del_res, marker='o', color=colors[idx], linestyle=styles[idx])



subplot(1, 5, 1)
ylabel('% of success')
xlim([0, 0.55]), ylim([-5, 105])
xlabel('Pattern freq')
xticks(array([0.1, 0.3, 0.5]))
yticks(array([0, 25, 50, 75, 100]))
title('A', weight='bold')

subplot(1, 5, 2)
xlim([-0.5, 6.5]), ylim([-5, 105])
xlabel('Jitter SD [ms]')
yticks(array([0, 50, 100]))
xticks(array([0, 2, 4, 6]))
yticks(array([]))
title('B', weight='bold')

subplot(1, 5, 3)
xlim([0.175, 0.625]), ylim([-5, 105])
xlabel('Prop. pat. neur.')
yticks(array([0, 50, 100]))
xticks(array([0.2, 0.4, 0.6]))
yticks(array([]))
title('C', weight='bold')

subplot(1, 5, 4)
plot(array([-100, -200]), array([-5, -10]), marker='x', color='g', linestyle='--', label='orig paper')
xlim([0.25, 0.5]), ylim([-5, 105])
xlabel('Initial weight')
yticks(array([0, 50, 100]))
xticks(array([0.3, 0.4, 0.5]))
yticks(array([]))
title('D', weight='bold')
legend(loc=3)

subplot(1, 5, 5)
xlim([-0.025, 0.325]), ylim([0-5, 105])
xlabel('Spike deletion')
yticks(array([0, 50, 100]))
xticks(array([0, 0.1, 0.2, 0.3]))
yticks(array([]))
title('E', weight='bold')

tight_layout()

savefig('%sfig5_from_cluster.pdf' % figdir)


###############################
# ### rename files

# runstr_load = 'newsave_rerune4'
# runstr_save = 'newsave'
# datadir_load = '/Users/ph/Documents/pf/data/'
# for seed_idx in arange(1500, 1800, 1):
# 	with open('%sresult_seedrunid=%s_%s.pickle' % (datadir_load, seed_idx, runstr_load), 'rb') as handle:
# 		run_to_rename = pickle.load(handle)
#
# 	with open('%sresult_seedrunid=%s_%s.pickle' % (datadir_load, seed_idx, runstr_save), 'wb') as handle:
# 		pickle.dump(run_to_rename, handle, protocol=pickle.HIGHEST_PROTOCOL)
