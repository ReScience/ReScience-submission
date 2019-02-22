import numpy as np
import argparse
from utils.history import History
from os.path import join
import matplotlib.pyplot as plt

plot_args = dict(
    standard=dict(
        label='Trained RNN and VAE',
        k='crimson'
    ),
    untrainedrnn=dict(
        label='Untrained RNN, trained VAE',
        k='navy'
    ),
    untrainedrnnvae=dict(
        label='Unrained RNN and VAE',
        k='forestgreen')
    )

def gen_plot(histories, out):

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    for exp_name, history in histories:

        h = history._hist['Return']

        x = np.array(list(h.keys()))
        maxi = - np.array([np.min(r_list) for r_list in h.values()])
        mean = - np.array([np.mean(r_list) for r_list in h.values()])
        # mini = - np.array([np.max(r_list) for r_list in h.values()])
        std = np.array([np.std(r_list) for r_list in h.values()])
        ax.plot(x, mean, c=plot_args[exp_name]['k'], label=plot_args[exp_name]["label"])
        # ax.plot(x, mini, c=plot_args[exp_name]['k'], linestyle=':', alpha=0.5)
        ax.plot(x, maxi, c=plot_args[exp_name]['k'], linestyle='--', alpha=0.8)
        ax.fill_between(x, mean-std, mean+std, facecolor=plot_args[exp_name]['k'], alpha=0.2)

    ax.set_xlim(0, 195)
    ax.set_xlabel("Number of epochs")
    ax.set_ylabel("Return")
    ax.set_ylim(-120, 920)
    ax.set_yticks([100*x for x in range(-1, 10)], minor=False)

    ax.legend()
    plt.savefig(out, format='png', dpi=500)
                

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--experiments-folder', type=str, help='folder in which the experiment directories are stored')
    parser.add_argument('--experiments-names', nargs='*', type=str, help='log directories of the experiments')
    parser.add_argument('--out', type=str, default='out.png', help='Filename of the figure')

    
    args = parser.parse_args()


    histories = [(expname, History(join(args.experiments_folder, expname, 'ctrl', 'history.pkl'))) \
                 for expname in args.experiments_names]

    gen_plot(histories, out=args.out)
    
