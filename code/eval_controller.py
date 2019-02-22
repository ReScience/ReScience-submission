"""
Training a linear controller on latent + recurrent state
with CMAES.

This is a bit complex. num_workers slave threads are launched
to process a queue filled with parameters to be evaluated.
"""
import argparse
import sys
from os.path import join, exists
from os import mkdir, unlink, listdir, getpid
from time import sleep
import torch

from models import Controller
from tqdm import tqdm
import numpy as np
from utils.misc import RolloutGenerator, ASIZE, RSIZE, LSIZE

import logging


import pdb
################################################################################
#                           Thread routines                                    #
################################################################################
def slave_routine(p_queue, r_queue, e_queue, p_index, logdir):
    """ Thread routine.

    Threads interact with p_queue, the parameters queue, r_queue, the result
    queue and e_queue the end queue. They pull parameters from p_queue, execute
    the corresponding rollout, then place the result in r_queue.

    Each parameter has its own unique id. Parameters are pulled as tuples
    (s_id, params) and results are pushed as (s_id, result).  The same
    parameter can appear multiple times in p_queue, displaying the same id
    each time.

    As soon as e_queue is non empty, the thread terminate.

    When multiple gpus are involved, the assigned gpu is determined by the
    process index p_index (gpu = p_index % n_gpus).

    :args p_queue: queue containing couples (s_id, parameters) to evaluate
    :args r_queue: where to place results (s_id, results)
    :args e_queue: as soon as not empty, terminate
    :args p_index: the process index
    """
    # init routine
    
    # device = torch.device('cpu')
    # redirect streams
    tmp_dir = join(logdir, 'tmp', str(getpid()))
    sys.stdout = open(tmp_dir + '.out', 'a')
    sys.stderr = open(tmp_dir + '.err', 'a')


    handler = logging.FileHandler(tmp_dir + '.out')    
    logger = logging.getLogger('main_'+str(getpid()))
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)

    gpu = p_index % torch.cuda.device_count()
    device = torch.device('cuda:{}'.format(gpu) if torch.cuda.is_available() else 'cpu')
    try:
        with torch.no_grad():
            time_limit = 1000
            r_gen = RolloutGenerator(logdir, device, time_limit, logger)
            while e_queue.empty():
                if p_queue.empty():
                    sleep(.1)
                else:
                    s_id, params = p_queue.get()
                    X = r_gen.rollout(params)
                    r_queue.put((s_id, X))
    except Exception:
        logger.error(f"Fatal error in process {p_index}", exc_info=True)


################################################################################
#                           Evaluation                                         #
################################################################################
def evaluate(rollouts, p_queue, r_queue):
    """ Give current controller evaluation.

    Evaluation is minus the cumulated reward averaged over rollout runs.

    :args solutions: CMA set of solutions
    :args results: corresponding results
    :args rollouts: number of rollouts

    :returns: minus averaged cumulated reward
    """    
    results = []
    for s_id in range(rollouts):
        p_queue.put((s_id, None))
    for _ in tqdm(range(rollouts)):
        while r_queue.empty():
            sleep(.1)
        results.append(r_queue.get()[1])
    return np.mean(results), np.std(results)


def main(logdir):
    p_queue = mp.Queue()
    r_queue = mp.Queue()
    e_queue = mp.Queue()

    # multiprocessing variables
    tmp_dir = join(logdir, 'tmp')
    if not exists(tmp_dir):
        mkdir(tmp_dir)
    else:
        for fname in listdir(tmp_dir):
            unlink(join(tmp_dir, fname))

    ctrl_dir = join(logdir, 'ctrl')
    assert exists(ctrl_dir)

    
    p_list = []
    for p_index in range(args.max_workers):
        p_list.append(mp.Process(target=slave_routine, args=(p_queue, r_queue, e_queue, p_index, logdir)).start())

    controller = Controller(LSIZE, RSIZE, ASIZE)  # dummy instance

    ctrl_file = join(ctrl_dir, 'best.tar')
    assert exists(ctrl_file)
    state = torch.load(ctrl_file, map_location={'cuda:0': 'cpu'})
    controller.load_state_dict(state['state_dict'])
    best, std_best = evaluate(args.rollouts_eval, p_queue, r_queue)
    
    print(f"Evaluation: {best:.3f}+-{std_best:.3f}")
    for p in p_list:
        if p is not None:
            p.terminate()

if __name__ == '__main__':
    # parsing
    # torch.multiprocessing.set_start_method("spawn")
    mp = torch.multiprocessing.get_context('spawn')
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--logdir', type=str, help='Where everything is stored.')
    parser.add_argument('--display', action='store_true', help="Use progress bars if "
                        "specified.")
    parser.add_argument('--max-workers', type=int, help='Maximum number of workers.',
                        default=32)
    parser.add_argument('--rollouts-eval', type=int, help='Number of rollouts used in evaluation', default=100)
    args = parser.parse_args()


    main(args.logdir)
    
    