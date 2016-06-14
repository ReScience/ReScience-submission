import augment
import argparse
import multiprocessing
import time
from statistics import median

def unsigned_int(string):
    value = int(string)
    if value <= 0:
        raise argparse.ArgumentTypeError("must be > 0")
    return value

def process(task_class):
    task = task_class()
    network = augment.Network(task.input_size, task.output_size)
    convergence_time = task.train(network)
    return convergence_time

if __name__ == "__main__":
    tasks = {"saccade": augment.SaccadeTask, "saccadenoshaping": augment.SaccadeNoShapingTask, "probabilistic": augment.ProbabilisticTask}
    
    parser = argparse.ArgumentParser(description="Train AuGMEnT networks and print median convergence time")
    parser.add_argument("-n", "--networks", default=10, type=unsigned_int, required=False, help="number of networks to train")
    parser.add_argument("-t", "--task", default="saccade", type=str, choices=tasks.keys(), required=False, help="task to run")
    parser.add_argument("-p", "--processors", default=multiprocessing.cpu_count(), type=unsigned_int, required=False, help="number of processors to use")
    args = parser.parse_args()
    
    pool = multiprocessing.Pool(min(args.processors, args.networks, multiprocessing.cpu_count()))
    results = []
    
    print("Training", args.networks, "network(s) for", tasks[args.task].description)
    
    start_time = time.time()
    
    for i in range(args.networks):
        pool.apply_async(process, (tasks[args.task],), callback=results.append)
    pool.close()
    
    while len(results) != args.networks:
        time.sleep(1.)
        print("\r{0:.1f}% completed in {1:.1f}s".format(100.0 * len(results) / args.networks, time.time() - start_time), end="")
    
    print("\n", end="")
    
    convergence_times = [result for result in results if result is not None]
    print("Success rate: {0:.2f}%".format(100.0 * len(convergence_times) / args.networks))
    print("Median convergence time:", median(convergence_times), "trials")