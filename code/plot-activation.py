import augment
import argparse
import matplotlib.pyplot as plt

def plot_activity(task_class, trial_types, setup_function, titles, colors):
    task = task_class()
    network = augment.Network(task.input_size, task.output_size)
    
    print("Training network for", task.description)
    
    convergence = None
    while convergence is None:
        convergence = task.train(network)
    
    network.beta = 0.
    network.epsilon = 0.
    step = network.step
    def decorated_step(input, reward, end):
        output = step(input, reward, end)
        record["left"].append(network.qvalue_units.rates[0])
        record["fixate"].append(network.qvalue_units.rates[1])
        record["right"].append(network.qvalue_units.rates[2])
        record["state"].append(task.state)
        return output
    network.step = decorated_step
    
    print("Plotting activity")
    
    figure, axes = plt.subplots(nrows=len(trial_types), sharex=False, sharey=True)
    
    for i, trial_type in enumerate(trial_types):
        record = {"left": [],
                  "fixate": [],
                  "right": [],
                  "state": []}
        setup_function(task, trial_type)
        task.run(network)
        axes[i].set_title(titles[i], size=12)
        axes[i].plot(record["fixate"], colors[0])
        axes[i].plot(record["left"], colors[1])
        axes[i].plot(record["right"], colors[2])
        axes[i].set_xticks(range(len(record["state"])))
        axes[i].set_xticklabels(record["state"])
    
    figure.legend(axes[0].get_lines(), ("fixate", "left", "right"), loc="lower center", ncol=3, prop={"size": 12})
    figure.subplots_adjust(bottom=0.12, top=0.95, hspace=0.55)
    
    print("Displaying graph")
    
    plt.show()

def setup_saccade(task, trial_type):
    task.fixations = (trial_type[0],)
    task.cues = (trial_type[1],)

def setup_probabilistic(task, trial_type):
    task.targets_types = (trial_type,)
    task.input_symbols = (8,)
    task.sequence_length = 4

if __name__ == "__main__":
    tasks = {"saccade": {"task_class": augment.SaccadeTask,
                         "trial_types": (("pro", "left"), ("pro", "right"), ("anti", "left"), ("anti", "right")),
                         "setup_function": setup_saccade,
                         "titles": ("Pro-Saccade/Left Cue", "Pro-Saccade/Right Cue", "Anti-Saccade/Left Cue", "Anti-Saccade/Right Cue"),
                         "colors": ("b", "g", "r")},
            "probabilistic": {"task_class": augment.ProbabilisticTask,
                              "trial_types": ("green_red", "red_green"),
                              "setup_function": setup_probabilistic,
                              "titles": ("Green/Red", "Red/Green"),
                              "colors": ("b", "m", "y")}}
    
    parser = argparse.ArgumentParser(description="Plot Q-value units' activation of a trained network")
    parser.add_argument("-t", "--task", default="saccade", type=str, choices=tasks.keys(), required=False, help="task to run")
    args = parser.parse_args()
    
    plot_activity(**tasks[args.task])