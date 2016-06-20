import augment
import matplotlib.pyplot as plt
from itertools import product

class ObservableNetwork(augment.Network):
    def __init__(self, n_sensory_units, n_qvalue_units):
        super(ObservableNetwork, self).__init__(n_sensory_units, n_qvalue_units)
        self.recording = False
        self.record = {"left": [], "fixate": [], "right": []}
    
    def select_action(self):
        if self.recording:
            self.record["left"].append(self.qvalue_units.rates[0])
            self.record["fixate"].append(self.qvalue_units.rates[1])
            self.record["right"].append(self.qvalue_units.rates[2])
        super(ObservableNetwork, self).select_action()

if __name__ == "__main__":
    task = augment.SaccadeTask()
    network = ObservableNetwork(task.input_size, task.output_size)
    
    print("Training network...")
    
    convergence = None
    while convergence is None:
        convergence = task.train(network)
    
    print("Done.")
    
    print("Plotting activity...")
    
    network.beta = 0.
    network.epsilon = 0.
    network.recording = True
    
    fig, axes = plt.subplots(nrows=4, sharex=True, sharey=True)
    indices = {"proleft": 0,
               "proright": 1,
               "antileft": 2,
               "antiright": 3}
    titles = {"proleft": "Pro-Saccade/Left Cue",
              "proright": "Pro-Saccade/Right Cue",
              "antileft": "Anti-Saccade/Left Cue",
              "antiright": "Anti-Saccade/Right Cue"}
    
    for fixation, cue in product(task.fixations, task.cues):
        task.fixations = (fixation,)
        task.cues = (cue,)
        # ugly
        network.record = {"left": [], "fixate": [], "right": []}
        task.run(network)
        type = fixation + cue
        axes[indices[type]].set_title(titles[type], size=12)
        axes[indices[type]].plot(network.record["fixate"], "b")
        axes[indices[type]].plot(network.record["left"], "g")
        axes[indices[type]].plot(network.record["right"], "r")
    
    fig.suptitle("Q-value Units Activation Traces", size=12)
    fig.legend(axes[0].get_lines(), ("fixate", "left", "right"), loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.01))
    # plt.xticks(range(len(task.record["state"])), task.record["state"])
    
    print("Displaying graphs (Close the window to end the program).")
    
    plt.show()