# Implementation of the AuGMEnT network according to
# Rombouts JO, Bohte SM, Roelfsema PR (2015)
# How Attention Can Create Synaptic Tags for the Learning of Working Memories in Sequential Tasks.
# PLoS Comput Biol 11: e1004060.

import numpy as np
import random
from functools import reduce
from itertools import product
from collections import deque
from sys import exit

def sigmoid(x, theta):
    # from equation (5)
    return 1. / (1. + np.exp(theta - x))

def derivative_sigmoid(x):
    # from equation (6)
    return x * (1. - x)

class Population(object):
    """
    Base class for a population of neurons/units
    """
    
    def __init__(self, n, add_bias=False):
        if add_bias:
            self.biased_rates = np.zeros((n + 1,))
            self.biased_rates[0] = 1.
            self.rates = self.biased_rates[1:]
        else:
            self.rates = np.zeros((n,))
        
        self.preceding_connections = []
    
    def get_synaptic_inputs(self):
        return reduce(np.add, [connection.get_synaptic_outputs() for connection in self.preceding_connections])

class InstantaneousPopulation(Population):
    """
    Population of instantaneous sensory neurons/units (noted xi)
    """
    
    def __init__(self, n):
        super(InstantaneousPopulation, self).__init__(n, add_bias=True)
    
    def compute_rates(self, input):
        self.biased_rates[1:] = input

class TransientPopulation(Population):
    """
    Population of transient sensory neurons/units (noted x+ and x-)
    """
    
    def __init__(self, n):
        super(TransientPopulation, self).__init__(n)
        self.previous_input = np.zeros((n//2,))
    
    def compute_rates(self, input):
        # from equation (1)
        self.rates[:len(self.rates)//2] = (input - self.previous_input).clip(0.)
        # from equation (2)
        self.rates[len(self.rates)//2:] = (self.previous_input - input).clip(0.)
        self.previous_input = input

class AssociationPopulation(Population):
    """
    Base class for a population of association neurons/units (noted y)
    """
    
    def __init__(self, n, theta, add_bias=False):
        super(AssociationPopulation, self).__init__(n, add_bias)
        self.theta = theta
        self.preceding_feedback_connections = []
    
    def get_synaptic_feedbacks(self):
        return reduce(np.add, [connection.get_synaptic_feedbacks() for connection in self.preceding_feedback_connections])

class RegularPopulation(AssociationPopulation):
    """
    Population of regular association neurons/units (noted yR)
    """
    
    def __init__(self, n, theta):
        super(RegularPopulation, self).__init__(n, theta, add_bias=True)
    
    def compute_rates(self):
        # from equations (3) and (4)
        self.biased_rates[1:] = sigmoid(self.get_synaptic_inputs(), self.theta)

class MemoryPopulation(AssociationPopulation):
    """
    Population of memory association neurons/units (noted yM)
    """
    
    def __init__(self, n, theta):
        super(MemoryPopulation, self).__init__(n, theta)
        self.traces = np.zeros((n,))
    
    def compute_rates(self):
        # from equations (7) and (8)
        self.traces += self.get_synaptic_inputs()
        self.rates = sigmoid(self.traces, self.theta)
    
    def reset(self):
        self.traces = np.zeros(self.traces.shape)

class QValuePopulation(Population):
    """
    Population of Q-Value neurons/units (noted q)
    """
    
    def __init__(self, n):
        super(QValuePopulation, self).__init__(n)
        self.z = np.zeros((n,))
    
    def compute_rates(self):
        # from equation (11)
        self.rates = self.get_synaptic_inputs()

class Synapses(object):
    """
    Base class for synapse connections between two population of neurons/units
    """
    
    def __init__(self, presynaptic_neurons, postsynaptic_neurons, weights_range, beta, lambda_, gamma, add_bias=False):
        self.presynaptic_neurons = presynaptic_neurons
        self.postsynaptic_neurons = postsynaptic_neurons
        self.postsynaptic_neurons.preceding_connections.append(self)
        
        if add_bias:
            self.weights = np.random.uniform(*weights_range, size=(len(presynaptic_neurons.biased_rates), len(postsynaptic_neurons.rates)))
            self.tags = np.zeros((len(presynaptic_neurons.biased_rates), len(postsynaptic_neurons.rates)))
        else:
            self.weights = np.random.uniform(*weights_range, size=(len(presynaptic_neurons.rates), len(postsynaptic_neurons.rates)))
            self.tags = np.zeros((len(presynaptic_neurons.rates), len(postsynaptic_neurons.rates)))
        
        self.alpha = 1. - lambda_ * gamma
        self.beta = beta
    
    def get_synaptic_outputs(self):
        return np.dot(self.presynaptic_neurons.rates, self.weights)
    
    def update_weights(self, delta):
        # from equation (18)
        self.weights += self.beta * delta * self.tags
    
    def update_tags(self):
        pass
    
    def update(self, delta):
        self.update_weights(delta)
        self.update_tags()
    
    def reset(self):
        self.tags = np.zeros(self.tags.shape)
    
class InstantaneousRegularSynapses(Synapses):
    """
    Synapses between instantaneous sensory neurons and regular association neurons (noted vR)
    """
    
    def __init__(self, presynaptic_neurons, postsynaptic_neurons, weights_range, beta, lambda_, gamma):
        super(InstantaneousRegularSynapses, self).__init__(presynaptic_neurons, postsynaptic_neurons, weights_range, beta, lambda_, gamma, add_bias=True)
    
    def get_synaptic_outputs(self):
        return np.dot(self.presynaptic_neurons.biased_rates, self.weights)
    
    def update_tags(self):
        # from equation (14)
        self.tags += -self.alpha * self.tags
        self.tags += self.postsynaptic_neurons.get_synaptic_feedbacks() * derivative_sigmoid(self.postsynaptic_neurons.rates) * self.presynaptic_neurons.biased_rates[:,np.newaxis]
        # the numpy syntax "[:,np.newaxis]" creates a column view of a row vector
        # element-wise operations between a vector and a matrix are repeated on all rows or columns

class TransientMemorySynapses(Synapses):
    """
    Synapses between transient sensory neurons and memory association neurons (noted vM)
    """
    
    def __init__(self, presynaptic_neurons, postsynaptic_neurons, weights_range, beta, lambda_, gamma):
        super(TransientMemorySynapses, self).__init__(presynaptic_neurons, postsynaptic_neurons, weights_range, beta, lambda_, gamma)
        self.traces = np.zeros((len(presynaptic_neurons.rates), len(postsynaptic_neurons.rates)))
    
    def update_tags(self):
        # from equation (15)
        self.traces += self.presynaptic_neurons.rates[:,np.newaxis]
        # from equation (16)
        self.tags += -self.alpha * self.tags
        self.tags += self.postsynaptic_neurons.get_synaptic_feedbacks() * derivative_sigmoid(self.postsynaptic_neurons.rates) * self.traces
        
    def reset(self):
        super(TransientMemorySynapses, self).reset()
        self.traces = np.zeros(self.traces.shape)

class RegularQValueSynapses(Synapses):
    """
    Synapses between regular association neurons and Q-Value neurons (noted wR)
    """
    
    def __init__(self, presynaptic_neurons, postsynaptic_neurons, weights_range, beta, lambda_, gamma):
        super(RegularQValueSynapses, self).__init__(presynaptic_neurons, postsynaptic_neurons, weights_range, beta, lambda_, gamma, add_bias=True)
        presynaptic_neurons.preceding_feedback_connections.append(self)
    
    def update_tags(self):
        # from equation (13)
        self.tags += -self.alpha * self.tags
        self.tags += self.presynaptic_neurons.biased_rates[:,np.newaxis] * self.postsynaptic_neurons.z
    
    def get_synaptic_outputs(self):
        return np.dot(self.presynaptic_neurons.biased_rates, self.weights)
    
    def get_synaptic_feedbacks(self):
        return self.weights[1:,np.argmax(self.postsynaptic_neurons.z)]

class MemoryQValueSynapses(Synapses):
    """
    Synapses between memory associationn neurons and Q-Value neurons (noted wM)
    """
    
    def __init__(self, presynaptic_neurons, postsynaptic_neurons, weights_range, beta, lambda_, gamma):
        super(MemoryQValueSynapses, self).__init__(presynaptic_neurons, postsynaptic_neurons, weights_range, beta, lambda_, gamma)
        self.presynaptic_neurons.preceding_feedback_connections.append(self)
    
    def update_tags(self):
        # from equation (13)
        self.tags += -self.alpha * self.tags
        self.tags += self.presynaptic_neurons.rates[:,np.newaxis] * self.postsynaptic_neurons.z
    
    def get_synaptic_feedbacks(self):
        return self.weights[:,np.argmax(self.postsynaptic_neurons.z)]

class Network(object):
    """
    AuGMEnT (Attention-Gated MEmory Tagging) neural network
    """
    
    def __init__(self,
                 n_sensory_units,
                 n_qvalue_units,
                 n_regular_units = 3,
                 n_memory_units  = 4,
                 beta            = 0.15,
                 lambda_         = 0.2,
                 gamma           = 0.9,
                 epsilon         = 0.025,
                 theta           = 2.5,
                 weights_range   = (-0.25, 0.25)):
        """
        Constructor
        
        n_sensory_units -- number of sensory input neurons
        n_regular_units -- number of association regular neurons
        n_memory_units  -- number of association memory neurons
        n_qvalue_units  -- number of Q-Value output neurons
        beta            -- learning rate
        lambda_         -- tag decay rate
        gamma           -- discount factor
        epsilon         -- exploration rate
        theta           -- sigmoid function parameter
        weights_range   -- initial synapses' weights range
        """
        self.gamma = gamma
        self.epsilon = epsilon
        self.delta = 0.
        
        self.predicted_value = 0.
        self.previous_predicted_value = 0.
        
        self.learning = True
        self.exploration = True
        
        self.instantaneous_units = InstantaneousPopulation(n_sensory_units)
        self.transient_units = TransientPopulation(n_sensory_units * 2)
        self.regular_units = RegularPopulation(n_regular_units, theta)
        self.memory_units = MemoryPopulation(n_memory_units, theta)
        self.qvalue_units = QValuePopulation(n_qvalue_units)
        
        self.instantaneous_regular_synapses = InstantaneousRegularSynapses(self.instantaneous_units, self.regular_units, weights_range, beta, lambda_, gamma)
        self.transient_memory_synapses = TransientMemorySynapses(self.transient_units, self.memory_units, weights_range, beta, lambda_, gamma)
        self.regular_qvalue_synapses = RegularQValueSynapses(self.regular_units, self.qvalue_units, weights_range, beta, lambda_, gamma)
        self.memory_qvalue_synapses = MemoryQValueSynapses(self.memory_units, self.qvalue_units, weights_range, beta, lambda_, gamma)
    
    def feedforward(self, input):
        self.instantaneous_units.compute_rates(input)
        self.transient_units.compute_rates(input)
        self.regular_units.compute_rates()
        self.memory_units.compute_rates()
        self.qvalue_units.compute_rates()
        self.select_action()
    
    def feedback(self, reward, end):
        self.update_delta(reward, end)
        self.memory_qvalue_synapses.update(self.delta)
        self.regular_qvalue_synapses.update(self.delta)
        self.instantaneous_regular_synapses.update(self.delta)
        self.transient_memory_synapses.update(self.delta)
    
    def update_delta(self, reward, end):
        if end:
            self.predicted_value = 0.
        # from equation (17)
        self.delta = reward + self.gamma * self.predicted_value - self.previous_predicted_value
    
    def select_action(self):
        if (random.random() < self.epsilon) and self.exploration:
            selected = Network.random_from_boltzmann_distribution(self.qvalue_units.rates)
        else:
            max_value = self.qvalue_units.rates.max()
            selected = np.where(self.qvalue_units.rates == max_value)[0]
            if len(selected) != 1:
                selected = np.random.choice(selected)
        
        self.previous_predicted_value = self.predicted_value
        self.predicted_value = self.qvalue_units.rates[selected]
        self.qvalue_units.z = np.zeros(self.qvalue_units.z.shape)
        self.qvalue_units.z[selected] = 1.
    
    @staticmethod
    def random_from_boltzmann_distribution(a):
        # from equation (12)
        max_value = a.max()
        # for numerical stability
        exp_values = np.exp(a - max_value)
        exp_sum = np.sum(exp_values)
        
        i = 0
        tmp = exp_values[0]
        rand = random.random()
        
        while tmp / exp_sum < rand:
            i += 1
            tmp += exp_values[i]
        
        return i
    
    def step(self, input, reward, end):
        """
        Execute a step in the network
        Return the output of the network
        
        input  -- new sensory input
        reward -- reward given for the previous action selected
        end    -- signal the end of a trial
        """
        self.feedforward(input)
        if self.learning:
            self.feedback(reward, end)
        if end:
            self.reset()
        
        return self.qvalue_units.z
    
    def reset(self):
        """
        Reset the network for a new trial
        """
        self.memory_units.reset()
        self.instantaneous_regular_synapses.reset()
        self.transient_memory_synapses.reset()
        self.regular_qvalue_synapses.reset()
        self.memory_qvalue_synapses.reset()

class Task(object):
    """
    Base class for a task executable on an AuGMEnT network
    """
    
    def __init__(self):
        self.states = {"end": None}
        self.verbose = False
    
    def reset(self):
        self.t = 0
        self.current_state_t = 0
        self.next_reward = 0.
        self.success = False
        self.state = None
    
    def get_new_input(self):
        pass
    
    def step(self, output):
        self.states[self.state](output)
        return np.array(self.get_new_input())
    
    def run(self, network):
        """
        Execute a trial
        Return the success of the trial
        
        network -- AuGMEnT network on which to execute the task
        """
        self.reset()
        output = []
        
        while self.state != "end":
            input = self.step(output)
            output = network.step(input, self.next_reward, self.state == "end")
            
            self.current_state_t += 1
            self.t += 1
            self.next_reward = 0.
        
        return self.success
    
    def change_state(self, state):
        self.state = state
        self.current_state_t = 0

class FixationTask(Task):
    """
    Base class for saccade and probabilistic tasks
    """
    
    output_size = 3
    
    def __init__(self):
        super(FixationTask, self).__init__()
        self.fixation_reward = 0.2
        self.end_reward = 1.5
        self.max_t_init = 1
        self.max_t_wait = 10
        self.max_t_fixate = 1
        self.max_t_cue = 1
        self.max_t_delay = 2
        self.max_t_go = 8
        self.target = None
        self.states["init"] = self.state_init
        self.states["wait"] = self.state_wait
        self.states["fixate"] = self.state_fixate
        self.states["cue"] = self.state_cue
        self.states["delay"] = self.state_delay
        self.states["go"] = self.state_go
    
    def reset(self):
        super(FixationTask, self).reset()
        self.rewarded_target = None
        self.expected_target = None
        self.state = "init"
    
    @staticmethod
    def get_action(output):
        return next(action for i, action in enumerate(("left", "fixate", "right")) if output[i] == 1)
    
    def state_init(self, output):
        if self.current_state_t >= self.max_t_init:
            self.change_state("wait")
    
    def state_wait(self, output):
        if self.current_state_t >= self.max_t_wait:
            self.change_state("end")
        elif FixationTask.get_action(output) == "fixate":
            self.change_state("fixate")
    
    def state_fixate(self, output):
        if FixationTask.get_action(output) != "fixate":
            self.change_state("end")
        elif self.current_state_t >= self.max_t_fixate:
            self.next_reward = self.fixation_reward
            self.change_state("cue")
    
    def state_cue(self, output):
        if FixationTask.get_action(output) != "fixate":
            self.change_state("end")
        elif self.current_state_t >= self.max_t_cue:
            self.change_state("delay")
    
    def state_delay(self, output):
        if FixationTask.get_action(output) != "fixate":
            self.change_state("end")
        elif self.current_state_t >= self.max_t_delay:
            self.change_state("go")
    
    def state_go(self, output):
        if self.current_state_t >= self.max_t_go:
            self.change_state("end")
        elif FixationTask.get_action(output) != "fixate":
            self.next_reward = self.end_reward if FixationTask.get_action(output) == self.rewarded_target else 0.
            self.success = FixationTask.get_action(output) in self.expected_target
            self.change_state("end")

class SaccadeTask(FixationTask):
    """
    Saccade/Antisaccade task as described in paper
    """
    
    description = "saccade/antisaccade task"
    input_size = 4
    
    def __init__(self):
        super(SaccadeTask, self).__init__()
        # to make the task pro-saccades only : self.fixations = ("pro",)
        self.fixations = ("pro", "anti")
        self.cues = ("left", "right")
    
    def reset(self):
        super(SaccadeTask, self).reset()
        
        self.fixation = random.choice(self.fixations)
        self.cue = random.choice(self.cues)
        
        self.rewarded_target = self.cue if self.fixation == "pro" else next(cue for cue in ("left", "right") if cue != self.cue)
        self.expected_target = self.rewarded_target
    
    def get_new_input(self):
        if self.state in ("init", "go", "end"):
            return [0, 0, 0, 0]
        elif self.state in ("wait", "fixate", "delay"):
            return [1 if action in (self.fixation) else 0 for action in ("left", "pro", "anti", "right")]
        elif self.state == "cue":
            return [1 if action in (self.fixation, self.cue) else 0 for action in ("left", "pro", "anti", "right")]
    
    def evaluate(self, network):
        """
        Evaluate the network on all trial types
        Return True if the network successfuly learned the task
        
        network -- network to evaluate
        """
        learning = network.learning
        exploration = network.exploration
        fixations = self.fixations
        cues = self.cues
        
        network.learning = False
        network.exploration = False
        
        success = True
        for fixation, cue in product(fixations, cues):
            self.fixations = (fixation,)
            self.cues = (cue,)
            if not self.run(network):
                success = False
                break;
        
        network.learning = learning
        network.exploration = exploration
        self.fixations = fixations
        self.cues = cues
        
        return success
    
    def train(self, network):
        """
        Run trials until convergence or trial limit
        Return the convergence time or None if the network did not converge
        
        network -- network to train
        """
        n_max_trials = 25000
        success_goal = 0.9
        sample_size = 50
        samples = {}
        success_rates = {}
        for fixation, cue in product(self.fixations, self.cues):
            type = fixation + cue
            samples[type] = deque([0] * sample_size)
            success_rates[type] = 0.
        i = 0
        
        while (i < n_max_trials) and any(success_rates[key] < success_goal for key in success_rates.keys()):
            success = self.run(network)
            type = self.fixation + self.cue
            
            samples[type].popleft()
            samples[type].append(1 if success else 0)
            success_rates[type] = float(sum(samples[type])) / sample_size
            
            i += 1
        
        return i if self.evaluate(network) else None

class SaccadeNoShapingTask(SaccadeTask):

    description = "saccade/antisaccade task without shaping strategy"
    
    def __init__(self):
        super(SaccadeNoShapingTask, self).__init__()
        self.fixation_reward = 0.

class ProbabilisticTask(FixationTask):
    """
    Proabilistic decision making task as defined in paper
    """
    
    description = "probabilistic decision making task"
    input_size = 45
    infinity = 999.
    max_weight = 0.9 * 4
    shape_weights = [infinity, -infinity, 0.3, -0.3, 0.9, -0.9, 0.7, -0.7, 0.5, -0.5]
    difficulties = {1: {"n_symbols": 2, "seq_length": 1, "sample_size": 1000},
                    2: {"n_symbols": 4, "seq_length": 1, "sample_size": 1500},
                    3: {"n_symbols": 6, "seq_length": 1, "sample_size": 2000},
                    4: {"n_symbols": 8, "seq_length": 1, "sample_size": 2500},
                    5: {"n_symbols": 10, "seq_length": 1, "sample_size": 3000},
                    6: {"n_symbols": 10, "seq_length": 2, "sample_size": 10000},
                    7: {"n_symbols": 10, "seq_length": 3, "sample_size": 10000},
                    8: {"n_symbols": 10, "seq_length": 4, "sample_size": 20000}}
    
    def __init__(self):
        super(ProbabilisticTask, self).__init__()
        self.difficulty = 1
        self.targets_types = ("red_green", "green_red")
    
    def reset(self):
        super(ProbabilisticTask, self).reset()
        
        self.targets_type = random.choice(self.targets_types)
        self.shapes = np.random.randint(self.difficulties[self.difficulty]["n_symbols"], size=self.difficulties[self.difficulty]["seq_length"])
        self.locations = random.sample(range(4), len(self.shapes))
        self.max_t_cue = len(self.shapes)
        
        probability = self.probability()
        if probability == 0.5:
            self.expected_target = ("left", "right")
        if self.targets_type == "red_green":
            self.expected_target = "left" if probability > 0.5 else "right"
            self.rewarded_target = "left" if random.random() < probability else "right"
        elif self.targets_type == "green_red":
            self.expected_target = "right" if probability > 0.5 else "left"
            self.rewarded_target = "right" if random.random() < probability else "left"
    
    def get_new_input(self):
        input = [0] * self.input_size
        
        if self.state in ("init", "end"):
            return input
        elif self.state in ("wait", "fixate", "cue", "delay"):
            input[0] = 1
            if self.targets_type == "green_red":
                input[1] = 1
                input[2] = 1
            elif self.targets_type == "red_green":
                input[3] = 1
                input[4] = 1
        if self.state == "cue":
            for i in range(self.current_state_t + 1):
                input[5 + self.locations[i] * len(self.shape_weights) + self.shapes[i]] = 1
        
        return input
    
    def probability(self):
        """
        Return the probability of rewarding the red target
        """
        weights_sum = sum([self.shape_weights[shape] for shape in self.shapes])
        if weights_sum < -self.max_weight:
            return 0.
        elif weights_sum > self.max_weight:
            return 1.
        power = np.power(10., weights_sum)
        return power / (1. + power)
    
    def train(self, network):
        """
        Run trials until convergence or trial limit
        Return the convergence time
        
        network -- network to train
        """
        n_max_trials = 500000
        success_goal = 0.85
        i = 0
        
        for key in self.difficulties.keys():
            self.difficulty = key
            sample_size = self.difficulties[key]["sample_size"]
            samples = deque([0] * sample_size)
            success_rate = 0.
            
            while (i < n_max_trials) and (success_rate < success_goal):
                success = self.run(network)
                
                samples.popleft()
                samples.append(1 if success else 0)
                success_rate = float(sum(samples)) / sample_size
                
                i += 1
        
        return i if i < n_max_trials else None