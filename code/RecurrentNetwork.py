# Main network class as presented in:
#
#   Laje, R. and Buonomano, D.V. (2013). Robust timing and motor patterns by taming chaos in recurrent neural networks. Nat Neurosci.
#
# Author: Julien Vitay (julien.vitay@informatik.tu-chemnitz.de)
# Licence: MIT
from __future__ import print_function
import numpy as np

class RecurrentNetwork(object):
    """
    Class implementing a recurrent network with read-out weights and RLS learning rules.

    **Parameters:**

    * Ni : Number of input neurons
    * N : Number of recurrent neurons
    * No : Number of read-out neurons
    * tau : Time constant of the neurons
    * g : Synaptic strength scaling
    * pc : Connection probability
    * Io : Noise variance
    * delta : Initial value of the P matrix
    * P_plastic : Percentage of neurons receiving plastic synapses
    """
    def __init__(self, Ni=2, N=800, No=1, tau=10.0, g=1.5, pc=0.1, Io=0.001, delta=1.0, P_plastic=0.6):
        # Copy the parameters
        self.Ni = Ni
        self.N = N
        self.No = No
        self.tau = tau
        self.g = g
        self.pc = pc
        self.Io = Io
        self.delta = delta
        self.P_plastic = P_plastic
        self.N_plastic = int(self.P_plastic*self.N) # Number of plastic cells = 480

        # Build the network
        self.build()

    def build(self):
        """
        Initializes the network including the weight matrices.
        """
        # Input
        self.I = np.zeros((self.Ni, 1))

        # Recurrent population
        self.x = np.random.uniform(-1.0, 1.0, (self.N, 1))
        self.r = np.tanh(self.x)

        # Read-out population
        self.z = np.zeros((self.No, 1))

        # Weights between the input and recurrent units
        self.W_in = np.random.randn(self.N, self.Ni)

        # Weights between the recurrent units
        self.W_rec = np.random.randn(self.N, self.N) * self.g/np.sqrt(self.pc*self.N)

        # The connection pattern is sparse with p=0.1
        connectivity_mask = np.random.binomial(1, self.pc, (self.N, self.N))
        connectivity_mask[np.diag_indices(self.N)] = 0
        self.W_rec *= connectivity_mask

        # Store the pre-synaptic neurons to each plastic neuron
        self.W_plastic = [list(np.nonzero(connectivity_mask[i, :])[0]) for i in range(self.N_plastic)]

        # Inverse correlation matrix of inputs for learning recurrent weights
        self.P = [1./self.delta*np.identity(len(self.W_plastic[i])) for i in range(self.N_plastic)]

        # Output weights
        self.W_out = np.random.randn(self.No, self.N) / np.sqrt(self.N)

        # Inverse correlation matrix of inputs for learning readout weights
        self.P_out = [1./self.delta*np.identity(self.N) for i in range(self.No)]

    def simulate(self, stimulus, noise=True, trajectory=np.array([]), learn_start=-1, learn_stop=-1, learn_readout=False, verbose=True):
        """
        Simulates the recurrent network for the given duration, with or without plasticity.

        * `stimulus`: np.array for the inputs. Determines the duration.
        * `noise`: if noise should be added to the recurrent units (default: True)
        * `trajectory`: during learning, defines which target function should be learned (default: no learning)
        * `learn_start`: time when learning should start.
        * `learn_stop`: time when learning should stop.
        * `learn_readout`: defines whether the recurrent (False) or readout (True) weights should be learned.
        * `verbose`: defines if the loss should be printed (default: True)
        """
        # Arrays for recording
        record_r = []
        record_z = []

        # Get the stimulus shape to know the duration
        nb, dummy, duration = stimulus.shape

        # Reset the recurrent population
        self.x = np.random.uniform(-1.0, 1.0, (self.N, 1))
        self.r = np.tanh(self.x)

        # Reset loss term
        self.loss = 0.0

        # Simulate for the desired duration
        for t in range(duration):

            # Update the neurons' firing rates
            self.update_neurons(stimulus[:, :, t], noise)

            # Recording
            record_r.append(self.r)
            record_z.append(self.z)

            # Learning
            if trajectory.size > 0 and t>=learn_start and t<learn_stop and t%2==0:
                if not learn_readout:
                    self.rls_recurrent(trajectory[t, :, :])
                else:
                    self.rls_readout(trajectory[:, :, t])

        # Print the loss at the end of the trial
        if trajectory.size > 0 and verbose:
            print('\tLoss:', self.loss/(learn_stop-learn_start)*2.)

        return np.array(record_r), np.array(record_z)

    def update_neurons(self, stimulus, noise):
        """
        Updates neural variables for a single simulation step.
        """
        # Inputs are set externally
        self.I = stimulus
        # Noise can be shut off
        I_noise = self.Io * np.random.randn(self.N, 1) if noise else 0.0
        # tau * dx/dt + x = I + sum(r) + I_noise
        self.x += (np.dot(self.W_in, self.I) + np.dot(self.W_rec, self.r) + I_noise - self.x)/self.tau
        # r = tanh(x)
        self.r = np.tanh(self.x)
        # z = sum(r)
        self.z = np.dot(self.W_out, self.r)

    def rls_recurrent(self, target):
        """
        Applies the RLS learning rule to the recurrent weights.
        """
        # Compute the error of the recurrent neurons
        error = self.r - target
        self.loss += np.mean(error**2)
        # Apply the FORCE learning rule to the recurrent weights
        for i in range(self.N_plastic): # for each plastic post neuron
            # Get the rates from the plastic synapses only
            r_plastic = self.r[self.W_plastic[i]]
            # Multiply with the inverse correlation matrix P*R
            PxR = np.dot(self.P[i], r_plastic)
            # Normalization term 1 + R'*P*R
            RxPxR = (1. + np.dot(r_plastic.T,  PxR))[0, 0]
            # Update the inverse correlation matrix P <- P - ((P*R)*(P*R)')/(1+R'*P*R)
            self.P[i] -= np.dot(PxR, PxR.T)/RxPxR
            # Learning rule W <- W - e * (P*R)/(1+R'*P*R)
            self.W_rec[i, self.W_plastic[i]] -= error[i, 0] * (PxR[:, 0]/RxPxR)

    def rls_readout(self, target):
        """
        Applies the RLS learning rule to the readout weights.
        """
        # Compute the error of the output neurons
        error = self.z - target
        self.loss += np.mean(error**2)
        # Apply the FORCE learning rule to the readout weights
        for i in range(self.No): # for each readout neuron
            # Multiply the rates with the inverse correlation matrix P*R
            PxR = np.dot(self.P_out[i], self.r)
            # Normalization term 1 + R'*P*R
            RxPxR = (1. + np.dot(self.r.T,  PxR))[0, 0]
            # Update the inverse correlation matrix P <- P - ((P*R)*(P*R)')/(1+R'*P*R)
            self.P_out[i] -= np.dot(PxR, PxR.T)/RxPxR
            # Learning rule W <- W - e * (P*R)/(1+R'*P*R)
            self.W_out[i, :] -= error[i, 0] * (PxR[:, 0]/RxPxR)

    def save(self, filename):
        """
        Saves the network into a .npz file.
        """
        np.savez(
            filename,
            Ni = self.Ni,
            N = self.N,
            No = self.No,
            tau = self.tau,
            g = self.g,
            pc = self.pc,
            Io = self.Io,
            delta = self.delta,
            P_plastic = self.P_plastic,
            I = self.I,
            x = self.x,
            r = self.r,
            z = self.z,
            W_in = self.W_in,
            W_rec = self.W_rec,
            W_out = self.W_out,
            W_plastic = self.W_plastic,
            P = self.P,
            P_out = self.P_out
        )

    def load(self, filename):
        """
        Loads the network from a .npz file.
        """
        net = np.load(filename)

        self.Ni = int(net['Ni'])
        self.N = int(net['N'])
        self.No = int(net['No'])
        self.tau = float(net['tau'])
        self.g = float(net['g'])
        self.pc = float(net['pc'])
        self.Io = float(net['Io'])
        self.delta = float(net['delta'])
        self.P_plastic = float(net['P_plastic'])
        self.N_plastic = int(self.P_plastic*self.N)
        self.I = net['I']
        self.x = net['x']
        self.r = net['r']
        self.z = net['z']
        self.W_in = net['W_in']
        self.W_rec = net['W_rec']
        self.W_out = net['W_out']
        self.W_plastic = net['W_plastic']
        self.P = net['P']
        self.P_out = net['P_out']
