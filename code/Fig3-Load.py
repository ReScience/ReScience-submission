# Script allowing to reproduce Fig. 3 of:
#
#   Laje, R. and Buonomano, D.V. (2013). Robust timing and motor patterns by taming chaos in recurrent neural networks. Nat Neurosci.
#
# Author: Julien Vitay (julien.vitay@informatik.tu-chemnitz.de)
# Licence: MIT
from __future__ import print_function
import numpy as np

# List of delays
delays = [250, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]

# Load the results from a previous run
data = np.load('../data/timingcapacity.npz')
pearsons = data['r']

# Visualization
import matplotlib.pyplot as plt
correlation_mean = np.mean(pearsons**2, axis=1)
correlation_std = np.std(pearsons**2, axis=1)
plt.errorbar(np.array(delays)/1000., correlation_mean, correlation_std/np.sqrt(10), linestyle='-', marker='^')
plt.xlim((0., 8.5))
plt.ylim((-0.1, 1.1))
plt.xlabel('Interval (s)')
plt.ylabel('Performance ($R^2$)')
plt.show()
