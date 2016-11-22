# Requirements

 - Python 2.7 or 3.4
 - NumPy 1.8
 - Matplotlib 1.3

# Description

The `augment.py` module contains classes describing the AuGMEnT model,
the saccade/anti-saccade task with or without shaping strategy and the probabilistic decision making task.

The `simulation.py` executable runs multiple networks training and prints the success rate and median convergence time.
**_Warning: simulations can be heavily time and processing power consuming._**

```
usage: python simulation.py [-h] [-n NETWORKS]
                     [-t {saccade,saccadenoshaping,probabilistic}]
                     [-p PROCESSORS]

Train AuGMEnT networks and print some statistics

optional arguments:
  -h, --help            show this help message and exit
  -n NETWORKS, --networks NETWORKS
                        number of networks to train
  -t {saccade,saccadenoshaping,probabilistic}, --task {saccade,saccadenoshaping,probabilistic}
                        task to run
  -p PROCESSORS, --processors PROCESSORS
                        number of processors to use
```

The `plot-activation.py` executable plots the activity of Q-value units in trained networks.

```
usage: python plot-activation.py [-h] [-t {probabilistic,saccade}]

Plot Q-value units activation of a trained network

optional arguments:
  -h, --help            show this help message and exit
  -t {probabilistic,saccade}, --task {probabilistic,saccade}
                        task to run
```