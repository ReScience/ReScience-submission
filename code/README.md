# Description

The `augment.py` module contains classes describing the AuGMEnT model,
the saccade/anti-saccade task with or without shaping strategy and the probabilistic decision making tasks.

The `simulation.py` executable runs multiple networks training and prints the success rate and median convergence time.

This implementation uses Python 3 and NumPy.

*Warning: running a simulation can be heavily time and processing power consuming.*

```
usage: simulation.py [-h] [-n NETWORKS]
                     [-t {saccade,saccadenoshaping,probabilistic}]
                     [-p PROCESSORS]

Train AuGMEnT networks and print median convergence time

optional arguments:
  -h, --help            show this help message and exit
  -n NETWORKS, --networks NETWORKS
                        number of networks to train
  -t {saccade,saccadenoshaping,probabilistic}, --task {saccade,saccadenoshaping,probabilistic}
                        task to run
  -p PROCESSORS, --processors PROCESSORS
                        number of processors to use
```