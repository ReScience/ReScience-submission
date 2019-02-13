# Introduction

A reference implementation of
"*Fast-Activating Voltage- and Calcium-Dependent Potassium (BK) Conductance
Promotes Bursting in Pituitary Cells: A Dynamic Clamp Study*,
J. Tabak, M. Tomaiuolo, A. Gonzalez-Iglesias,  L. Milescu and R. Bertram,
Journal of Neuroscience 31.46 (2011), 10.1523/JNEUROSCI.3235-11.2011"



# Docker environment

We have created a [Docker](https://www.docker.com/) environment
with all dependencies installed.
This Docker environment can be started, and the `code` and `article` directory mounted
by running the bash script `run_docker.sh` from within this directory.
All results have been created in this Docker environment.


# Dependencies

The required dependencies are:

* `numpy`
* `matplotlib`
* `uncertainpy`
* `chaospy`
* `tqdm`
* `NEURON`

These can be installed with:

```
pip install numpy
pip install matplotlib
pip install uncertainpy
pip install chaospy
pip install tqdm
```

Additionaly the [Neuron](https://www.neuron.yale.edu/neuron/download) simulator
with the Python interface is required. NEURON must be manually installed
by the user.

# Running the code

To create Figure 1 and Figure 2 in Tabak et al. run:

```
python analysis.py
```

This runs in parallel and takes around 10 hours on a workstation computer.

This reproduction can be speed up by increasing the time step `dt` in `tabak.py`.
Setting `dt = 0.25` gives results similar to the results in the paper,
and only takes around 23 minutes.

To perform the uncertainty quantification and sensitivity analysis of the model
run:

```
python uq.py
```

This takes around 20 minutes on a workstation computer.

To examine the effect of different timestep and area run:

```
python difference.py
```

This takes around 60 hours on a workstation computer.


# Content

The content of this folder is:

* `tabak.py` - contains the model implementation.
* `*.mod` - NEURON files that implements the various ion channels.
* `burstines.py` - contains the functions for calculating the burstiness.
* `analysis.py` - contains the original analysis of the model and recreates Figure 1 and 2 from the original publication.
* `uq.py` - contains the uncertainty analysis and parameter exploration.
* `difference.py` - contains the examination of the timestep and area.


# Platform and package specifications

All results have been generated inside a Docker environment with:

```
Platform: linux
Python: 3.7.0 (default, Jun 28 2018, 13:15:42)
[GCC 7.2.0]
Machine and architecture x86_64 64bit
NumPy: 1.15.2
SciPy: 1.1.0
matplotlib: 3.0.0
Chaospy: 2.3.5
Uncertainpy: 1.1.4
NEURON: 7.6.2-3-g9f36b13
tqdm: 4.28.1
```