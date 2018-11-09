# Reproduction of: Fast-Activating Voltage- and Calcium-Dependent Potassium (BK) Conductance Promotes Bursting in Pituitary Cells: A Dynamic Clamp Study

A reference implementation of
"*Fast-Activating Voltage- and Calcium-Dependent Potassium (BK) Conductance
Promotes Bursting in Pituitary Cells: A Dynamic Clamp Study*,
J. Tabak, M. Tomaiuolo, A. Gonzalez-Iglesias,  L. Milescu and R. Bertram,
Journal of Neuroscience 31.46 (2011), 10.1523/JNEUROSCI.3235-11.2011"



## Docker environment

We have created a [Docker](https://www.docker.com/) environment
with all dependencies installed.
This Docker environment can be started, and the `code` and `article` directory mounted
by running the bash script `run_docker.sh` from within the this directory.
All results have been created in this Docker environment.


## Dependencies

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

## Content

The content is:

* `article/` - Contains the description of the reproduction.
* `code/` - Contains the code for reproducing the results.


## Running the code

The code is found in the `code/` directory.
To create Figure 1 and Figure 2 in Tabak et al. run from `code/`:

```
python analysis.py
```

This runs in parallel and takes around 10 hours on a workstation computer.

This reproduction can be speed up by increasing the time step `dt` in `tabak.py`.
Setting `dt = 0.25` gives results similar to the results in the paper,
and only takes around 23 minutes.

To perform the uncertainty quantification and sensitivity analysis of the model
run from `code/`:

```
python uq.py
```

This takes around 20 minutes on a workstation computer.

