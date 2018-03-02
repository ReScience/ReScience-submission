## Introduction

This is a reference implementation of the following model:

  T.C. Potjans and M. Diesmann. “The Cell-Type Specific Cortical Microcircuit: 
  Relating Structure and Activity in a Full-Scale Spiking Network Model”. In: 
  Cerebral Cortex 24.3 (2014), pp. 785–806.

## Platform information

**Platform:** Ubuntu

**gcc (GCC):** 5.4.0

**Python:** 2.7.5 

**Brian:** 2.0.1 

**Matplotlib:** 2.0.0

**NumPy:** 1.11.1

**Pandas:**  0.19.2

**SciPy:** 0.19.0

**Joblib:** 0.11

**Multiprocessing:** 0.70a1

**Weave** 0.15.0

**Cython** 0.23.4

The network model is implemented using the simulator for spiking neural networks Brian 2, which is written in Python.
For data processing and visualization we use Matplotlib, NumPy, Pandas and Scipy.
To run multiple simulations in parallel we use Joblib and Multiprocessing.

## Packages installation

**To assure that all data and figures will be generated as presented in the article, we recommend to follow the instructions below and use the modules in the same versions described here. Although the code also works in newer versions, for Pandas v.0.2.0 (or newer) some colors are not displayed as presented in the article (e.g. figure 2 barplots might be displayed only in blue).**

### Python installation
The network simulation is implemented with Python (v.2.7.5). Although, there is no restriction for more recent versions (Tested in Python v.3.6.3).

To install Python 2, type in console:

```
$sudo apt-get update 
$sudo apt-get install python2.7
```

### Installing pip

We use pip, a package management system, to install the Python modules described above.
To install pip in Ubuntu type in a terminal:

```
sudo apt-get install python-pip
```

Upgrade pip to the latest version:

```
pip install --upgrade pip
```

Installation of packages using pip can be done with the following command:

```
pip install --user PYTHON_PACKAGE_NAME
```

#### Python modules installation using pip (recommended)

To install the required packages type in terminal:

```
pip install --user brian2
pip install --user matplotlib==2.0.0
pip install --user numpy==1.11.1
pip install --user pandas==0.19.2
pip install --user scipy==0.19.0
pip install --user joblib==0.11
pip install --user weave==0.15.0
pip install --user cython==0.23.4
```
or

```
pip install --user brian2 matplotlib==2.0.0 numpy==1.11.1 pandas==0.19.2 scipy==0.19.0 joblib==0.11
```

All software packages are also available in the anaconda distribution, and in brian-team channel (see below).

### Alternative installation (using Anaconda)

Alternatively, you can install the scientific Python modules with the Anaconda data science platform.

For 32 bits system use the download link: https://repo.continuum.io/archive/Anaconda2-5.0.1-Linux-x86.sh
For 64 bits system use the download link: https://repo.continuum.io/archive/Anaconda2-5.0.0-Linux-ppc64le.sh

To install open a terminal in the folder containing the downloaded file and type:

```
bash Anaconda2-5.0.1-Linux-x86.sh
```

for 32 bits.

or

```
bash Anaconda2-5.0.0-Linux-ppc64le.sh
```

for 64 bits.

#### Python modules installation using Anaconda

Brian 2 and Joblib are not included in Anaconda by default.

To install brian-team channel, necessary to install Brian 2, open a terminal and type:

```
conda install -c brian-team brian2
```

For further information access the website: http://brian2.readthedocs.io/en/stable/introduction/install.html

For Joblib install:

```
conda install joblib=0.11
```

The other useful packages can be installed with:

```
conda install PYTHON_PACKAGE_NAME
```


## Code repository

This folder contains seven Python codes:
  *  **netRun.py:** Main script to run the simulations. Here is defined the protocols to run the experiments.
  *  **netPD.py:** Assemble parameters of the network model.
  *  **netParams.py:** Networks structure parameters.
  *  **netModels.py:** Function responsible to make connections between different populations of neurons and to connect stimulation (background or DC currents).
  *  **neuronModels.py:** Neuron model equations and parameters.
  *  **figures.py:** Functions to build figures 2, 5 and 6 from the article.
  *  **ks_test.py:** Script to build figures 3 and 4 from the article.


## Running the scripts

The main script used to simulate the network and generate figures [], can be run typing the following command in a terminal:

```
python netRun.py ARG1 ARG2
```

where ARG1 and ARG2 are command line arguments, passed to the script to specify the experimental protocol and the simulation time in seconds, respectively.
The experimental protocols are specified below:

protocol = 0:   Spontaneous activity (figure 2)
protocol = 1:   DC input and layer-independent experiments (figures 5A and 5B)
protocol = 2:   Layer-independent randomized to generate histograms in figure 5C
protocol = 3:   Dependence of network activity on the background firing rate and the relative inhibitory synaptic strength (figure 6)
protocol = 4:   comparison of spontaneous activity using equations 3 or 4 from paper to calculate the number of synapses between populations
protocol = 5:   response to transient thalamic input (figure 7)

After running the simulations, every data generated will be saved in the data folder (../data/) and the figures will be saved in the figures folder (./figures/).

**WARNING**: It is required around 13GB of RAM for each network to be simulated.

Running the scripts takes approx. 4 minutes to simulate 1 second of network activity.

To run protocol 0 for 60 seconds type in the console:

```
python netRun.py 0 60.0
```

**Obs**: there is the possibility to run the code in parallel using openmp implemented by Brian 2. To do so uncomment lines 32-33. Nevertheless, all simulations were done without it.

To run protocol 1 for 60 seconds type in the console:

```
python netRun.py 1 60.0
```

To run protocol 2 for 5 seconds type in the console:

```
python netRun.py 2 5.0
```

To run protocol 3 for 10 seconds type in the console:

```
python netRun.py 3 10.0
```

Observe that protocols 2 and 3 vary over several parameters and this way cost some time.
Instead of running the simulation in serial which is default, one can run the range of parameters in parallel by simply attributing **False** to the variable **serial** (line 83 in netRun.py).
The number of cores can be set in the variable **num_cores** (line 84 in netRun.py).
**Obs**: in this case the memory cost is multiplied by num_cores.

To run protocol 4 for 60 seconds type in the console:

```
python netRun.py 4 60.0
```

**Obs**: to compare data from NEST it is necessary to run the original code available at http: //www.opensourcebrain.org/projects/potjansdiesmann2014.

To run protocol 5 for 100 seconds type in the console:

```
python netRun.py 5 100.0
```
