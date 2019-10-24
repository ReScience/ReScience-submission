### Code repository

Author: René Larisch (rene.larisch@informatik.tu-chemnitz.de)

License: GNU GPL v2

#### Dependencies
The reimplementation are tested on the following software packages:

* Python v2.7 and Python v3.6
* Numpy v1.11.0
* Matplotlib v1.5.1
* ANNarchy v4.6.8.1

##### Install ANNarchy

The model reimplementation requires the neuro-simulator ANNarchy (Artificial Neural Networks architect).
The code is open source (GPLv2+) and can be found online (https://bitbucket.org/annarchy/annarchy).

ANNarchy requires the following packages.

* g++ >= 4.8
* make >= 3.0
* python == 2.7 or >= 3.3
* cython >= 0.19
* setuptools >= 0.6
* numpy >= 1.8
* sympy >= 0.7.4
* scipy >= 0.12
* matplotlib >= 1.3.0

ANNarchy is available for the Python package manager `pip`, to install the last stable version:

```
pip install ANNarchy
```

One can also download the source code from bitbucket and install it with:

```
python setup.py install
```

or

```
python setup.py install
```
to install it in the home directory

ANNarchy is available for GNU/Linux distributions and MacOS X (with limitations).

The developers provide an extensive documentation with different examples (https://annarchy.readthedocs.io/en/latest/index.html).

#### Python Scripts

**Network definition**

* **network.py:** Model description with the ANNarchy framework.

**Individual experiments**

* **Fig1_clamp.py:** Protocol for the voltage clamp experiment.
* **Fig1_window.py:** Protocol for the classical STDP learning window.
* **Fig1_pairing.py:** Protocol pairing repetition experiment.
* **Fig2_burst.py:** Protocols for the burst spiking experiments.
* **Fig3_temporalCode.py:** Protocol for the connection patterns, depending on spiking order.
* **Fig3_rateCode_Poisson.py:**: Protocol for the connection patterns, depending on firing rate.
* **Fig3_temporalCode_stand.py:** Protocol for the connection patterns, depending on spiking order with a standard STDP learning rule.
* **Fig3_rateCode_Poisson_stand.py:**: Protocol for the connection patterns, depending on firing rate with a standard STDP learning rule.
* **Fig4_stableW.py:** Protocol for emergence of stable weights by presenting random input.
* **Fig4_RF.py:** Protocol for emergence of receptive fields, similar to them of simple cells in the primary visual cortex, by presenting natural scenes. Please not, that the image data set from Olshausen and Fields (1996) is required for this experiment (can be found here: https://www.rctn.org/bruno/sparsenet/IMAGES.mat). The Matlab file with the images must be in the same directory as this python script.
* **startAnalysis.py:** Script to call all the protocol-scripts after each other.
Please note, at the end of each python script, a image is showing.
Close the image to terminate the script and to start the next one.

**Utilities**

* **cmap.py:** Define of a custom matplotlib colormap for the connectivity experiments, similar to the one used in the original publication.

To test the reimplementation of the network, start one of the experiment protocols, for example:

```
python ann_burst.py
```
or run all experiment protocols after each other:
```
python startAnalysis.py
```
