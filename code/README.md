### Code repository

Author: René Larisch (rene.larisch@informatik.tu-chemnitz.de)

License: GNU GPL v2

#### Dependencies

* Python v2.7 or >= v3.4
* Numpy v1.11.0
* Matplotlib v1.5.1
* ANNarchy >= v4.6.4

##### Install ANNarchy
The model reimplementation requires the neuron simulator ANNarchy.
The source can be found online (https://bitbucket.org/annarchy/annarchy).

ANNarchy needs following packages.
* g++ >= 4.6 (4.7 or above is recommended)
* make >= 3.0
* python == 2.7 or >= 3.3
* cython >= 0.19
* setuptools >= 0.6
* numpy >= 1.8
* sympy >= 0.7.4
* scipy >= 0.12
* matplotlib >= 1.3.0

ANNarchy is available for the Python package manager PIP, to install it easy:

```
pip install ANNarchy

```

Or download the actual source code and install with administrator permissions:
```
sudo python setup.py install

```
Or in the home directory with:
```
python setup.py install --user

```

ANNarchy is available for GNU/Linux distributions and MacOS X (with limitations).


The developers provide a good documentation with different examples (https://annarchy.readthedocs.io/en/latest/index.html).

#### Python Scripts

  * **net_fix.py:** Model description with the ANNarchy framework and a static homeostatic mechanism.
  * **net_homeostatic.py:** Model description with the ANNarchy framework and a dynamic homeostatic mechanism.
  * **ann_clamp.py:** Protocol for the voltage clamp experiment.
  * **ann_window.py:** Protocol for the classical STDP learning window.
  * **ann_pairing.py:** Protocol pairing repetition experiment.
  * **ann_burst.py:** Protocols for the burst spiking experiments.
  * **ann_temporalCode.py:** Protocol for the connection patterns, depending on spiking order.
  * **ann_rateCode_Poisson.py:**: Protocol for the connection patterns, depending on firing rate.
  * **ann_stable.py:** Protocol for emergence of stable weights by presenting random input.
  * **ann_RF.py:** Protocol for emergence of receptive fields, similar to them of simple cells in the primary visual cortex, by presenting natural scenes. Please not, that the image data set from Olshausen and Fields (1996) is required for this experiment. The Matlab file with the images must be in the same directory as this python script.
  * **cmap.py:** Define of a custom matplotlib colormap for the connectivity experiments, similar to them in the original publication.

  To test the reimplementation of the network, start one of the experiment protocols, for example:
  ```
  python ann_burst.py

  ```
