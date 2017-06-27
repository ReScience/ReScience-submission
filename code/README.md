## Code repository

### Overview

The Jupyter notebook replication.ipynb contains all scripts to produce the figures. 

Since the simulations take some time, the file fixation_probabilities.dat contains saved data that can be loaded in the notebook to quickly plot and inspect results without running the simulations. 

The cython file algorithms.pyx contains the core algorithms for the simulations. 

The two Matplotlib style files style1plot.mplstyle and style2plots.mplstyle can be used to make figures with the right dimensions for the ReScience article. By default, they are not loaded.


### Requirements

* Python 2 or 3
* Numpy
* Matplotlib
* NetworkX
* ipyparallel (recommended, to run simulations in parallel)
* ipywidgets (recommended, to display progress bars during simulations)
* Cython (only needed if you modify the cython-code)


### Compilation

To build the python module ''algorithms'' run:
```python setup.py build_ext --inplace ```

By default, this will use the provided algorithms.cpp. If you modify the cython code in algorithms.pyx, set USE_CYTHON = True in setup.py and run the command above. This will regenerate algorithms.cpp and build the python module.    



