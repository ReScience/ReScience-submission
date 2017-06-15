## Code repository

The Jupyter notebook replication.ipynb contains all scripts to produce the figures. 

Since the simulations take some time, the file fixation_probabilities.dat contains saved data that can be loaded in the notebook to quickly plot and inspect results without running the simulations. 

The two Matplotlib style files style1plot.mplstyle and style2plots.mplstyle can be used to make figures with the right dimensions for the ReScience article. By default, they are not loaded.

The cython file algorithms.pyx contains the core algorithms for the simulations. 
The first code cell in the notebook will try to compile and import the cython-code automatically. If this fails with an error like "DistutilsExecError: command 'gcc' failed with exit status 1" try to compile it manually:
* Generate the C-code with: ```cython algorithms.pyx```
* Compile it with:
   ```
    gcc -shared `python3-config --cflags` `python3-config --ldflags` \
         -I`python -c "import numpy; print(numpy.get_include())"` \
         -o algorithms.so algorithms.c
   ```
* In the first cell of the Jupyter notebook, remove the line 
    "import pyximport; pyximport.install()".


## Requirements

* Python 2 or 3
* Numpy
* Matplotlib
* NetworkX
* ipywidgets (recommended, but not necessary)
* Cython
