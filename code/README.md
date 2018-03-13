# Reproducing the results

## Installing dependencies

The simulation code is written based on the neural simulator [Brian 2](http://brian2.readthedocs.io) (which in turn depends on the packages ``numpy``, ``sympy``, ``cython``, ``pyparsing``,``scipy``,``jinja2``) and makes use of the [joblib package](https://pythonhosted.org/joblib/) to parallelize parameter explorations and to save intermediate results (see *Technical notes* below). The figures are created with [matplotlib](http://matplotlib.org/).

The code should run both on Python 2.7 and Python 3 (versions >= 3.4), on Linux, OS X, and Windows (the results included in the paper where run on a Linux machine).

Below, we describe how the dependencies can be installed for users with environments based on the [Anaconda distribution](https://www.anaconda.com/download/#linux), or using Python's standard installation tool [pip](https://pip.pypa.io/en/stable/installing/).

### Using the Anaconda distribution

#### Installing latest versions

```console
$ conda install -c brian-team 'brian2>=2.1.2' 'joblib>=0.11' 'matplotlib>=2.1.0' 'scipy>=1.0'
```

#### Installing frozen "known-to-work" versions

```console
$ conda env create -f environment.yml
$ source activate ReScience_CazeStimbergGirard  # on Linux / OS X
$ activate ReScience_CazeStimbergGirard  # on Windows
```

### Using pip

#### Installing latest versions

```console
$ pip install 'brian2>=2.1.2' 'cython>=0.26.1' 'joblib>=0.11' 'matplotlib>=2.1.0' 'scipy>=1.0'
```

#### Installing frozen "known-to-work" versions

```console
$ pip install -r requirements.txt
```

## Running the code

To run a reduced set of simulations (with coarser parameter exploration and less trials per parameter setting), navigate to the ``code`` directory and run the ``main.py`` script:

```console
$ python main.py
```

This will run the simulations and produce the figures in the `../article/Figs` directory. It will also report the progress of the simulation. To disable any output, run the script with the `-q` or `--quiet` option.

To run the full set of simulations that were used in the paper (which will take several hours), run:

```console
$ python main.py --long
```

## Technical notes

### Changes to Brian's default settings

Two of Brian's default settings where changed:

* Brian's preference `codegen.target` was explicitly set to `cython`. While this is the default code generation target on Python 3 (if Cython is available), on Python 2 Brian defaults to `weave` (if available). However, `weave`  is not compatible with multiple simulations run in parallel, as used for this paper (see Brian's ["Known Issues page"](http://brian2.readthedocs.io/en/2.1.2/introduction/known_issues.html#parallel-brian-simulations-with-the-weave-code-generation-target))
* By default, Brian logs debugging information to a file that will be deleted when the Python interpreter exists, provided that no error occurred during the run. Given the large number of simulations that are run for the results of this paper, this log file can become very large (several GB) because it repeatedly logs detailed debug output for each of the models/simulations. To avoid this, we added a `brian_preferences` file (see [Brian docs](http://brian2.readthedocs.io/en/2.1.2/advanced/preferences.html#preference-files)) that switches off logging (as explained [here](http://brian2.readthedocs.io/en/2.1.2/advanced/logging.html#preferences))

### Result caching

The code uses joblib's "lightweight pipelining" approach (see [joblib documentation](https://pythonhosted.org/joblib/)), i.e. the code is written as if plotting the results for a figure meant running all the simulations/analysis functions again every time (as opposed to having a separate phases which each phase storing the results to disk). In reality, results from certain functions (annotated with ``@mem.cache`` in the code) will be transparently cached to disk and automatically reused. Running the full exploration a second time (e.g. to change some of the analysis or to change the way the results are plotted) will therefore only take a short time. Note that this cache will be stored in the `code` directory under `__joblib_cache__` and will take up around 5GB of disk space for the full parameter exploration.

### Parallelization

Because the parameter explorations presented in the paper require the run of a large number of short and independent simulations, we decided to make use of the multiple CPU cores a modern computer provide, by parallelizing the simulations (over parameters and/or trials for the same parameters). Such "embarissingly parallel" problems can be easily distributed over independent processes (which will end up on separate CPU cores) using joblib's `Parallel` helper class (see [their documentation](https://pythonhosted.org/joblib/parallel.html)). The script will use a number of parallel jobs equal to the number of CPU cores minus one. To change this number, change the `n_jobs` argument of `Parallel`in the code.
