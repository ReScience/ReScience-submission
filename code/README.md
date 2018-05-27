## Running the code: quick guide

To create a new environment with all necessary packages called HathwayGoodman, in a terminal and in the *code* directory type

```bash
$ conda update conda
$ conda env create -f environment.yml
```

Activate the environment with

```bash
$ source activate HathwayGoodman
```

Create most figures, but load data for Figure 5 and 6 to save time (requires the files in the data folder), which should take ca. 5min. The figures are saved in article/figures for easy comparison with the figures in the text.

```bash
$ python main.py 
```

Create all figures from scratch, but use less repetitions per parameter combination (10 instead of 100) than in the paper, which should take ca. 8 hours.  The figures are saved in article/figures for easy comparison with the figures in the text.

```bash
$ python main.py --new True
```

Running all 25 sets of 100 repetitions for Fig 5 (22 sets) and for Fig 6 (3 sets) is best done on a cluster since running all sets with all 100 repetitions would take more than a day even when running on multiple cores. More details are provided below.



Notes:

Figure 7CD might not look like in the paper because the all-to-all rule is not as reliable at finding the pattern. Therefore it might be necessary to run figure_7CD.py a few times to see a similar figure as in the text.

 



## Running the code: in detail

This implementation is written in Python and uses the spiking neural network simulator Brian 2, which also requires `numpy`, `sympy`, `cython`, `pyparsing`,`scipy`,`jinja2`, as well as numba and matplotlib.

The packages required to run this implementation and the versions on which they were tested on are:

- brian2=2.1.2
- matplotlib=2.2.2
- numba=0.38.0
- numpy=1.13.3
- python=3.6.5



#### Installing packages

We advise to create a seperate environment to make sure that package versions match and the code runs successfully. 

To install packages manually, one can use Anaconda 

```bash
$ conda install -c brian-team 'brian2>=2.1.2' 'matplotlib>=2.2.2' 'numba>=0.38.0' 'numpy>=1.13.3'
```

or the Python installer pip

```bash
$ pip install 'brian2>=2.1.2' 'matplotlib>=2.2.2' 'numba>=0.38.0' 'numpy>=1.13.3'
```



The source code was tested on 

- Ubuntu 16.4
- macOS 10.13.4 (High Sierra)



#### Available files

All figures can be created independently by running

**main.py** : this script creates all figures. if run without the "--new True" flag, Figures 5 and 6 will be created using the saved data in the data folder.
Usage: python main.py
Usage: python main.py --new True

**figure_1.py** : creates figure 1, can be run independently from other figures

**figure_2.py** : creates figure 2, can be run independently from other figures

**figure_3_4_7AB.py**  : cannot be run independently, will be called by run_simulation to create figures 3, 4 and 7AB

**figure_5.py** : creates figure 5, by default running this script will run 10 repetitions of all 22 parameter combinations. can be run independently from other figures

**figure_6.py** : creates figure 6, by default running this script will run 10 repetitions of 3 different time step sizes. can be run independently from other figures

**figure_7CD.py** : creates figure 1, can be run independently from other figures. Note that this figure might differ from the one in the text due to the all-to-all STDP rule beind less reliable at finding the pattern.

**figure_8.py** : creates figure 1, can be run independently from other figures

**run_simulation.py** : runs the main pattern finding algorithm and creates figures 3, 4, 7AB, can be run independently from other figures

**create_input.py** : is called by run_simulation.py to create the spike trains. can be run independently

**write_param_inputfile.py** : writes the parameter file that stores the combinations of parameters used for Figure 5.

The parameter combinations are stored in para.npy, which has the structure

```
row     win     jit     # pat neur      pat freq        % deleted
0       0.475   1       1000            0.25            0.0

1       0.275   1       1000            0.25            0.0
2       0.325   1       1000            0.25            0.0
3       0.375   1       1000            0.25            0.0
4       0.425   1       1000            0.25            0.0

5       0.475   0       1000            0.25            0.0
6       0.475   2       1000            0.25            0.0
7       0.475   3       1000            0.25            0.0
8       0.475   4       1000            0.25            0.0
9       0.475   5       1000            0.25            0.0
10      0.475   6       1000            0.25            0.0

11      0.190   1       400             0.25            0.0
12      0.285   1       600             0.25            0.0
13      0.380   1       800             0.25            0.0
14      0.570   1       1200            0.25            0.0

15      0.475   1       1000            0.05            0.0
16      0.475   1       1000            0.10            0.0
17      0.475   1       1000            0.20            0.0
18      0.475   1       1000            0.50            0.0

19      0.4275  1       1000            0.25            0.1
20      0.380   1       1000            0.25            0.2
21      0.3325  1       1000            0.25            0.3
```



#### Running on a cluster 

The files that were used to run on the Imperial College HPC cx1 cluster with a PBS queue system are provided under code/cluster (main_cluster.py and run_main_cluster.sh). Please note that the paths in both files (line 9 in run_main_cluster.sh and lines 14-15 in main_cluster.py) need to be changed. 

The runs are initiated with 

```bash
qsub -J 0-2199 run_main_cluster.sh
```

which runs the parameter of para.npy 100 times with different random seeds.







