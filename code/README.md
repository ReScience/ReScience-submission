## Introduction

A reference implementation of the simple model of Modeling GABA alterations in schizophrenia: a link between impaired inhibition and altered gamma and beta range auditory entrainment, Vierling-Claassen, D., Siekmeier, P., Stufflebeam, S., & Kopell, N., Journal of Neurophysiology, 99(5), 2656-2671.

## Dependencies

It requires python, numpy and matplotlib:
```
$ pip install numpy
$ pip install matplotlib
```

## Running the model
To run a single example model type:
```
python example_run.py
```
This will run a single simulations (see run_example.py for details on the parameters) and produce
some plts of the data.

To run the main set of simulations for the article figures, type:
```
python run_main_simulations.py
```
Note: This might take 15-20 minutes.

Afterwards, in order to replicate the figures from the article, type:
```
python plot_main_figures.py
``` 

To run the set of simulations for the noise explorations, type:
```
python run_noise_exploration_sims.py
```
Note: This might take 15-20 minutes.

Afterwards, in order to replicate the figures from the article, type:
```
python plot_exploration_figures.py
``` 

Alternatively, you can have a look at the notebooks folder and the IPython notebooks there.
