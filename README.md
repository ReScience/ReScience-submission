
##  [Re] Spike Synchronization and Rate Modulation Differentially Involved in Motor Cortical Function

**Authors:** Vahid Rostami, Junji Ito, Michael Denker, Sonja Gr&uuml;n

**Corresponding author:** Vahid Rostami, v.rostami@fz-juelich.de

**A reference implementation of:**

Spike synchronization and rate modulation differentially involved in motor cortical function. Alexa Riehle, Sonja Gr&uuml;n, Markus Diesmann, and Ad Aertsen (1997) Science 278:1950-1953. DOI:10.1126/science.278.5345.1950


In this paper we illustrate the successful reproduction of the results shown in Riehle et al (1997) using our new Python implementation of the Unitary Events (UE) method. Our implementation of the UE analysis is available in the `unitary_event_analysis` module of the the Electrophysiology Analysis Toolbox [Elephant](https://github.com/NeuralEnsemble/elephant).

## Dependencies

- Elephant>=0.5.0
- neo>=0.4.0
- quantities>=0.9.0
- numpy>=1.6.2
- matplotlib>=1.5.1

## Installation
### On Ubuntu/Debian:
- sudo apt-get install python-numpy python-matplotlib python-pip ipython
- pip install quantities
- pip install elephant



## Structure of this repository

### article
This folder contains the accompanying text in markdown/pdf/tex format and all the reproduced figures (used in our paper).

### code
This folder contains all the plotting functions, and the functions for loading and converting (from `gdf` format to `Neo`) the data (provided in `utils.py`). 

The source code of our UE implementation is pull requested and accepted after peer-review [https://github.com/NeuralEnsemble/elephant/pull/64](https://github.com/NeuralEnsemble/elephant/pull/64) in [Elephant](https://github.com/NeuralEnsemble/elephant).

### data
This folder contains the preprocessed versions of the spike train data used in the original publication (provided by Dr. Alexa Riehle, CNRS-AMU, Marseille).

### notebook
An jupyter notebook is provided here to plot the figures (Figure 1, Figure 2 and Figure 5) of our paper.

how to use the jupyter notebook:
- https://jupyter.readthedocs.io/en/latest/index.html

how to install the jupyter notebook:
- pip install jupyter


