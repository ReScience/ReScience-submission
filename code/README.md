## Introduction

This is a reference implementation of the following model:

  A. Compte, M.V. Sanchez-Vives, D.A. McCormick, X.-J. Wang, "Cellular and 
  network mechanisms of slow oscillatory activity (<1 Hz) and wave propagations
  in a cortical network model", J. Neurophysiol. (2003), 2707-2725.



## NEST installation

For neuronal network simulation, we use the publicly available simulation platform NEST (v.2.8.0).
Extensive documentation can be found on the official website (http://nest-initiative.org/).

We here describe how to set up your NEST installation with the re-implemented neuron model of Compte et al. (2003) in Ubuntu:  
1. Download and unpack NEST source code from http://www.nest-simulator.org/download/ to '<path-to-nest>/nest_2.8.0/' folder  
2. Copy (and replace, when needed) the following files to '<path-to-nest>/nest-2.8.0/models/' folder:  
    'compte2003_ex.cpp'   
    'compte2003_ex.h'   
    'compte2003_in.cpp'   
    'compte2003_in.h'   
    'modelsmodule.cpp'  
    'Makefile.am'  
  
3. Copy and replace the following files to '<path-to-nest>/nest-2.8.0/nestkernel/' folder:  
    'nest_names.cpp'  
    'nest_names.h'  

4. Install NEST as follows (see http://www.nest-simulator.org/installation/):  
    a. From '<path-to-nest>/nest-2.8.0/' run ./bootstrap.sh  
    b. Create a build directory '<path-to-nest>/nest-2.8.0-build/'  
    c. Type in console from build directory:  
```
$../nest-2.8.0/configure --prefix=$HOME/<path-to-nest>/install-2.8.0
$make
$make install
```
    
5. To make your NEST installation visible from Python by default from the console, add the following lines to your .bashrc file:  
export PATH=$PATH:<path-to-nest>/install-2.8.0/bin  
export PYTHONPATH=$PYTHONPATH:<path-to-nest>/install-2.8.0/lib/python2.7/site-packages  


## Running the scripts
The network model and further analysis are implemented with Python (v.2.7.6). Scripts were tested with numpy (v.1.10.1) and matplotlib (v.1.5.0).
To run the model and perform the analysis, type from the console:
```
$python Compte2003.py
```

To see the neuronal response to DC injection into the soma, type:
```
$python Compte2003_bm_neuron.py
```
To see the comparison of original and simplified synaptic implementation, type:
```
$python Compte2003_bm_syn.py
```
To see the analysis of neuronal membrane resistance, type:
```
$python Compte2003_IV.py
```
To see a schematic representation of the virtual hyperpolarization method, type:
```
$python Compte2003_R_scheme.py
```
