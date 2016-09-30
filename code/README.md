## Introduction

This is a reference implementation of the following model:

  A. Compte, M.V. Sanchez-Vives, D.A. McCormick, X.-J. Wang, "Cellular and 
  network mechanisms of slow oscillatory activity (<1 Hz) and wave propagations
  in a cortical network model", J. Neurophysiol. (2003), 2707-2725.


## NEST installation

For neuronal network simulation, we use the publicly available simulation platform NEST (v.2.8.0).
Extensive documentation can be found on the official website (http://nest-initiative.org/).

We here describe how to set up your NEST installation with the re-implemented neuron model of Compte et al. (2003) in Ubuntu.


1. Download and unpack NEST source code from https://zenodo.org/record/32969#.V44c4-3nhC1 to '\<path-to-nest\>/nest-2.8.0/' folder  
2. Copy (and replace, when needed) the following files to '\<path-to-nest\>/nest-2.8.0/models/' folder:  
    'compte2003_ex.cpp'   
    'compte2003_ex.h'   
    'compte2003_in.cpp'   
    'compte2003_in.h'   
    'modelsmodule.cpp'  
    'Makefile.am'  
  
3. Copy and replace the following files to '\<path-to-nest\>/nest-2.8.0/nestkernel/' folder:  
    'nest_names.cpp'  
    'nest_names.h'  

4. NEST needs a few third party tools and libraries to work. Make sure the required NEST dependencies are 
installed (see section "standard configuration" in http://www.nest-simulator.org/installation/ )

5. Install NEST as follows (see http://www.nest-simulator.org/installation/):  
    a. From '\<path-to-nest\>/nest-2.8.0/' run ./bootstrap.sh  
    b. In parallel, create a build directory '\<path-to-nest\>/nest-2.8.0-build/'  
    c. Type in console from build directory:  
```
$../nest-2.8.0/configure --prefix=<path-to-nest>/install-2.8.0
$make
$make install
```
   

6. To make your NEST installation visible from Python by default from the console, add the following lines to your .bashrc file:  
export PATH=$PATH:\<path-to-nest\>/install-2.8.0/bin  
export PYTHONPATH=$PYTHONPATH:\<path-to-nest\>/install-2.8.0/lib/pythonX.Y/site-packages  

where pythonX.Y is your python version.


## Python installation
The network model and further analysis are implemented with Python (v.3.5.2). Scripts were tested with numpy (v.1.11.1) and matplotlib (v.1.5.1).  


To install Python 3, type in console:  
```
$sudo apt-get update 
$sudo apt-get install python3
```


In case of having readline issues while running the code with anaconda, type:  
```
$conda remove --force readline && pip install readline
```

One might also need to install texlive-latex-extra and dvipng packages, if not installed:
```
$sudo apt-get update
$sudo apt-get install texlive-latex-extra
$sudo apt-get install dvipng
```


## Running the scripts
To run the model and perform the analysis in the regime, where activity emerges
spontaneously, type from the console:
```
$python Compte2003_spont.py
```
Simulation takes approx. 11 minutes. Also, one can set "flag_ItoE = True" or 
"flag_EtoI = True" to enchance i-to-E or reduce E-to-I synaptic weights by 10% 
correspondingly.  


To run the model and perform the analysis in the regime, where activity is 
evoked by external stimulation, type from the console:
```
$python Compte2003_stim.py
```
Simulation takes approx. 5 minutes  


To run the model and perform the analysis in the regimes, when AMPA, NMDA, or 
GABA channels are blocked, type from the console:
```
$python Compte2003_block.py
```
Simulation takes approx. 25 minutes. The network size is half of original 
to reduce the simulation time.  



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
$python Compte2003_R_method.py
```

