[![Binder](http://mybinder.org/badge.svg)](https://beta.mybinder.org/v2/gh/ThierryMondeel/ReScience-submission-Shlomi2009/Mondeel-Ogundipe-Westerhoff-2017)

## Replication of T. Shlomi, M.N. Cabili, E. Ruppin (2009) "Predicting metabolic biomarkers of human inborn errors of metabolism"

This is a repository with the code and article submitted to the [ReScience journal](http://rescience.github.io).

The reference for the original article:
T. Shlomi, M.N. Cabili, E. Ruppin, Predicting metabolic biomarkers of human inborn errors of metabolism, Mol. Syst. Biol. 5 (2009) 263. [doi:10.1038/msb.2009.22](http://doi.org/10.1038/msb.2009.22).


## Repository structure
The article folder contains the description of our replication of the simulations underlying Figures 1 and 2 of the original article and in addition a brief discussion on the sensitivity of the approach to changes of settings in the method. [Click here to see the PDF](article/Mondeel_Ogundipe_Westerhoff-2017.pdf).

The code folder is empty because we provide all our code in the notebook folder since nearly all relevant code is made available in Jupyter notebooks.
In total we provide 4 notebooks that produce the illustrative network originally considered in Figure 1, analyse the illustrative network, analyse the human metabolic map with respect to 17 amino acid disorders (original Figure 2) and analyse the sensitivity of the approach. 

Our implementation of the biomarker prediction algorithm is also available from the notebook folder as: findBiomarkers.py. 

## How to run the simulations yourself?

The code has been developed on a Mac operating system in Python 3.6. To execute the code on your own computer you will need to install the following Python modules:
COBRApy, NumPy, SciPy, XlsxWriter, tqdm, lxml and xlrd.

Alternatively, use the MyBinder service, through the button at the top of this ReadMe file, to open an anonymous, personal cloud environment where all software is automatically installed and you will be able to run through all the notebooks there. After clicking the "Launch binder" button, you will end up on a page where Binder will attempt to build an environment for you with all required packages. Be patient this might take a minute or so. When completed, a Jupyter interface will be launched with access to all folders of the repository. Navigate to the notebook folder and open any of the notebooks to test and investigate our reproduction. 
