### Code repository

functions.py - containing all reference implementations of the free water DTI model

run_simulations_1.py - code to generate the article's figure 1 simulation
                    (a notebook version of this file can be found in the notebook folder 
                     and it is entitled as run_simulations_1.ipynb)

run_simulations_2.py - code to generate the article's figure 2 simulation
                    (a notebook version of this file can be found in the notebook folder 
                     and it is entitled as run_simulations_2.ipynb)

run_simulations_3.py - code to generate the article's figure 3 simulation
                    (a notebook version of this file can be found in the notebook folder 
                     and it is entitled as run_simulations_3.ipynb)

run_invivo_data.py - code to generate the article's real data analysis
                    (a notebook version of this file can be found in the notebook folder 
                     and it is entitled as run_data.ipynb)

## Dipy's instalation

As the main article mention, to run these procedures you need a Dipy installation.
For this you can follow the steps below:

1) If you are an OSX users, you have first to install the Apple [Xcode](https://developer.apple.com/xcode/)
developer tools. 

2) To install Dipy's dependencies, we suggest to use [Anaconda](https://www.continuum.io/downloads) distribution and install the
nibabel library by typing in the terminal "pip install nibabel" (step applied to
Windows, OSX and Linux users).

3) After installing the dependencies of step 3, Dipy can be
installed by typing in the terminal "pip install dipy".

More detailed information on Dipy's installation can be found in [dipy's website](http://nipy.org/dipy/installation.html).
If you still have problems on Dipy's installation, after following the information in
[dipy's website](http://nipy.org/dipy/installation.html), please send an e-mail to Dipy's team
to the [nipy mailing list](https://mail.python.org/mailman/listinfo/neuroimaging)
with the subject line starting with [dipy].