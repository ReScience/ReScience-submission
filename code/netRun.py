# ----------------------------------------------------------------------------
# Contributors: Renan O. Shimoura
#               Nilton L. Kamiji
#               Rodrigo F. O. Pena
#               Vinicius L. Cordeiro
#               Cesar C. Ceballos
#               Cecilia Romaro
#               Antonio C. Roque
# ----------------------------------------------------------------------------
# References:
#
# *The Cell-Type Specific Cortical Microcircuit: Relating Structure and Activity
# in a Full-Scale Spiking Network Model*,
# Tobias C. Potjans and Markus Diesmann,
# Cerebral Cortex, 24(3):785-806, 2014.
# ----------------------------------------------------------------------------
# File description:
#
# Definition of different sets of variables to run the experiments presented in
# the article.
# ----------------------------------------------------------------------------
import matplotlib
matplotlib.use('Agg')

import netPD as netPD
import figures as fig
from brian2 import *

import sys

# Remove these comments to set the simulation to run with openmp
# set_device('cpp_standalone', directory='parallel_PD')
# prefs.devices.cpp_standalone.openmp_threads = 10

# Function to define filenames
def filename(g, bgrate, suffix = ''):
    return '../data/data_raster_g' + str(g) + '_bgrate' + str(bgrate) +  suffix + '.dat'

# Function used in multiple runs: define different random seeds and clean memory
def runParamsParallel(s=1, g=4, bg_type=0, bg_freq=8.0, stim=0, tsim=1.0, filename=None):
    seed(s*1000)
    netPD.runParams(tsim=tsim, bg_type=bg_type, stim=stim, g=g, bg_freq=bg_freq, filename=filename)
    gc.collect()

###############################################################################
'''
Protocols parameters:

# background type
bg_type = 0:    layer-specific
          1:    layer-independent
          2:    layer-independent-random

# stimulation
stim = 0:       turn on the background noise
       1:       DC current experiment

# max number of synapses between populations
nsyn_type = 0:  no approximation (eq. 3)
nsyn_type = 1:  approximation (eq. 5)

# networks variables changed for different protocols:
g    => inhibitory weight balance
bg   => background rate

# thalamic input
thal = "ON"     turn on the thalamic input transient
thal = "OFF"    turn off the thalamic input transient
'''
###############################################################################

protocol = int(sys.argv[1])
tsim = float(sys.argv[2])   # time of simulation

g = 4.0                     # default value for inhibitory weight balance
bg = 8.0                    # default value for background rate

# default seed of pseudo random numbers to test reproducibility
s = 1000
seed(s)

# choose serial = False to run multiple simulations in parallel
serial = True
num_cores = 8              # number of cores to run in parallel

###############################################################################
# Simulation protocols
###############################################################################
'''
protocol = 0:   spontaneous activity (figure 2)

protocol = 1:   DC input and layer-independent experiments (figures 5A and 5B)

protocol = 2:   layer-independent randomized to generate histograms in figure 5C

protocol = 3:   dependence of network activity on the background firing rate (bg)
                and the relative inhibitory synaptic strength (g) (figure 6)

protocol = 4:   comparison of spontaneous activity using equations 3 or 4 from
                paper to calculate the number of synapses between populations

protocol = 5:   response to transient thalamic input

You can find below the sets of parameters defined to each different experiment.
'''

if protocol==0:
    bg_type = 0
    stim = 0
    netPD.runParams(tsim=tsim, bg_type=bg_type, stim=stim, filename=filename(g, bg, 'default'))
    gc.collect()    #garbage collector to clean memory
    fig.createfig2(tsim, filename(g, bg, 'default'))

elif protocol==1:
    bg_type = 0
    stim = 1
    netPD.runParams(tsim=tsim, bg_type=bg_type, stim=stim, filename=filename(g, bg, '_DC'))
    gc.collect()    #garbage collector to clean memory

    bg_type = 1
    stim = 0
    netPD.runParams(tsim=tsim, bg_type=bg_type, stim=stim, filename=filename(g, bg, '_layer-independent'))
    gc.collect()    #garbage collector to clean memory
    fig.createfig5(tsim, [filename(g, bg, '_DC'), filename(g, bg, '_layer-independent')])

elif protocol==2:
    bg_type = 2
    stim = 0
    if serial == True:
        for s in np.arange(100):
            seed(s*1000)
            netPD.runParams(tsim=tsim, bg_type=bg_type, stim=stim, filename=filename(g, bg, '_bg_random'+str(s)))
            gc.collect()    #garbage collector to clean memory
    else:
        from joblib import Parallel, delayed
        import multiprocessing
        Parallel(n_jobs=num_cores)(delayed(runParamsParallel)(s=s, tsim=tsim, bg_type=bg_type,\
         stim=stim, filename=filename(g, bg, '_bg_random'+str(s)))  for s in range(100))
    fig.createfig5_hist(tsim)

elif protocol==3:
    bg_type = 0
    stim = 0

    g_values = np.arange(2.0, 11.0, 0.5);            # relative inh. synaptic strength values
    bg_values = np.arange(3.0, 15.5,0.5);  # background rate values
    if serial == True:
        for g in g_values:
            for bg in bg_values:
                netPD.runParams(tsim=tsim, bg_type=bg_type, stim=stim, g=g, \
                bg_freq=bg, filename=filename(g, bg))
                gc.collect()    #garbage collector to clean memory
    else:
        from joblib import Parallel, delayed
        import multiprocessing
        Parallel(n_jobs=num_cores)(delayed(runParamsParallel)(tsim=tsim, bg_type=bg_type,\
         stim=stim, g=g, bg_freq=bg, filename=filename(g, bg)) for g in g_values for bg in bg_values)

    fig.datafig6(tsim)
    fig.createfig6()

elif protocol==4:
    bg_type = 0
    stim = 0
    netPD.runParams(tsim=tsim, bg_type=bg_type, stim=stim, nsyn_type=0 ,filename=filename(g, bg, '_noapprox'))
    gc.collect()    #garbage collector to clean memory

    netPD.runParams(tsim=tsim, bg_type=bg_type, stim=stim, nsyn_type=1 ,filename=filename(g, bg, '_approx'))
    gc.collect()    #garbage collector to clean memory

    fig.ks_test(tsim)

elif protocol==5:
    bg_type = 0
    stim = 1
    netPD.runParams(tsim=tsim, bg_type=bg_type, stim=stim, thal='ON', filename=filename(g, bg, '_thal'))
    gc.collect()    #garbage collector to clean memory
    fig.createfig7(tsim, filename(g, bg, '_thal'))
