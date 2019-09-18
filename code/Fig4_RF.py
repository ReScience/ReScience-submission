"""
Python script to reproduce the stable receptive field
learning as in Clopath et al. 2010.
The inputs are from the Olshausen (1996) image data set.
Recive it from https://www.rctn.org/bruno/sparsenet/IMAGES.mat .
The network consists of a Poisson input layer and one post synaptic neuron.
The size of the input layer is determine by the input. As mentioned in the
original publication, the input is a 16x16 pixel size patch, cutted out
randomly of the input dataset. Devided in an ON- and OFF- part.
So the input layer contains 16x16x2 =  512 neurons.

See Fig. 7d in Cloath et al. (2010)
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from ANNarchy import *
setup(dt=1,seed=1001)
import scipy.io as sio
import os
from network import *

# Global parameters
duration = 200#200 #ms presentation time per input patch
s_Patch = 16 # patchsize in pixel
n_patches = 500000 # number of patches to train
maxFR = 70.0 # maximum firering rate

# Presynaptic neuron model
"""
Because of the learning rule, we need a additional layer, that contains the
neccessary variables for the learning. This population is one to one connected
with the Poisson input layer and spiked for ever corresponding neuron in the
Poisson layer.
"""
params = """
EL = -70.4      :population
VTrest = -50.4  :population
taux = 15.0     :population  """
eqs = """
dg_vm/dt = EL/1000 : min = EL, init=-70.4
Spike = if state == 1: 1.0 else: 0.0
dReset/dt = if state == 1: +1 else: -Reset
dxtrace/dt = if state == 1: +1/taux else: -xtrace/taux : init = 0.0
state = if state >0: -1 else: 0"""
neuron = Neuron(parameters = params,
                equations = eqs,
                reset = """ g_vm = EL
                            state = 1""",
                spike = """g_vm > VTrest""")

#Populations
inputPop = PoissonPopulation(geometry=(s_Patch,s_Patch,2),rates=maxFR)
prePop = Population(geometry=(s_Patch,s_Patch,2),neuron=neuron)
postN = Population(geometry=1,neuron=AdExNeuron)

# Projections
projInpt_Pre = Projection(
    pre = inputPop,
    post = prePop,
    target='vm'
).connect_one_to_one(weights = 30)

projInp_N = Projection(
    pre = prePop,
    post= postN,
    target='Exc',
    synapse = ffSyn
).connect_all_to_all(weights = Uniform(0.0,2.0))

projInp_N.set_fix = 0.0 # use the homeostatic mechanisms in the LTD term

def preprocessData(matData):
    # function to split the prewhitened images into on and off counterparts
    images = matData['IMAGES']
    w,h,n_images = np.shape(images)
    new_images = np.zeros((w,h,2,n_images))
    for i in range(n_images):
        new_images[images[:,:,i] > 0, 0, i] = images[images[:,:,i] > 0, i]
        new_images[images[:,:,i] < 0, 1, i] = images[images[:,:,i] < 0, i]*-1

    return(new_images)

def run():
    print('Presenting natural scenes to learn V1 simple cell receptive fields')

    # compile command to create the ANNarchy network
    compile()

    projInp_N.transmit = 1.0
    projInp_N.aLTP = 0.00016*0.3
    projInp_N.aLTD = 0.00014*0.3

    # load input data set
    matData = sio.loadmat('IMAGES.mat')
    images = preprocessData(matData)
    w,h,d,n_img = np.shape(images)

    monW = Monitor(projInp_N,'w',period=5000)

    for p in range(n_patches):
        # ever 20 s make a normalization
        if ((p*duration)%20000)==0:
            wFF = projInp_N.w
            onoff = np.reshape(wFF,(s_Patch,s_Patch,2))
            onNorm = np.sqrt(np.sum(onoff[:,:,0]**2))
            offNorm= np.sqrt(np.sum(onoff[:,:,1]**2))
            onoff[:,:,0] *= offNorm/onNorm
            wFF = np.reshape(onoff,(1,s_Patch*s_Patch*2))
            projInp_N.w = wFF
        # choose randomly an image and position to cout out the patch
        r_img = np.random.randint(n_img)
        xPos = np.random.randint(w-s_Patch)
        yPos = np.random.randint(h-s_Patch)
        maxV = np.max(images[:,:,:,r_img])
        patch = images[xPos:xPos+s_Patch,yPos:yPos+s_Patch,:,r_img]
        # make random vertical and horizontal flip
        if np.random.rand() < 0.5:
            patch = np.fliplr(patch)
        if np.random.rand() < 0.5:
            patch = np.flipud(patch)
        # set the rates for the Poission input population
        inputPop.rates = (patch/maxV)*maxFR
        simulate(duration)

    w = monW.get('w')

    # create the resulting RF out of the input weights
    ff_W = projInp_N.w
    rf = np.reshape(ff_W,(s_Patch,s_Patch,2))
    rf = rf[:,:,0] - rf[:,:,1]
    plt.figure()
    plt.imshow(rf,interpolation='none',cmap=plt.get_cmap('gray'))
    #plt.colorbar()
    plt.savefig('Fig4_RF.png')
    plt.show()
    print("Done with the experiment.")

if __name__ == "__main__":
    if os.path.isfile('IMAGES.mat'):
        run()
    else:
        print('No IMAGES.mat found, please download the file from: https://www.rctn.org/bruno/sparsenet/IMAGES.mat')
