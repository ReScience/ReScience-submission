
### ReScience submission repository

This is the submission repository for the [Re**Science** journal](https://rescience.github.io).

**Author:** René Larisch (rene.larisch@informatik.tu-chemnitz.de)

This is a replication of:

Connectivity reflects coding: a model of voltage-based STDP with homeostasis,
C. Clopath, L. Büsing, E. Vasilaki and W. Gerstner, In: Nature Neuroscience 13.3
(2010), pp. 344–352, doi= 10.1038/nn.2479

The reimplementation is written in Python (v2.7 and tested on v3.6) with the help of the neuro-simulator ANNarchy (v4.6.8.1),
Numpy (v1.11.0) and Matplotlib (v1.5.1).

#### Data

The experiment for the emergence of V1 simple-cell-like receptive fields needs the dataset of 10 natural scenes from Olshausen and Field (1996).
Dataset can be found on the webpage of Olshausen (http://www.rctn.org/bruno/sparsenet/). Please note, that the whitened natural scenes are required.

#### Model

The implementation of the model can be found in the **code** directory.
It exists a model version with a fixed and static homeostasic mechanisms (can be found in *code/net_fix.py*),
and a model version with a dynamic homeostasic mechanisms (can be found in *code/net_homeostatic.py*).


#### Results
![Weight change in the Voltage clamp experiment. The blue line presents the weight change with the parameter set for the visual cortex.
The red line presents the weight change with the parameter set for the hippocampus.](article/figures/Fig1_clamp.png)
![Classic STDP learning window. On the x-axis is the time of a postsynaptic spike in relation to the presynaptic spike presented.](article/figures/Fig1_window.png)
![Weight changes as a function of pair frequency repetition.](article/figures/Fig1_pairing.png)
![**Upper left**, weight change as a function of the numbers of postsynaptic spikes.
**Upper right**, weight change as a function of the frequency between three postsynaptic spikes.
**Down**, weight change as a function of the time between the first of three postsynaptic spikes and one presynaptic spike.](article/figures/Fig2_burst.png)
![Neurons with similar high firing rates develop strong bidirectional connections (a). Temporal order of activity is mirrored in the connectivity structure (c). The same as (a) for a standard STDP rule (b). The same as (c) for a standard STDP rule (d).](article/figures/Fig3.png)
![Weight value at the end of the current epoch from the presynaptic neuron to a single postsynaptic neuron by presenting a Gaussian curve input. Black are weight values around zero and weight are weight values around the maximum weight (a). Four different receptive fields by presenting natural scenes as input (b). Higher input firing rate on natural scene input leads to a smaller receptive field (c).](article/figures/Fig4_W.png)
