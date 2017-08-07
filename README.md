
# Chemical Signalling in the Neurovascular Unit

This is a replication of the following article:

Witthoft A, Karniadakis GE (2012) A bidirectional model for communication in the neurovascular unit. *Journal of Theoretical Biology* 311: 80-93.

where the authors present *"an interactional model of bidirectional signalling in the neurovascular unit"*. This model combines previous models of the astrocyte signalling cascade [1,2] with a model of smooth muscle cell contractions in small arteries [3]. *"This is the first computational model of astrocyte response to vascular function, making it the first model of a neurovascular unit to include a two-way communication path between the brain and vasculature."*


## References

[1] Bennett MR, Farnell L, Gibson W (2008) Origins of blood volume change due to glutamatergic synaptic activity at astrocytes abutting on arteriolar smooth muscle cells. *Journal of Theoretical Biology* 250: 172-185.

[2] Farr H, David T (2011) Models of neurovascular coupling via potassium and EET signalling. *Journal of Theoretical Biology* 286: 13-23.

[3] Gonzalez-Fernandez JM, Ermentrout B (1994) On the origin and dynamics of the vasomotion of small arteries. *Mathematical Biosciences* 119: 127-167.


## Pre-requisites

This replication has been written and tested on Linux Mint 18.1 Cinnamon using the
following packages:

 * Python 3.5.3
 * Numpy 1.12.1
 * Scipy 0.19.0
 * Matplotlib 2.0.2
 
Original data is in the data directory.


## Usage

Four different simulation scenarios have been implemented using the ODE describing chemical transport and blood vessel mechanics in the NVU. To run a "control" simulation containing stimulus presentation from the synaptic space run

```
 main.py ../data/parameter.cfg
```

To run a simulation of the NVU under manual stretching of the blood vessel run

```
 main_trpv.py ../data/parameter.cfg
```

To run a simulation of the NVU under administration of the vasodilatory drug pinacidil run

```
 main_k.py ../data/parameter.cfg
```

To run a simulation of the NVU with fixed values for the astrocyte membrane potential run

```
 main_Vk.py ../data/parameter.cfg
```
