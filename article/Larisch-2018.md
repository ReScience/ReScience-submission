---
Title: "Connectivity reflects coding: a model of voltage-based STDP with homeostasis"
Author:
  - name: René Larisch
    affiliation: 1
Address:
  - code:    1
    address: Professorship for Artificial Intelligence, Department of Computer Science, Chemnitz University of Technology, D-09107 Chemnitz, Germany
Contact:
  - rene.larisch@informatik.tu-chemnitz.de
Editor:
  - Name Surname
Reviewer:
  - Name Surname
  - Name Surname
Publication:
  received:  Feb,  1, 2018
  accepted:  Feb, 1, 2018
  published: Feb, 1, 2018
  volume:    "**1**"
  issue:     "**1**"
  date:      Feb 2018
  number: 1
Repository:
  article:   "http://github.com/rescience/rescience-submission/article"
  code:      "http://github.com/rescience/rescience-submission/code"
  data:      
  notebook:   
Reproduction:
  - "Connectivity reflects coding: a model of voltage-based STDP with homeostasis, C. Clopath, L. Büsing, E. Vasilaki and W. Gerstner, In: Nature Neuroscience 13.3 (2010), pp. 344–352, doi= 10.1038/nn.2479"
Bibliography:
  bibliography.bib

---

# Introduction

Since the first describing of spike timing-dependent plasticity (STDP) [@Bi1998],
different description of STDP are published to reproduce different experimental findings.
The early implementations, so called pair-based STDP learning rules,
failed on reproducing some experimental observations,
like from triplet or quadruplets experiments [@Pfister2006].

Here, we introduce a reimplementation of the @Clopath2010 STDP model,
what is enable to reproduce experimental findings of triplet studies.
The proposed model try to be more biological plausible than previous models with
a theoretical approach, to be a voltage based STDP model.
This means, that the occur of long term depression (LTD) or long term potentiation (LTP) depends on the
postsynaptic membrane voltage.
Clopath and colleagues could show that their learning rule
can develop stable weights, as it is necessary for learning receptive fields of V1 simple cells.
To achieve stable weights, they implemented a homeostatic mechanism to
adjust the amount of generated LTD, over the relation between the average postsynaptic
membrane potential and a reference value.
Furthermore, their model lead to two different connection structures,
depending on the spiking behavior of the neurons [@Clopath2010].
If neurons fire high at the same time, they build strong bidirectional connections (rate code).
If neurons fire in a specific temporary order, they connection structure follow that order (temporal code).
They used a adaptive exponential integrate-and-fire (AdEx) neuron model.



# Methods
## Overview
From the original model exists a matlab implementation to demonstrate the
stable learning of weights on modeldb (http://modeldb.yale.edu/144566).
This model reimplementation is written in Python (v2.7) and with help of the neuronal
simulator ANNarchy [@Vitay2015] (v4.6).
For the analysis and the figures we used numpy (v1.11.0) and matplotlib (v1.5.1).
Not only the voltage based STDP learning rule is reimplemented, even the AdEx neuron model.
In the supplementary material of the original publication is the matlab source code for a simple example
published.

The reimplementation is mainly orientated on the description of neuron model and learning rule in the original publication [@Clopath2010].
Because of the lack on a further description of the homeostatic mechanism and a more precise description of the neuron after a emitted spike,
the available matlab code is the second reference for this reimplementation.

## Model description

The neuron model is mainly borrowed from the description the matlab source code.
Different as mentioned in the publication, after a spike is the membrane potential
set to $29.4mV$, one millisecond later to $29.4mV + 3.462mV$ and another millisecond later to $E_{L} + 15mV + 6.0984mV$ .
This so called 'resolution trick' is to simulate the spike upswing for $2ms$ after a spike is emitted.
Beside this, the equations of the neuron model and the values of the parameters are equal to the description in the original paper and not presented here.

Here, only a short description about the learning dynamic is given. For further information read the original publication by @Clopath2010.
The discussed learning rule consists of a term for long term depression (LTD) (Eq. @eq:LTD) and long term potentiation (LTP) (Eq. @eq:LTP).

$$ LTP_{Term} = A_{LTP} \bar{x}_i (u - \theta_+)_+ (\bar{u}_+ - \theta_-)_+ $$ {#eq:LTP}

LTP occurs if the presynaptic spike trace ( $\bar{x}_i$) is above zero, the membrane potential $u$ over the threshold $\theta_+$ and the
membrane potential trace $\bar{u}_-$ is above $\theta_-$.
This happens if the postsynaptic neuron spikes short after the presynaptic neuron or if the membrane potential is long enough high,
that $\bar{u}_-$ exceed $\theta_-$.

$$ LTD_{Term} = A_{LTD} (\frac{\bar{\bar{u}}}{u_{ref}^2}) X_i (\bar{u}_{-} - \theta_{-})_+ $$ {#eq:LTD}

If the presynaptic neuron spikes, the spike counter ($X_i$) is set to one, otherwise it is zero.
Further, if the second trace of the postsynaptic membrane potential above $\theta_{-}$ LTD occurs.
This can happen when the presynaptic neurons spikes after the postsynaptic.
The strength of the LTD term, and with that the balance between LTP and LTD, is adjusted over the ratio between $\bar{\bar{u}}$ and the reference value ($u_{ref}^2$),
what implements a homeostatic mechanism.

$$ \tau_{\bar{\bar{u}}}\frac{d \bar{\bar{u}}}{dt} =  [(u-E_L)^2] - \bar{\bar{u}}$$ {#eq:homeo}

// Membrane potential plot and weight distribution plot

These mechanism is computed over the quadratic distance of the membrane potential and the resting potential $E_L$ (Eq. @eq:homeo).
Further, with a higher activity increases the $\bar{\bar{u}}$ and a higher amount of LTD occur and the weights decreases.
In contrast, a lower activity decreases the amount of LTD and the weights can increases.
Over the ratio with $u_{ref}^2$, this mechanism can enforce the connections to decrease down to the minimum weight bound or increase to the maximum
weight bound.
This make a hard upper and lower bound for the weights necessary and leads to a binomial distribution of the weights.

$$ \frac{dw_i}{dt} = -LTD_{Term} + LTP_{Term} $$ {#eq:STDP}

All parameters of the neuron model and the basis set of parameters for the learning rule are taken from the original publication.
Unfortunately, some parameters of the learning rule differs from experiment to experiment. Mainly the value of the homeostatic mechanism and the
maximum weight is different. A table with the different parameters for the different task is .


In the original publication, the authors reproduce some experimental observations,
made in the visual cortex of rats [@Sjoestroem2001].
Further, they investigated the emerged structure of the connectivity depending
on the input activity.


The methods section should explain how you replicated the original results:

* did you use paper description
* did you contact authors ?
* did you use original sources ?
* did you modify some parts ?
* etc.


# Results

Results should be compared with original results and you have to explain why
you think they are the same or why they may differ (qualitative result vs
quantitative result). Note that it is not necessary to redo all the original
analysis of the results.


# Conclusion

Conclusion, at the very minimum, should indicate very clearly if you were able
to replicate original results. If it was not possible but you found the reason
why (error in the original results), you should explain it.


Heading 1                          Heading 2
---------- ----------- ----------- ----------- ----------- -----------
cell1 row1 cell2 row 1 cell3 row 1 cell4 row 1 cell5 row 1 cell6 row 1
cell1 row2 cell2 row 2 cell3 row 2 cell4 row 2 cell5 row 2 cell6 row 2
cell1 row3 cell2 row 3 cell3 row 3 cell4 row 3 cell5 row 3 cell6 row 3
---------- ----------- ----------- ----------- ----------- -----------

Table: Table caption {#tbl:table}

A reference to table @tbl:table.
A reference to figure @fig:logo.
A reference to equation @eq:1.

![Figure caption for Free!](rescience-logo.pdf){#fig:logo}

$$ A = \sqrt{\frac{B}{C}} $$ {#eq:1}


# References
