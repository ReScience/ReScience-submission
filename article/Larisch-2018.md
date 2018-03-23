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
After a spike is the membrane potential
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

## Reimplemented Analysis

In the original publication, the authors mentioned the possibility to reproduce spike timing triplet experiments,
made in the visual cortex of rats [@Sjoestroem2001].
Further, they investigated the emerged structure of the connectivity depending
on the spiking behavior.

To validate the reimplementation, we reproduce the classical spike timing-dependent
learning window (Fig. 2a in [@Clopath2010]) and the frequency repetition task
to reproduce a triplet experiment (Fig.2b in [@Clopath2010]).
Further the influence of spiking order to connectivity
(Fig.4a,down and Fig.4b,down in [@Clopath2010]).
In the available matlab source code, they presenting stable learning of weights
in rate code by 500 input neurons. The firing rate of these neurons follow
a Gaussian distribution and the spike timing a Poisson process
(similar to Fig.5a in [@Clopath2010]). This task is reimplemented as well.
With this analysis, the functionality of the reimplementation is shown on the
main feature of this learning rule, reproduce pair based and triplet STDP data,
and the analysis of connection patterns, depending on the neuronal spiking behavior.

The experiment protocols based on the description on the publication of @Clopath2010.
The implementation of the learning rule was mainly orientated on the
available matlab source code.
Despite the effort to be so close as possible on
the original implementation and description, the internal processing of the equations by ANNarchy
can lead to a different execution order and with this, to other results.
Further, the chosen integration time step can have an influence of the computation result.
On all reimplemented analysis, a time step of $ dt= 1ms$ is chosen.
Because of this, adaption on some parameters was
necessary to reproduce the original results.
For the connectivity experiments, the homeostatic mechanism follow Eq. @eq:homeo.

To reproduce the STDP learning window, we create a list of discrete time points,
where the pre- or postsynaptic neuron is active. The presynaptic neuron spikes
every $50 ms$. The postsynaptic neurons spikes in a range from $1 ms$ to $15 ms$ before
or after the presynaptic neuron.
For the repetition frequency experiment, or the triplet experiment,
the number of pre- and postsynaptic spike pairs increases from a pair frequency of
$0.1 Hz$ to $50 Hz$. The time between a pre- and postsynaptic spike of a pair is
$10 ms$. For this experiments, it was necessary to hold $\bar{\bar{u}}$ fix.
The parameter changes in **Tab.** @tbl:table_FH.

To evaluate the connectivity over the number of spikes, a small network with ten neurons
connected to each other was build up.
Every neuron receives input from one additional neuron, with Poisson distributed
spike timing. Therefore, the PoissonPopulation population type of ANNarchy is used.
The firing rate of every Poisson neuron is increased from two to twenty and with that,
the firing rate of the ten corresponding neurons in the network.

For temporal order of spiking, the Poisson neurons was replaced with the SpikeSourceArray population type from ANNarchy.
With a given list of time steps it is possible to determine the exact spiking time point of the corresponding neuron.

Because the reimplementation of the model based mainly on the matlab source code from  modelDB,
the emergent of stable weights by presenting a Gaussian input over 500 presynaptic neurons and one postsynaptic neuron.
The presynaptic population is implemented with the PoissonPopulation from ANNarchy to achieve a Poisson distributed
spiking behavior, as mentioned in the original matlab implementation.
Equal to the matlab source code, the learning rates ($A_{LTP}$ and $A_{LTD}$) are
ten times faster to speed up the learning.

The experiments for stable weight learning or to show a specific connectivity pattern
require the original homeostatic mechanism.
Changes to the default parameters are shown in **Tab.** @tbl:table_VH.


Task                            Parameter Value             
------------------------------- --------- -------------
Rate based connectivity         $w_{max}$ $0.25 nA$
Temporal based connectivity     $w_{max}$ $0.30 nA$
Stable weight by Gaussian input $w_{max}$ $3.0 nA$
Stable weight by Gaussian input $A_{LTD}$ $1.4*10^{-3}$
Stable weight by Gaussian input $A_{LTP}$ $0.8*10^{-3}$
------------------------------- --------- --------------
Table: Changed parameter for connectivity experiments. {#tbl:table_VH}


Task                 Parameter       Value             
-------------------- --------------- ----------
STDP learning window $\bar{\bar{u}}$ $80 mV^2$
STDP learning window $\theta_{-}$    $-60.0 mV$
triplet experiment   $\bar{\bar{u}}$ $120 mV^2$
-------------------- --------------- ----------
Table: Changed parameter for weight change experiments. {#tbl:table_FH}

\pagebreak

# Results

\begin{figure}
\includegraphics[width=0.5\textwidth]{./figures/deltaW.png}
\includegraphics[width=0.5\textwidth]{./figures/pairing.png}
\caption{ \textbf{Reproduce of experimental findings.}
         \textbf{Left}, the classic STDP learning window. X-axis is the time of a postsynaptic spike in relation to the presynaptic spike.
         \textbf{Right}, weight changes as a function of pair frequency repetition.
         Therefore, pre-post pairs are the blue line and post-pre pairs the red line.}
\label{Fig_exp}
\end{figure}

The \textbf{Fig. \ref{Fig_exp} left} shows the classic pair based spike timing learning window.
If the postsynaptic neuron spikes before the presynaptic one, LTD occurs (red line).
If the postsynaptic neuron spikes after the presynaptic one, LTP occurs (blue line).
Therefore, the x-axis is the time difference between a pre- and postsynaptic spike, with reference to the postsynaptic spike.
The resulting graph is very similar to the presented one in the original publication.
A small different can be seen in the highest positive and negative change.
This can be caused by a different internal processing of ANNarchy.

The analysis of the pairing repetition frequency task is seen in \textbf{Fig. \ref{Fig_exp} right}.
On lower repetition frequency, post-pre pairs (red line) lead to LTD. At a repetition frequency around $30 Hz$,
the post-pre pairs are under the influence of the next post-pre pair and the post-pre-post triplets
lead to LTP. If the repetition frequency of post-pre pairs around $50 Hz$, the same amount of LTP emerges as
in pre-post pairs. Same results was shown in the original paper.

\begin{figure}
\includegraphics[width=0.5\textwidth]{./figures/rate_Code_Weights.png}
\includegraphics[width=0.5\textwidth]{./figures/temporal_Code.png}
\caption{ \textbf{Different connectivity patterns.}
         Depending on the activity emerge different connectivity patterns between the neurons.
         The color scheme is different from them in the original publication.
         Weak connections are yellow, strong unidirectional connections are light green
         and dark green are strong bidirectional connections.
         \textbf{Left}, between neurons with similar high firing rates develop strong connections.
         \textbf{Right}, connection pattern follow the temporal order of the occurred spikes.}
\label{Fig_con}
\end{figure}

Besides the replication of experimental findings of pair based and triplet STDP experiments,
@Clopath2010 presented how the connectivity, they emerge by the proposed learning rule,
is influenced by the spiking behavior of the neurons.
Therefore, \textbf{Fig. \ref{Fig_con} left} shows the connection structure, if the neuron firing with different frequency.
Here, the color scheme is different to them in the original publication.
Weak connection ( above $\frac{3}{4}$ of the maximum activity) are yellow.
Light green are strong unidirectional and dark green are strong bidirectional connections.
Neurons with similar high firing rates develop strong bidirectional connections, because they are often active at the same time point.
This suggest, that the learning is based on correlation between the neuronal activities.
This corresponds with the connection pattern in the original paper.

If the neurons firing in a specific temporal order, this is reflected in the connection pattern ( \textbf{Fig. \ref{Fig_con} left}).
The connections are weak to neurons they fired a long time after or before the neuron.
Meanwhile, the connections are strong to neurons they fired a short time after the neurons.
This is similar to the analysis in the original paper.

\begin{figure}
\includegraphics[width=1\textwidth]{./figures/weights_stable.png}
\caption{\textbf{Stable weights on Poisson distributed Input.} }
\label{Fig_stab}
\end{figure}
% was ist eine Epoche ?!

The \textbf{Fig. \ref{Fig_stab}} shows the development of the weights over time, by presenting
Gaussian shaped input. Around 500 epochs, stable weights begin to emerge.
Similar can be observed in running the original matlab source code.

# Conclusion

The here proposed reimplementation of voltage based STDP learning rule from @Clopath2010,
is be able to reproduce the experimental data and the emergent connectivity structures
as proposed in the original paper.
The description of the learning rule in the original publication is detailed to understand
the different components of the learning rule.
Certainly, two main components are described inadequate.
The 'resolution trick' of the membrane potential and the equation for the
homeostatic mechanism relevant $\bar{\bar{u}}$.
Stable learning is only accessible with a good working homeostasis mechanism,
what regulated the LTD amount. The dependency of $\bar{\bar{u}}$ make it necessary to
implement the right behavior of the membrane potential, what needs the 'resolution trick'.
So the reimplementation benefits from the release of the source code on modelDB.
Nonetheless, the reimplementation with the ANNarchy frame work was successful. 
# References
