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

Since the first description of spike timing-dependent plasticity (STDP) [@Bi1998],
different models of STDP have been published to reproduce different experimental findings.
Early implementations such as pair-based STDP learning rules failed to reproduce some experimental observations,
for example triplet or quadruplets experiments [@Pfister2006].

@Clopath2010 introduced a STDP model able to reproduce the experimental findings of triplet studies.
They propose a biologically-motivated model with a voltage-based learning rule, where the occurrence of long term depression (LTD) or long term potentiation (LTP) depends on the postsynaptic membrane voltage.
Clopath and colleagues could reproduce different STDP experiments. **TODO: summarize which ones**
Furthermore, they were able to show that their learning rule
can develop stable weights, as needed for learning the receptive fields of V1 simple cells.
Therefore, they implemented a homeostatic mechanism to
control the amount of generated LTD, based on the relationship between the average postsynaptic
membrane potential and a reference value.
Their model led to two different connectivity structures,
depending on the spiking behavior of the neurons: if the neurons fire trongly at the same time, they build strong bidirectional connections (as in rate-coded Hebbian learning). If they fire in a specific temporal order, their connectivity structure follows that order (temporal coding).

# Methods

## Overview

The original model was implemented in Matlab (<http://modeldb.yale.edu/144566>) to demonstrate the
stable learning of weights.
This model reimplementation is written in Python (v2.7) with the help of the neuro-simulator ANNarchy [@Vitay2015], numpy (v1.11.0) and matplotlib (v1.5.1).
The reimplementation is mainly based on the description of neuron model and learning rule in the original publication [@Clopath2010].
Because of the lack of further description of the homeostatic mechanism and the neural behavior after a emitted spike,
the Matlab code is used as the second reference for this reimplementation.

## Model description

### Neural model {-}

@Clopath2010 used an adaptive exponential integrate-and-fire (AdEx) neuron in their model. The exact neural model is mainly derived from the description in the Matlab source code (see **Eq.** @eq:memb):

$$ C\frac{du}{dt} = -g_L\, (u-E_L) + g_L \, \Delta_T \, e^{\frac{u-V_T}{\Delta_T}} - w_{ad} + z + I$$ {#eq:memb}

where $u$ is the membrane potential, $C$ the membrane capacitance, the leak conductance is $g_L$ and $E_L$ is the resting potential.
The slope factor ($\Delta_T$) and the spiking threshold ($V_T$) are describing the behavior of the exponential term.

If the membrane potential ($u$) is above $V_T$, the neuron spikes and the membrane potential increases exponentially.
To simulate the spike upswing for $2ms$ after a spike is emitted, they used the so called 'resolution trick' **TODO: what is it?**.  
The membrane potential is set to $29.4mV$,  one millisecond later to $29.4mV + 3.462mV$ and another millisecond later to $E_{L} + 15mV + 6.0984mV$ .

The depolarizing spike afterpotential is $z$ and the hyperpolarization current described by $w_{ad}$.
Together with the adaptive spiking threshold ($V_T$), they are changing over the time.
Their descriptions are equal to those in the original paper. **TODO: give the equation**
The values of the parameters of the neural model are taken from the original paper and are presented in the.

### Synaptic model {-}

Here, only a short description about the learning dynamic is given. For further information read the original publication by @Clopath2010.

**No, the article should be standalone, reproduce the equations**

The proposed learning rule consists of two terms: a long term potentiation (LTP) term (**Eq.** @eq:LTP) controls increases in synaptic efficiency:

$$ LTP_{Term} = A_{LTP} \, \bar{x}_i \, (u - \theta_+)_+ \, (\bar{u}_+ - \theta_-)_+ $$ {#eq:LTP}

With this term, LTP occurs when the presynaptic spike trace ($\bar{x}_i$ **TODO: how is it computed?**) is above zero, the membrane potential $u$ is over the threshold $\theta_+$ and the membrane potential trace $\bar{u}_-$ is above $\theta_-$. This happens when the postsynaptic neuron spikes shortly after the presynaptic neuron or if the membrane potential is high long enough, so that $\bar{u}_-$ exceeds $\theta_-$.

The long term depression (LTD) term (**Eq.** @eq:LTD) governs decreases in synaptic efficiency:

$$ LTD_{Term} = A_{LTD} \, (\frac{\bar{\bar{u}}}{u_{ref}^2}) \, X_i \, (\bar{u}_{-} - \theta_{-})_+ $$ {#eq:LTD}

The presynaptic spike counter ($X_i$) is set to one after a spike and  zero otherwise.
$\bar{u}_{-}$ is a second trace of the postsynaptic membrane potential.
If it exceeds the $\theta_{-}$ threshold, and a presynaptic spike is emitted, LTD occurs. This happens when the presynaptic neurons spikes after the postsynaptic one.

The amplitude of the LTD term, and with that the balance between LTP and LTD, is adjusted with respect to the ratio between $\bar{\bar{u}}$ and a reference value ($u_{ref}^2$), hence implementing a homeostatic mechanism (**Eq.** @eq:homeo).

$$ \tau_{\bar{\bar{u}}}\frac{d \bar{\bar{u}}}{dt} =  [(u-E_L)^2] - \bar{\bar{u}}$$ {#eq:homeo}

This mechanism is computed over the quadratic distance of the membrane potential and the resting potential $E_L$. **Sentence not very clear**
With a higher activity, $\bar{\bar{u}}$ increases, leading to a higher amount of LTD and a weight decrease.
In contrast, a lower activity decreases the amount of LTD and the weights can increase.
Through the ratio of $\bar{\bar{u}}$ with $u_{ref}^2$, this mechanism can enforce the connections to decrease down to the minimum weight bound or increase to the maximum weight bound.
This requires a hard upper and lower bound for the weights and leads to a binomial distribution of the weights.
The weight change over time depends on the positive LTP and the negative LTD terms (**Eq.** @eq:STDP):

$$ \frac{dw}{dt} = -LTD_{Term} + LTP_{Term} $$ {#eq:STDP}

All parameters of the neuron model and the basis set of parameters for the learning rule are taken from the original publication.
Some parameters of the learning rule differ from experiment to experiment, in particular the value of the homeostatic mechanism (**be more precise**) and the
maximum weight value.
A table with the different parameters for each task is presented in **Tab.** @tbl:table_FH and **Tab.** @tbl:table_VH.

## Reproduction of experiments

In the original publication, the authors reproduce spike timing triplet experiments in the visual cortex of rats [@Sjoestroem2001].
Furthermore, they investigate the emerged structure of the connectivity depending on the spiking behavior.

To validate the reimplementation, we reproduce the classical spike timing-dependent
learning window (**Fig.** 2a in [@Clopath2010]), the frequency repetition task
to reproduce a triplet experiment (**Fig.** 2b in [@Clopath2010]) and the influence of spiking order to connectivity
(**Fig.** 4a, down and **Fig.** 4b, down in [@Clopath2010]).
In the available matlab source code, they investigate the stable learning of weights
using 500 input neurons and one postsynaptic neuron.
The firing rate of these neurons follow
a Gaussian distribution and the spike timing a Poisson process
(similar to **Fig.** 5a in [@Clopath2010]). This task is reimplemented as well.
With this analysis, the functionality of the reimplementation is shown on the
main feature of this learning rule.
We reproduced pair-based and triplet STDP data and analyzed the obtained connection patterns depending on the neuronal spiking behavior.

**TODO: are there experiments you did not reproduce and why? It is important for the partial/full replication**

The experiment protocols are based on the description on the publication of @Clopath2010.
The implementation of the learning rule was mainly made according to the
available matlab source code.
Despite the effort to be as close as possible to
the original implementation and description, the internal processing of the equations by ANNarchy
can lead to a different execution order. Therefore, differing results occurred.

**Be more specific: what is different? why?**

Further, the chosen integration time step can have an influence of the computation result.
On all reimplemented analysis, a time step of $dt= 1ms$ is chosen. **Why? What is it in the original paper?**
Because of that, some parameters had to be adapted to reproduce the original results. **Which ones?**

To reproduce the STDP learning window **Which Fig?**, we create a list of discrete time points where the pre- or postsynaptic neurons should emit spikes. The presynaptic neuron spikes every $50 ms$. The postsynaptic neurons spikes in a range from $1 ms$ to $15 ms$ before
or after the presynaptic neuron.
For the repetition frequency experiment or the triplet experiment **Figs?**,
the number of pre- and postsynaptic spike pairs increases from a pair frequency of
$0.1 Hz$ to $50 Hz$. The time between a pre- and postsynaptic spike of a pair is
$10 ms$. To reproduce this experiments, it was necessary to set $\bar{\bar{u}}$ to a fixed value. **Is it different from the paper?**
The parameter changes are shown in **Tab.** @tbl:table_FH.

To analyze the connectivity depending the number of spikes, a small network with ten neurons
connected with each other is build.
Every neuron receives input from one additional neuron, with Poisson-distributed
spike patterns.
The firing rate of each Poisson neuron is increased from 2Hz to 20Hz, influencing the firing rate of the 10 corresponding neurons in the network.


Because the reimplementation of the model is mainly based on the Matlab source code from modelDB,
the emergence of stable weights was achieved by presenting a Gaussian input over 500 presynaptic neurons and one postsynaptic neuron (**contrary to what in the original paper?**.
For every trial ($125$ms) ten Gaussian patterns are created to determine the activity of the 500 input neurons.

As in the Matlab source code, the learning rates ($A_{LTP}$ and $A_{LTD}$) are
increased by a factor ten to speed up the learning.

The experiments for stable weight learning, to show a specific connectivity pattern, require the original homeostatic mechanism.
Changes to the default parameters for this tasks are shown in **Tab.** @tbl:table_VH.


Task                            Parameter Value             
------------------------------- --------- -------------
Rate based connectivity         $w_{max}$ $0.25 nA$
Temporal based connectivity     $w_{max}$ $0.30 nA$
Stable weight by Gaussian input $w_{max}$ $3.0 nA$
Stable weight by Gaussian input $A_{LTD}$ $1.4*10^{-3}$
Stable weight by Gaussian input $A_{LTP}$ $0.8*10^{-3}$
------------------------------- --------- --------------
Table: Changed parameters for connectivity experiments. {#tbl:table_VH}


Task                 Parameter       Value             
-------------------- --------------- ----------
STDP learning window $\bar{\bar{u}}$ $80 mV^2$
STDP learning window $\theta_{-}$    $-60.0 mV$
triplet experiment   $\bar{\bar{u}}$ $120 mV^2$
-------------------- --------------- ----------
Table: Changed parameters for weight change experiments. {#tbl:table_FH}

## Reimplementation

**What is missing is description of the reimplementation itself, i.e. present ANNarchy, which structures are used (I suppressed the references to PoissonPopulation and co in the previous part), what was hard, why ANNarchy and not pure Python?, etc**

**You can put code in the article:**

```python
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
```

# Results

**Perhaps structure more the results**

## Pair-based and triplet STDP experiments

\begin{figure}
\includegraphics[width=0.5\textwidth]{./figures/deltaW.png}
\includegraphics[width=0.5\textwidth]{./figures/pairing.png}
\caption{ \textbf{Replication of experimental findings.}
         \textbf{Left}, the classic STDP learning window. On the x-axis is the time of a postsynaptic spike in relation to the presynaptic spike presented.
         \textbf{Right}, weight changes as a function of pair frequency repetition.
         Pre-post pairs are the blue line and post-pre pairs the red line.}
\label{Fig_exp}
\end{figure}

The classic-pair based spike timing learning window is presented in \textbf{Fig. \ref{Fig_exp} left}.
If the postsynaptic neuron spikes before the presynaptic one, LTD occurs (red line).
If the postsynaptic neuron spikes after the presynaptic one, LTP occurs (blue line).
The x-axis represents the time difference between pre- and post-synaptic spikes, relative to the postsynaptic spike.
The resulting graph is similar to the presented one in the original publication.
A small difference can be seen in the higher positive and negative change (**From how much?**).
This could be caused by a different internal processing of ANNarchy. (**Be more precise...**)

The analysis of the pairing repetition frequency task is shown on \textbf{Fig. \ref{Fig_exp} right}.
With lower repetition frequency, post-pre pairs (red line) lead to LTD. At a repetition frequency around $30 Hz$,
the post-pre pairs are under the influence of the next post-pre pair and the post-pre-post triplets
lead to LTP. If the repetition frequency of post-pre pairs is around $50 Hz$, the same amount of LTP emerges as
in pre-post pairs. These results are similar to the original paper.

\begin{figure}
\includegraphics[width=0.5\textwidth]{./figures/rate_Code_Weights.png}
\includegraphics[width=0.5\textwidth]{./figures/temporal_Code.png}
\caption{ \textbf{Different connectivity patterns.}
         Depending on the spiking activity, different connectivity patterns emerge between the neurons.
         The color scheme is different from the original publication.
         Weak connections are yellow, strong unidirectional connections are light green
         and dark green are strong bidirectional connections. **Why not use the same color code, then?**
         \textbf{Left}, Neurons with similar high firing rates develop strong bidirectional connections.
         \textbf{Right}, connection pattern follows the temporal order of the occurred spikes.}
\label{Fig_con}
\end{figure}

## Connectivity analysis

In addition to the replication of experimental findings of pair-based and triplet STDP experiments,
@Clopath2010 presented how the synaptic connectivity, emerging from the proposed learning rule,
is influenced by the spiking behavior of the neurons.
\textbf{Fig. \ref{Fig_con} left} shows the obtained connection structure if neurons fire with different frequencies.
Here, the color scheme is different from the original publication (**why?**).
Weak connections ( above $\frac{3}{4}$ of the maximum activity) are yellow.
Light green represents strongly unidirectional connections while dark green represents strong bidirectional connections.
Neurons with similarly high firing rates develop strong bidirectional connections, because they are often active at the same time.
This suggests that learning is based on the correlation between the neuronal activities in a Hebbian manner.
This is in line with the connection pattern presented in the original paper. (**be more quantitative there**)

If the neurons fire in a specific temporal order, this sequence is reflected in the connection pattern (see \textbf{Fig. \ref{Fig_con} left}).
As in the original paper, the connections from neurons which are firing a long time after or before the post-synaptic neuron are weak, while they are strong to neurons which fired a short time after the neurons.

\begin{figure}
\includegraphics[width=1\textwidth]{./figures/weights_stable.png}
\caption{\textbf{Stable weights on Poisson distributed Input.}
    Colors show the weight value at the end of the current epoch from the presynaptic neuron to a single postsynaptic neuron.
    On the y-Axis is the presynaptic index indicated and the x-axis shows the number of epoch.}
\label{Fig_stab}
\end{figure}

The development of the weights over time, by presenting
Gaussian shaped inputs, is presented in \textbf{Fig. \ref{Fig_stab}},
After 500 epochs, stable weights begin to emerge, similar to the original matlab source code. **More details**

# Conclusion

Our reimplementation of voltage based STDP learning rule from @Clopath2010
is able to reproduce the experimental data and the emergent connectivity structures
proposed in the original paper.

The description of the learning rule in the original publication comprises enough details
to understand the different components and their interaction.
However, two main components have not been described adequately to allow a direct reimplementation: the 'resolution trick' of the membrane potential after spike emission and the equation for the
homeostatic mechanism ($\bar{\bar{u}}$).  
Stable learning is only possible with a working homeostasis mechanism, that regulates the amount of LTD. **Was it described in the paper, or did you have to guess it from the Matlab code?**
The dependency of $\bar{\bar{u}}$ on the homeostatic mechanism is necessary to
implement the right behavior of the membrane potential.
The reimplementation has greatly benefited from the release of the source code on modelDB.

## Acknowledgment

This work was supported by the European Social Fund (ESF) and the Freistaat Sachsen.


# References
