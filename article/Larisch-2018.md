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
Clopath and colleagues could reproduce how the occurrence of LTD and LTP depends
on the depolarizing of the postsynaptic membrane potential, as observed in voltage-clamp [@Ngezahayo2000] and stationary-depolarization experiments [@Artola1990].  
Further, they could reproduce experimental finding from spike pair repetition and triplet experiments [@Sjoestroem2001].
They were able to show that their learning rule
can develop stable weights, as needed for learning the receptive fields of V1 simple cells.
Therefore, they implemented a homeostatic mechanism to
control the amount of generated LTD, based on the relationship between the average postsynaptic
membrane potential and a reference value.
Their model led to two different connectivity structures,
depending on the spiking behavior of the neurons: if the neurons fire strongly at the same time, they build strong bidirectional connections (as in rate-coded Hebbian learning). If they fire in a specific temporal order, their connectivity structure follows that order (temporal coding).

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
To simulate the spike upswing for $2ms$ after a spike is emitted, they used the so called 'resolution trick'.
For this, they simulated once the complete change in the membrane potential through a spike with high precision and integrated.
For simulations, they used the integrated number and clamped the membrane potential for $2ms$.  
This means, the membrane potential is set to $29.4mV$ after a spike,  one millisecond later to $29.4mV + 3.462mV$ and another millisecond later to $E_{L} + 15mV + 6.0984mV$ .

The depolarizing spike afterpotential is $z$ and it decays over time to zero (see **Eq.** @eq:Z).
$$ \tau_{z} \frac{dz}{dt} =-z  $$ {#eq:Z}
The hyperpolarization current described by $w_{ad}$ (see **Eq.** @eq:wad).
After a spike, $w_{ad}$ is increased by the amount $b$ and decays down to the resting potential $E_{L}$ otherwise.
$$ \tau_{w_{ad}} \frac{dw_{ad}}{dt} = a(u-E_{L})-w_{ad}  $$ {#eq:wad}
The adaptive spiking threshold ($V_T$) is set to $V_{T_{max}}$ after the neuron spiked and will decay to $V_{T_{rest}}$ (see **Eq.** @eq:VT).
$$ \tau_{V_{T}} \frac{dV_{T}}{dt} =- (V_T - V_{T_{rest}})  $$ {#eq:VT}
The values of the parameters of the neural model are taken from the original paper.

### Synaptic model {-}

The proposed learning rule consists of two terms: a long term potentiation (LTP) term (**Eq.** @eq:LTP) controls increases in synaptic efficiency:

$$ LTP_{Term} = A_{LTP} \, \bar{x}_i \, (u - \theta_+)_+ \, (\bar{u}_+ - \theta_-)_+ $$ {#eq:LTP}

Therefore, $A_{LTP}$ is the learning rate for the LTP term.
The parameters $\theta_{+}$ and $\theta_{-}$ are plasticity thresholds for the membrane potential ($u$),
respectively for the over time averaged version ($\bar{u}_+$ (see **Eq.** @eq:barU)).
The original paper they did not mentioned the meaning of both thresholds.
With $\theta_{+} = -45.3mV$ it is above the spiking threshold.
This suggest, that the threshold avoid the occurrence of LTP if the postsynaptic neuron not spiked already.
With $\theta_{+} = E_{L}$, LTP only occurs if the membrane potential is above the resting potential.

$$ \tau_{+} \frac{d\bar{u}_+}{dt} = -\bar{u}_{+} + u $$ {#eq:barU}

On every spike of the presynaptic neuron, the spike train $\bar{x}_i$ is increased by one and will decay over time with a time constant $\tau_{x}$ (see **Eq.** @eq:xbar).
If the presynaptic neuron spikes to a time point $t$, the spike counter $X_i$ is $1$, otherwise $0$.

$$ \tau_{x}\frac{d \bar{x}_i}{dt} = -\bar{x}_i + X_i $$ {#eq:xbar}

With this term, LTP occurs when the presynaptic spike trace ($\bar{x}_i$) is above zero, the postsynaptic membrane potential $u$ is over the threshold $\theta_+$ and the membrane potential trace $\bar{u}_-$ is above $\theta_-$. This happens when the postsynaptic neuron spikes shortly after the presynaptic neuron or if the membrane potential is high long enough, so that $\bar{u}_-$ exceeds $\theta_-$.


The second term is the long term depression (LTD) term (**Eq.** @eq:LTD) and governs decreases in synaptic efficiency:

$$ LTD_{Term} = A_{LTD} \, (\frac{\bar{\bar{u}}}{u_{ref}^2}) \, X_i \, (\bar{u}_{-} - \theta_{-})_+ $$ {#eq:LTD}

The presynaptic spike counter ($X_i$) is set to one after a spike and  zero otherwise (as mentioned above).
$\bar{u}_{-}$ is a second trace of the postsynaptic membrane potential similar to $\bar{u}_{+}$, but with $\tau_{-} > \tau_{+}$.
If it exceeds the $\theta_{-}$ threshold, and a presynaptic spike is emitted, LTD occurs. This happens when the presynaptic neurons spikes after the postsynaptic one.

The amplitude of the LTD term, and with that the balance between LTP and LTD, is adjusted with respect to the ratio between $\bar{\bar{u}}$ and a reference value ($u_{ref}^2$), hence implementing a homeostatic mechanism (**Eq.** @eq:homeo).

$$ \tau_{\bar{\bar{u}}}\frac{d \bar{\bar{u}}}{dt} =  [(u-E_L)^2] - \bar{\bar{u}}$$ {#eq:homeo}

Therefore, the homeostatic variable $\bar{\bar{u}}$ is computed over the quadratic difference between the postsynaptic membrane potential and the resting potential ($E_L$).
With a higher activity, $\bar{\bar{u}}$ increases, leading to a higher amount of LTD and a weight decrease.
In contrast, a lower activity decreases the amount of LTD and the weights can increase.
Through the ratio of $\bar{\bar{u}}$ with $u_{ref}^2$, this mechanism can enforce the connections to decrease down to the minimum weight bound or increase to the maximum weight bound.
This requires a hard upper and lower bound for the weights and leads to a binomial distribution of the weights.
The weight change over time depends on the positive LTP and the negative LTD terms (**Eq.** @eq:STDP):

$$ \frac{dw}{dt} = -LTD_{Term} + LTP_{Term} $$ {#eq:STDP}

All parameters of the neuron model and the basis set of parameters for the learning rule are taken from the original publication.
Some parameters of the learning rule differ from experiment to experiment, in particular the reference value of the homeostatic mechanism ($u_{ref}^2$),
the learning rates for the LTP and LTD terms ($A_{LTP}$ and $A_{LTD}$), $\theta_{-}$ and the maximum weight value ($w_{max}$).
A table with the different parameters for each task is presented in **Tab.** @tbl:table_FH and **Tab.** @tbl:table_VH.

## Reproduction of experiments

In the original publication, the authors reproduce spike timing triplet experiments in the visual cortex of rats [@Sjoestroem2001].
Furthermore, they investigate the emerged structure of the connectivity depending on the spiking behavior.

To validate the reimplementation, we reproduce the voltage clamp experiment (**Fig.** 1h in [@Clopath2010]), the classical spike timing-dependent
learning window (**Fig.** 2a in [@Clopath2010]), the frequency repetition task
to reproduce a triplet experiment (**Fig.** 2b in [@Clopath2010]), the influence of spiking order to connectivity
(**Fig.** 4a, down and **Fig.** 4b, down in [@Clopath2010]) and the emergent of receptive fields by presenting natural scenes.
In the available matlab source code, they investigate the stable learning of weights
using 500 input neurons and one postsynaptic neuron.
The firing rate of these neurons follow
a Gaussian distribution and the spike timing a Poisson process.
This task is reimplemented as well.
With this analysis, the functionality of the reimplementation is shown on the
main feature of this learning rule.

**TODO: are there experiments you did not reproduce and why? It is important for the partial/full replication**

The experiment protocols are based on the description on the publication of @Clopath2010.
The implementation of the learning rule was mainly made according to the
available Matlab source code.
Despite the effort to be as close as possible to
the original implementation and description, the internal processing of the equations by ANNarchy
can lead to a different execution order. Therefore, differing results occurred.

**Be more specific: what is different? why?**

Further, the chosen integration time step can have an influence of the computation result.
In the original publication, no integration time step is mentioned. In the published Matlab source code is
a time step of $dt= 1ms$ chosen, which we chose too.

To reproduce the STDP learning window (see \textbf{Fig. \ref{Fig_exp}} left), we create a list of discrete time points where the pre- or postsynaptic neurons should emit spikes. The presynaptic neuron spikes every $50 ms$. The postsynaptic neurons spikes in a range from $1 ms$ to $15 ms$ before
or after the presynaptic neuron.
For the repetition frequency experiment or the triplet experiment (see \textbf{Fig. \ref{Fig_exp}} right),
the number of pre- and postsynaptic spike pairs increases from a pair frequency of
$0.1 Hz$ to $50 Hz$. The time between a pre- and postsynaptic spike of a pair is
$10 ms$. To reproduce this experiments, it was necessary to set $\bar{\bar{u}}$ to a fixed value as mentioned in the original publication [@Clopath2010].
The parameter changes are shown in **Tab.** @tbl:table_FH.

To analyze the connectivity depending the number of spikes, a small network with ten neurons
connected with each other is build.
Every neuron receives input from one additional neuron, with Poisson-distributed
spike patterns.
The firing rate of each Poisson neuron is increased from 2Hz to 20Hz, influencing the firing rate of the 10 corresponding neurons in the network.

The reimplementation of the model is mainly based on the Matlab source code from modelDB.
Besides of experiments in the original publication, we reimplemented the experiment for the emergent of stable weights out of the Matlab source code.
The emergence of stable weights was achieved by presenting a Gaussian input over 500 presynaptic neurons and one postsynaptic neuron.
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
The reimplementation was done with Python 2.7 and the neural-simulator ANNarchy [@Vitay2015] (v.4.6.4).
With ANNarchy, it is possible to implement the description of the learning and the neuronal behavior
by define the mathematical equations and it cares about the temporal execution of the equations.
This supports the implementation of more complex models and bigger neuronal networks.
Therefore, ANNarchy supports rate based and spiking learning rules, and it provides a way to combine both kinds of neuronal networks.
The definition of learning rules, the neurons and the network structure is done in the easy understanding Python language.
To archive a good performance on the execution of the model, the Python code is compiled in C++.
With that, ANNarchy make it easy to implement the a neuronal network model and show good performances [@Vitay2015].

As mentioned in the original publication, to reproduce the voltage-clamp experiment, the pairing repetition task and the STDP learning window we set $\bar{bar{u}} = u_{ref}$. The implementation of the network for this tasks can be found in **net_fix.py**.
For the connectivity experiments and for the emergent of stable weights, the homeostatic mechanism dynamic as described
in the original publication [@Clopath2010]. The network with a dynamic homeostatic mechanism can be found in **net_homeostatic.py**.
The following explanation of the network implementation is from the **net_homeostatic.py**.

### Network implementation
``` python
neuron_eqs = """
dvm/dt  = if state>=2:+3.462 else:
          if state==1: -(vm+51.75)+
          1/C*(Isp - (wad+b))+ g_Exc else:
          1/C * ( -gL * (vm - EL) + gL * DeltaT *
          exp((vm - VT) / DeltaT) - wad + z )+
          g_Exc: init = -70.6
dvmean/dt = (pos(vm - EL)**2 - vmean)/taumean    :init = 0.0
dumeanLTD/dt = (vm - umeanLTD)/tauLTD : init=-70.0
dumeanLTP/dt = (vm - umeanLTP)/tauLTP : init =-70.0
dxtrace /dt = (- xtrace )/taux
dwad/dt = if state ==2:0 else:
          if state==1:+b/tauw else:
          (a * (vm - EL) - wad)/tauw : init = 0.0
dz/dt = if state==1: -z+Isp-10 else:
        -z/tauz  : init = 0.0
dVT/dt = if state==1: +(VTMax - VT)-0.4 else:
         (VTrest - VT)/tauVT  : init=-50.4
dg_Exc/dt = -g_Exc/tau_gExc
state = if state > 0: state-1 else:0
Spike = 0.0  """
```

The code above shows the definition of neuron model equations in ANNarchy.
All necessary equations are typed in one string variable ('neuron_eqs').
The variable 'vm' describes the membrane potential $u$, 'vmean' the homeostatic variable $\bar{bar{u}}$,
'umeanLTD' and 'umeanLTP' are the equations for $u_{-}$, respectively $u_{+}$.
The variable 'xtrace' describes $\bar{x}$, 'wad' is $w_{ad}$, 'z' is $z$, 'g_Exc' is the input current and 'Spike' is the spike counter ($X$).
To implement the 'resolution trick', we use a extra discrete variable 'state'.
With that, we control the behavior of the different variables after a spikes to recreate the behavior of the variables as in the Matlab source file.
The neuron spikes only if the membrane potential exceeds the threshold and if the 'state' variable is equal to zero.

``` python
spkNeurV1 = Neuron( parameters = params,
                    equations=neuron_eqs,
                    spike="""(vm>VT) and (state==0)""",
                    spike="""(vm>VT) and (state==0)""",
                    reset="""vm = 29.0
                          state = 2.0
                          VT = VTMax
                          Spike = 1.0
                          xtrace+= 1/taux""")
```

To define a neuron model, ANNArchy provides the __Neuron__ object, what expects a
string object for the parameters ('parameters'), the equations that describes the neuronal behavior ('equations'),
a string that define the conditions to release a spike ('spike') and a string that defines the changes in the variables after a spike ('reset').
With that, we reproduce the behavior of the neuron as described in the published Matlab source code.


``` python
equatSTDP = """
    ltdTerm = if w>wMin : (aLTD*(post.vmean/urefsquare)*
              pre.Spike * pos(post.umeanLTD - thetaLTD)) else : 0.0
    ltpTerm = if w<wMax : (aLTP * pos(post.vm - thetaLTP)*
              (pre.xtrace)* pos(post.umeanLTP - thetaLTD)) else : 0.0
      deltaW = ( -ltdTerm + ltpTerm)
        dw/dt = deltaW :min=0.0"""
```

As for the neuron model, the equations for the spiking learning are defined by strings of the differential equations.
The 'ltdTerm' describes the $LTD_{Term}$ and the 'ltpTerm' the $LTP_{Term}$.
Variables of the pre- or post-synaptic neuron, they are define in the neuron model, can be addressed with the prefix 'pre.', respectively 'post.'.
With the 'if w>wMin' statement in the 'ltdTerm', the weight only decreases if the weight is above the lower bound.
In the 'ltpTerm' a analogous term is implement to avoid, that weights exceeds the upper bound.
The parameters are defined in a string, analogous to the parameters of the neuron model.
Therefore, the parameter 'urefsquare' is the homeostatic reference parameter $u^{2}_{ref}$.
The learning rates $A_{LTD}$ and $A_{LTP}$ are 'aLTD', respectively 'aLTP'.
And the threshold $\theta_{-}$ is defined by 'thetaLTD' and $\theta_{+}$ by 'thetaLTP'.
The parameter 'transmit' is zero or one, dependent if the synaptic current for the experiment should transmit or not.

``` python

ffSyn = Synapse( parameters = parameterFF,
    equations= equatSTDP,
    pre_spike='''g_target += w*transmit''')
```

ANNarchy provides a __Synapse__ object, what expects a 'parameters' argument,
the 'equations' and a description if the pre-synaptic neuron spikes ('pre_spike').
After a pre-synaptic spike, the input current of the post-synaptic neurons increases by the value of the synaptic weight.
Therefore, 'g_target' is the target variable on the post-synaptic side what should be increased.
The target variable is defined in the __Projection__ object (see below).
Additionally, a description for a post-synaptic spike is possible.

### Implementation of the Experiments

The implementation of the different experiments are in the current python files.
To perform an experiment, the network with the neuron populations and the weights between them must be
initialized.
To create a population, ANNarchy provides the __Population__ object.
The 'geometry' argument expects a tuple or a integer and defines the spatial geometry, respectively the number of neurons in the population.
The 'neuron' arguments expects a __Neuron__ object. It defines the used model for population neurons.
The population can be give a unique name, optionally.
Besides that, it exist a set of predefined __Population__ objects in ANNarchy,
for example the __PoissonPopulation__.
This object provides a population of spiking neurons, which spiking behavior follows a Poisson distribution.
As for the __Population__ object with the 'geometry' argument is the size or the spatial geometry of the population.
The argument 'rates' defines the mean firing rate of the population neurons.
A name can be given, optionally.

``` python
poisPop = PoissonPopulation(geometry=10, rates=100.0)
pop_Ten = Population(geometry=10, neuron=spkNeurV1, name="pop_Ten")
```
To connect two neuron populations and define the weight matrix between them, ANNarchy provides the __Projection__ object.
The 'pre' argument define the pre-synaptic population and the 'post' argument the post-synaptic population.
Both arguments expects a __Population__ object.
The 'target' argument defines the target variable of the post-synaptic neuron, which is increased by the weight value after a pre-synaptic spike.

``` python
projInp_Ten = Projection(
    pre = poisPop,
    post= pop_Ten,
    target='Exc'
).connect_one_to_one(weights = 30.0)

projTen_Ten = Projection(
    pre= pop_Ten,
    post= pop_Ten,
    target= 'Exc',
    synapse= ffSyn               
).connect_all_to_all(weights = 0.1,allow_self_connections=True)
```

A further description of the experiment implementations can be found in the corresponding python files.

### Recording variables

With the **Monitor** object provides ANNarchy a easy possibility to record variables from projections and populations.

``` python
    dendrite = projV1_V1.dendrite(0)
    m_d = Monitor(dendrite, ['w','deltaW','ltdTerm','ltpTerm'])
```

# Results

**Perhaps structure more the results**
## Voltage-Clamp experiment

\begin{figure}
\centering
\includegraphics[width=0.5\textwidth]{./figures/W_hippo.png}
\caption{TestCapture}
\label{Fig_hipo}
\end{figure}

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
A small difference can be seen in the higher positive and negative change.
In the original publication is the normalized weight at a time difference of $-10ms$ around $70 \%$.
In our result, the weight is around $80 \%$.
This could have two reasons. First, we use
Second, this could be caused by a different internal processing of ANNarchy. (**Be more precise...**)
In the Matlab source code, the equations are calculated step by step for each integration step.
In ANNArchy, the changes are calculated separately for the neurons and synapses variables, then the changes added, then
it will be checked if a postsynaptic event appeared, and than the recording is happening.

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
         The color scheme is similar to them in the original publication.
         Weak connections are blue, strong unidirectional connections are yellow
         and red are strong bidirectional connections.
         \textbf{Left}, Neurons with similar high firing rates develop strong bidirectional connections.
         \textbf{Right}, connection pattern follows the temporal order of the occurred spikes.}
\label{Fig_con}
\end{figure}

## Connectivity analysis

In addition to the replication of experimental findings of pair-based and triplet STDP experiments,
@Clopath2010 presented how the synaptic connectivity, emerging from the proposed learning rule,
is influenced by the spiking behavior of the neurons.
\textbf{Fig. \ref{Fig_con} left} shows the obtained connection structure if neurons fire with different frequencies.
Here, the color scheme is similar to the original publication.
Weak connections ( above $\frac{3}{4}$ of the maximum activity) are blue.
Yellow represents strongly unidirectional connections while red represents strong bidirectional connections.
Neurons with similarly high firing rates develop strong bidirectional connections, because they are often active at the same time.
This suggests that learning is based on the correlation between the neuronal activities in a Hebbian manner.
Weak connections are emerged to neurons with a low firing under $5 Hz$ rate.
If the postsynaptic neuron firing with a rate above $5 Hz$, strong unidirectional weights emerge.
This is in line with the connection pattern presented in the original paper (see **Fig. 4** in [@Clopath2010]).

If the neurons fire in a specific temporal order, this sequence is reflected in the connection pattern (see \textbf{Fig. \ref{Fig_con} left}).
As in the original paper, the connections from neurons which are firing a long time after or before the post-synaptic neuron are weak, while they are strong to neurons which fired a short time after the neurons.

\begin{figure}
\includegraphics[width=1\textwidth]{./figures/weights_stable.png}
\caption{\textbf{Stable weights on Poisson distributed Input.}
    Colors show the weight value at the end of the current epoch from the presynaptic neuron to a single postsynaptic neuron. Blue are weight values around zero and red are weight values around the maximum weight value of $3$.
    On the y-Axis is the presynaptic index indicated and the x-axis shows the number of epoch.}
\label{Fig_stab}
\end{figure}

In the published Matlab source code, they demonstrated the emergent of stable weights.
Therefore, they presented a one dimensional input with 500 input values.
At every time step, a subset of near to each other lying neurons are active.
The emergent stable weights are shown in \textbf{Fig. \ref{Fig_stab}}.
After 500 epochs, a spatial related subset of weights increased to the maximum and the other weight values decrease down to zero.
This leads to a specific selectivity for the postsynaptic neuron.


**Was war für die Implementierung des Modelles in ANNArchy wichtig ? Was für Probleme gab es ?
Welche Schritte waren notewendig um es zu implementieren? Was für die implementierung relevantes
stand im Paper und was musste selbst 'entschlüsselt' werden?**

# Conclusion

Our reimplementation of voltage based STDP learning rule from @Clopath2010
is able to reproduce the experimental data and the emergent connectivity structures
proposed in the original paper.

The description of the learning rule in the original publication comprises enough details
to understand the different components and their interaction.
However, two main components have not been described adequately to allow a direct reimplementation:
the 'resolution trick' of the membrane potential after spike emission and the equation for the
homeostatic mechanism ($\bar{\bar{u}}$).  
The emergent of stable weights with high and low values depends on a functional homeostatic mechanism, as mentioned in the original publication [@Clopath2010].
But the formula to calculate the homeostatic variabel $\bar{\bar{u}}$ is not described in the publication.
The dependency of $\bar{\bar{u}}$ on the homeostatic mechanism is necessary to
implement the right behavior of the membrane potential.
The reimplementation has greatly benefited from the release of the source code on modelDB.

## Acknowledgment

This work was supported by the European Social Fund (ESF) and the Freistaat Sachsen.
# References
