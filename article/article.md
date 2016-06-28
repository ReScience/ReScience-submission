---            
Title: "Cellular and Network Mechanisms of Slow Oscillatory Activity (<1 Hz) 
        and Wave Propagations in a Cortical Network Model"
Author:
  - name: Andrei Maksimov
    affiliation: 1
  - name: Sacha J. van Albada
    affiliation: 1
  - name: Markus Diesmann
    affiliation: 1,2,3
Address:
  - code:    1
    address: Institute of Neuroscience and Medicine (INM-6) and Institute for Advanced Simulation (IAS-6) 
             and JARA BRAIN Institute I, Jülich Research Centre, 52428 Jülich, Germany
  - code:    2
    address: Department of Psychiatry, Psychotherapy and Psychosomatics, Medical Faculty, 
             RWTH Aachen University, 52062 Aachen, Germany
  - code:    3
    address: Department of Physics, Faculty 1, RWTH Aachen University, 52062 Aachen, Germany
    
Contact:
  - maksimov.andrei7@gmail.com
Editor:
  -
Reviewer:
  - 
  - 
Publication:
  received:  
  accepted:  
  published: 
  volume:   
  issue:    
  date:     
Repository:
  code:      
Reproduction:   
  - "*Cellular and network mechanisms of slow oscillatory activity (<1 Hz) 
  and wave propagations in a cortical network model*, A. Compte, M.V. Sanchez-Vives, 
  D.A. McCormick, X.-J. Wang, Journal of Neurophysiology, 2707--2725, 2003"
Bibliography:
  article.bib

---


# Introduction

We provide an implementation of the model of [@Compte03_2707], which reproduces single-neuron and collective network
behaviors during slow-wave oscillations in vitro in control conditions and under pharmacological manipulations.
In particular, we focus on the authors' model results that include: (a) neuronal membrane potentials oscillating
between Up and Down states at $<\!\!1\,\mathrm{Hz}$; (b) characteristic membrane resistance behavior and activation of neuronal
ion channels with proportional excitation and inhibition during Up states; (c) spontaneous initiation and further
wave-like propagation of population spiking activity. The original implementation is in
C\nolinebreak[4]\hspace{-0.05em}\raisebox{.4ex}{\tiny{\bf ++}}, but the source code is not publicly available. 
The implementation we propose is coded in the NEST [@Gewaltig_07_11204] framework, one of the modern actively 
developed simulation platforms that is publicly available. The code uses the Python interface [@Eppler09_12] for 
legibility. The model and analysis scripts are implemented using Python 2.7.6.


# Methods

We use the description of the model given in the original study, with the exception of the synaptic kinetics, which
is greatly simplified (see model description below). In the original model, the majority of parameters are given
per unit membrane surface area. For simplicity, we combine these parameters with the surface areas of the corresponding 
compartments in the following description. The resulting parameters are denoted with a tilde.


The model for excitatory neurons contains a somatic and a dendritic compartment. The somatic compartment includes fast
$\mathrm{Na}^{+}$ and $\mathrm{K}^{+}$ spiking currents ($\widetilde{I}_{\mathrm{Na}}$, $\widetilde{I}_{\mathrm{K}}$), 
a leak current $\widetilde{I}_{\mathrm{L}}$, a fast A-type $\mathrm{K}^{+}$ current $\widetilde{I}_{\mathrm{A}}$, a 
non-inactivating slow $\mathrm{K}^{+}$ current $\widetilde{I}_{\mathrm{KS}}$, and a $\mathrm{Na}{}^{+}$-dependent 
$\mathrm{K}^{+}$ current $\widetilde{I}_{\mathrm{KNa}}$. The dendritic compartment includes a high-threshold 
$\mathrm{Ca}{}^{2+}$ current $\widetilde{I}_{\mathrm{Ca}}$, a $\mathrm{Ca}{}^{2+}$-dependent $\mathrm{K}^{+}$ current 
$\widetilde{I}_{\mathrm{KCa}}$, a non-inactivating (persistent) $\mathrm{Na}{}^{+}$ current $\widetilde{I}_{\mathrm{NaP}}$, 
and an inward rectifier (activated by hyperpolarization) non-inactivating $\mathrm{K}^{+}$ current $\widetilde{I}_{\mathrm{AR}}$. 
Somato-dendritic coupling is implemented through the axial dendritic conductance $g_{\mathrm{sd}}$. Non-synaptic currents 
are modeled using the Hodgkin-Huxley formalism $\widetilde{I}\left(t\right)=\tilde{g}m^{k}h^{l}\left(V-E_{\mathrm{rev}}\right)$, 
where gating variables $m$ and $h$ are calculated using a first-order activation scheme, 
$\frac{dx}{dt}=\phi\left[\alpha\left(V\right)\left(1-x\right)-\beta\left(V\right)x\right]=\left(x_{\infty}\left(V\right)-x\right)/\tau\left(V\right)$
with $x_{\infty}=\frac{\alpha}{\alpha+\beta}$, $\tau=\frac{1}{\phi\left(\alpha+\beta\right)}$, and $\phi$ being the temperature factor 
(constant). In cases where time dependence is neglected (due to rapid activation or inactivation), gating variables are 
substituted by their saturation levels $m_{\infty}$ or $h_{\infty}$. The concentration of intrinsic neuronal $\mathrm{Ca}{}^{2+}$
and $\mathrm{Na}{}^{+}$ ions (in mM) is drawn from first-order differential equations 
$\frac{d\left[\mathrm{Ca}{}^{2+}\right]}{dt}=\alpha_{\mathrm{\mathrm{Ca}}}\widetilde{I}_{\mathrm{Ca}}-\left[\mathrm{Ca}{}^{2+}\right]/\tau_{\mathrm{Ca}}$
and $\frac{d\left[\mathrm{Na}{}^{+}\right]}{dt}=\alpha_{\mathrm{Na}}(\widetilde{I}_{\mathrm{Na}}+\widetilde{I}_{\mathrm{NaP}})-R_{\mathrm{pump}}
\left(\left[\mathrm{Na}{}^{+}\right]^{3}/\left(\left[\mathrm{Na}{}^{+}\right]^{3}+15^{3}\right)-\left[\mathrm{Na}{}^{+}\right]_{\mathrm{eq}}^{3}/
\left(\left[\mathrm{Na}{}^{+}\right]_{\mathrm{eq}}^{3}+15^{3}\right)\right)$
where $\alpha_{\mathrm{Ca}}$, $\alpha_{\mathrm{Na}}$, $\tau_{\mathrm{Ca}}$, $R_{\mathrm{pump}}$, $\left[\mathrm{Na}{}^{+}\right]_{\mathrm{eq}}$
are constants. The model for inhibitory neurons just consists of a somatic compartment with only $\widetilde{I}_{\mathrm{Na}}$, 
$\widetilde{I}_{\mathrm{\mathrm{K}}}$, and $\widetilde{I}_{\mathrm{L}}$ currents. AMPA and NMDA synaptic inputs target the 
dendritic compartment for excitatory and the somatic compartment for inhibitory neurons. GABA inputs always target the somatic compartment. 

Because of the complexity of the single-neuron model (with multiple voltage-dependent channels), we verify that our 
implementation behaves like the original according to the information given in the paper. Similar to the original work, $500\,\mathrm{ms}$
injection of $250\,\mathrm{pA}$ current into the soma (Fig.$\,$@fig:bm-neuron) results in an adapting firing pattern in pyramidal (PY) 
and non-adapting firing in fast-spiking (FS) neurons with average firing rates of $22$ and $76$ spikes/s, respectively. 
Distribution of the membrane leak conductance with 10\% (PY) and 2.5\% (FS) standard deviation around the mean value leads 
to a small fraction of pyramidal (but not fast-spiking) neurons showing pacemaking activity. 


![**Model response to 250 pA current injection for 500 ms into the soma of a pyramidal (PY) (A) and a fast-spiking (FS) (B) neuron.** 
The respective mean rates of 22 spikes/s for the PY neuron and 76 spikes/s for the FS neuron over the input period match those 
in the original implementation of Compte et al. [@Compte03_2707].](../code/fig1.png){#fig:bm-neuron}


The notable difference in our implementation is the synaptic kinetics. In the original model, synaptic gating variables 
depend on the presynaptic membrane potential. Implementing such a dependence is at present highly nontrivial in NEST due 
to optimizations for distributed computing (see Kunkel et al. [@Kunkel14_78], Hahne et al. [@Hahne15_00022] for background). 
For this reason, we simplify the synapse models while preserving the amplitude and shape of postsynaptic conductances. 
Specifically, the first-order kinetics of AMPA and GABA channels is replaced by a simple exponential decay 
$\frac{\partial g}{\partial t}=\underset{i}{\sum}W\delta\left(t-t_{i}\right)-\frac{g}{\tau}$, where g is the synaptic conductance, 
$\tau$ is the synaptic decay time, $W$ is the synaptic weight, and $i$ indexes the incoming spike times. The second-order kinetics 
of NMDA channels is substituted by a difference of slow and fast exponential components, 
$g^{\mathrm{NMDA}}=g_{\mathrm{slow}}^{\mathrm{NMDA}}-g_{\mathrm{fast}}^{\mathrm{NMDA}}$ with 
$\frac{dg_{\mathrm{slow/fast}}}{dt}=\underset{i}{\sum}W\delta\left(t-t_{i}\right)-\frac{g_{\mathrm{slow/fast}}}{\tau_{\mathrm{slow/fast}}}$. 
These simplifications are justified by the stereotyped trajectory of the presynaptic membrane potential during an action potential. 
Further, these simplifications allow us to merge gating variables $s$ of all synapses of one type (AMPA, NMDA or GABA) into a single 
postsynaptic time-dependent conductance (see Table @tbl:syn-values). While multiple activations of the 
same synapse lead to linear summation of postsynaptic AMPA and GABA conductances in the original model, NMDA conductances saturate at high 
input rates. Such behavior is modeled using the short-term plasticity formalism suggested by [@Markram98]. In short, the amplitude 
of the postsynaptic current $PSC=A\cdot R\cdot u$ is proportional to the fraction of available synaptic efficacy $R$ and 
utilization of synaptic efficacy $u$. Spike-triggered synaptic activation leads to a reduction of the available synaptic efficacy 
(corresponding to short-term depression) together with an increase in the utilization of synaptic efficacy (corresponding 
to short-term facilitation). In the time period $\Delta t$  between subsequent spikes, $R$ and $u$ recover to corresponding 
resting-state values ($R_{0}=1$ and $u_{0}=U$) with time constants $\tau_{\mathrm{rec}}$ and $\tau_{\mathrm{facil}}$. 
Our implementation closely follows the behavior of the original synapse model (Fig.$\,$@fig:bm-syn). Note that the equation 
describing the dynamics of the NMDA gating variable $s$ in the original paper contains a misprint and should read 

$$ \frac{ds}{dt}=\alpha(1-s)x-s/\tau. $$ {#eq:1}
 

![**Simplified synaptic kinetics (dashed curve) for AMPA, NMDA, and GABA conductances closely reproduces the behavior 
of the original model (solid curve) in a wide range of spiking rates.** DC currents, simultaneous injected into the soma of a pre-synaptic 
excitatory and an inhibitory neuron for $500\,\mathrm{ms}$, lead to synaptic activation at rates of $4$, $14$, $22$ spikes/s (top to bottom) 
for AMPA and NMDA, and $14$, $36$, $52$ spikes/s for GABA channels. These synaptic activation rates span the range of excitatory 
and inhibitory neuronal firing rates observed during Up states in the model. Curves for the original model obtained with the help of a 
reimplementation of the single-neuron and synapse dynamics in Python. To construct this figure, synaptic inputs with a 
weight of $W=1\,\mathrm{nS}$ were provided to AMPA (A), NMDA (B), and GABA (C) channels for the original implementation. Note that 
due to gating variable $s$ in the original implementation, the effective amplitude of a single synaptic input differs from 
the synaptic weight. For the implementation with simplified synaptic kinetics, synaptic weights are chosen to match the 
amplitude of the initial postsynaptic response in the original implementation.](../code/fig2.png){#fig:bm-syn}

The network architecture represents a chain of excitatory and inhibitory neurons (of length $l_{\mathrm{chain}}$) equidistantly 
distributed over $5\,\mathrm{mm}$. Each neuron projects a given number (drawn from a Gaussian distribution) of outgoing connections to 
target neurons of each type. The probability of a connection between any two neurons decays with inter-neuronal distance 
according to a Gaussian with characteristic scales $\lambda_{e}$ and $\lambda_{i}$ for excitatory and inhibitory pre-synaptic 
neuron types, respectively. Multiple connections can exist for a given pair of neurons. In Tables$\;$@tbl:Model-summary$\;\!$--@tbl:neuron-values 
we provide the description of the model.

The simulations are performed with NEST 2.8.0 [@Nest280] and combine an adaptive step size for the single-neuron solver with 
communication between neurons at a step size of $0.1\,\mathrm{ms}$. The time resolution of all recordings is $0.1\,\mathrm{ms}$. In Fig.$\,$@fig:result-intra$\,\!$E 
excitatory and inhibitory conductances are filtered with a $20\,\mathrm{ms}$ rectangular kernel. This kernel width is chosen to yield 
average input levels, rather than individual synaptic events. We find that filtering is required to reproduce the proportionality 
of excitatory and inhibitory conductances, shown in the original Figure 6, although no corresponding information is given in the original paper. 

To estimate the neuronal membrane conductance according to Eq.$\,$@eq:G-exp-def in the model at a time $t_{0}$ near the membrane potential $V_{0}$,
we use a procedure we refer to as the “virtual hyperpolarization method” (schematically shown in Fig.$\,$@fig:R-scheme). First, 
an isolated copy of a neuron model instance is created with its state identical to that of the original neuron at the time of interest $t_{0}$. 
Then, the neuron is allowed to relax to its equilibrium state, while being clamped to the corresponding membrane potential $V_{0}$. 
In this steady-state configuration, the total cross-membrane current $I_{0}$ is calculated as the sum of all channel currents 
according to Eq.$\,$@eq:I-as-sum. Then the neuron is slightly hyperpolarized to potential $V_{1}$ and the corresponding steady-state 
current $I_{1}$ is calculated. The ratio $\Delta{I}/\Delta{V}$ then gives an estimate of the membrane conductance at time $t_{0}$. 
This approach can be applied to the case where the membrane is approximately isopotential across the whole neuron (true for the present model). 
In case of more complex non-isopotential multi-compartment neurons, the corresponding equilibrium currents $I_{0}$ and $I_{1}$ should be 
estimated through direct simulation. Note that during this procedure the state of the synaptic input should be fixed at the level of time $t_{0}$. 

The relaxation time constant of certain ion concentrations (e.g., $80$--$350\,\mathrm{ms}$ for intracellular $\mathrm{Ca^{2+}}$ [@Markram1995_1], [@Schiller1995_583]) 
can be much larger than the membrane time constant (on average $<\!\!20\,\mathrm{ms}$ at rest, e.g., [@Levy2012_5609]) in cortical neurons. Therefore, the 
typical duration of DC pulse injection ($\sim\!\!80$ ms [@Contreras96], [@Waters06_8267]) used in experiments with active networks is enough 
to overcome the transient phase of prominent capacitive currents, while certain ion concentrations might be not equilibrated. 
In the present model, however, the membrane relaxation time constant reaches hundreds of milliseconds (see relaxation phase 
in Fig.$\,$@fig:IV), which makes it difficult to separate these two effects. Therefore, in our application of the virtual 
hyperpolarization method to the model we consider limit cases of both instantaneous (frozen) and equilibrated ion concentrations.

![**Schematic representation of the virtual hyperpolarization method.** An isolated copy of a neuron model instance is created with its 
state identical to that of the original neuron at the time of interest $t_{0}$, including the level of synaptic input. Then the 
cross-membrane current $I_{0}$, required to keep the neuron clamped to membrane potential $V_{0}$ in the steady-state scenario, 
is determined. After that, the neuron is slightly hyperpolarized to the potential $V_{1}$ and the corresponding steady-state current 
$I_{1}$ is estimated. The resulting ratio $\Delta I/\Delta V$ then gives the conductance estimate at time $t_{0}$. 
In case of an isopotential neuron model, the corresponding currents $I_{0}$ and $I_{1}$ can be calculated using 
Eq.$\,$@eq:I-as-sum .](../code/fig3.png){#fig:R-scheme}


\pagebreak



+--------------+---------------------------------------------------------------------------+
|Populations   |one excitatory and one inhibitory cortical population                      |
+--------------+---------------------------------------------------------------------------+
|Topology	   |one-dimensional (chain); Gaussian spatial connectivity profile             |
+--------------+---------------------------------------------------------------------------+
|Connectivity  |random connections with outdegree drawn from a Gaussian distribution       |
+--------------+---------------------------------------------------------------------------+
|Neuron model  |single or multi-compartment Hodgkin-Huxley-type model with multiple channel types|
+--------------+---------------------------------------------------------------------------+
|Channel model |Hodgkin-Huxley formalism                                                   |
+--------------+---------------------------------------------------------------------------+
|Synapse model |single- or double-exponential-shaped postsynaptic currents                 |
+--------------+---------------------------------------------------------------------------+
|Plasticity    |presynaptic short-term plasticity                                          |
+--------------+---------------------------------------------------------------------------+
|External input|None                                                                       |
+--------------+---------------------------------------------------------------------------+
|Recordings	   |Spike times from all neurons; membrane potential, $\mathrm{Na}^{+}$ and $\mathrm{Ca^{2+}}$  |
|              |concentrations, and all conductances and currents from a subset of neurons in both populations |
+--------------+---------------------------------------------------------------------------+

Table: **Model summary.** {#tbl:Model-summary}


\pagebreak


+----------------+--------+--------------------------------------------------------------------------------+
|Connectivity    |        |number (drawn from Gaussian distribution) of outgoing connections are randomly  |
|paradigm        |        |distributed across the target population with probability drawn from a Gaussian |
|                |        |inter-somatic distance-dependent profile; multiple connections for the same     |
|                |        |pair of neurons are allowed; autapses are forbidden                             |
+----------------+--------+--------------------------------------------------------------------------------+
|Synaptic weights|        |same for all connections of the same type                                       |
+----------------+--------+--------------------------------------------------------------------------------+
|Synaptic delays |        |same for all connections                                                        |
+----------------+--------+--------------------------------------------------------------------------------+
|Synaptic model	 |        |static synapse for connections to AMPA and GABA channels:            \          |
|                |        |$PSC_{i}=W$                                                          \          | 
|                |        |synapse with short-term plasticity according to [@Markram98] for NMDA channels:\| 
|	             |        |$PSC_{i+1}=W\cdot R_{i+1}\cdot u_{i+1}$                              \          |
|                |        |$R_{i+1}=1+\left(R_{i}-R_{i}u_{i}-1\right)\cdot\exp\left(-\delta t/\tau_{\mathrm{rec}}\right)$ \|
|                |        |$u_{i+1}=U+u_{n}\left(1-U\right)\exp\left(-\delta t/\tau_{\mathrm{fac}}\right)$ |
+----------------+--------+--------------------------------------------------------------------------------+
|$N_{e}$         |1024    |number of excitatory neurons                                                    |
+----------------+--------+--------------------------------------------------------------------------------+
|$N_{i}$         |256     |number of inhibitory neurons                                                    |
+----------------+--------+--------------------------------------------------------------------------------+
|$\lambda_{e}$   |250     |characteristic scale of spatial connectivity decay ($\mathrm{\mu m}$) for excitatory connections |
+----------------+--------+--------------------------------------------------------------------------------+
|$\lambda_{i}$   |125     |characteristic scale of spatial connectivity decay ($\mathrm{\mu m}$) for inhibitory connections |
+----------------+--------+--------------------------------------------------------------------------------+
|$l_{\mathrm{    |5000    |chain length ($\mathrm{\mu m}$)                                                 |
|chain}}$        |        |                                                                                |
+----------------+--------+--------------------------------------------------------------------------------+
|$\mathrm{       |$20\pm5$|mean and standard deviation of Gaussian distribution used to determine numbers of outgoing connections for each neuron|
|outdegree}$     |        |                                                                                |
+----------------+------+----------------------------------------------------------------------------------+

Table: **Network topology and synapse model.** {#tbl:Topo-syn}


\pagebreak


+---------------------------------+-----------------------------------------------------------------------------------------------------+
|Neuron model                     |Hodgkin-Huxley-type model with multiple channel types and exponential-based                          |
|                                 |synaptic conductances; two compartments for excitatory and one compartment for inhibitory neurons.   |
+---------------------------------+-----------------------------------------------------------------------------------------------------+
|Subthreshold\                    |soma:                                                                                                |
|dynamics \                       |$C_{\mathrm{m}}\frac{dV}{dt} = -\left(\widetilde{I}_{\mathrm{L}} + \widetilde{I}_{\mathrm{Na}} + \widetilde{I}_{\mathrm{K}} + \widetilde{I}_{\mathrm{A}} + \widetilde{I}_{\mathrm{KS}} + \widetilde{I}_{\mathrm{KNa}} + g_{\mathrm{sd}}\left(V_{\mathrm{d}} - V_{\mathrm{s}}\right) + I_{\mathrm{GABA}}\right)$ \   |
|for PY   \                       |dendrite:                                                                                                                                                                                                                                                                                                           |
|neurons  \                       |$C_{\mathrm{m}}\frac{dV}{dt}=-\left(\widetilde{I}_{\mathrm{Ca}} + \widetilde{I}_{\mathrm{KCa}} + \widetilde{I}_{\mathrm{NaP}} + \widetilde{I}_{\mathrm{AR}} + g_{\mathrm{sd}}\left(V_{\mathrm{s}} - V_{\mathrm{d}}\right) + I_{\mathrm{AMPA}} + I_{\mathrm{NMDA}}\right)$                                           |
+---------------------------------+-----------------------------------------------------------------------------------------------------+
|Subthreshold   \                 |soma:                                                                                                |
|dynamics for   \                 |$C_{\mathrm{m}}\frac{dV}{dt} = -\left(\widetilde{I}_{\mathrm{L}} + \widetilde{I}_{\mathrm{Na}} + \widetilde{I}_{\mathrm{K}} + I_{\mathrm{AMPA}} + I_{\mathrm{NMDA}} + I_{\mathrm{GABA}}\right)$ \|
|FS neurons     \                 |                                                                                                     |
+---------------------------------+-----------------------------------------------------------------------------------------------------+
|Spike detection                  |A spike is detected when the somatic membrane potential rises above $0\,\mathrm{mV}$                 |
|                                 |and its derivative becomes negative: $\left(V_{s}>0\right)\;\land\;\left(\frac{dV_{s}}{dt}<0\right)$.| 
|                                 |After that, the neuron becomes refractory and no spike emission is allowed during a fixed time.      |
+---------------------------------+-----------------------------------------------------------------------------------------------------+
|Post-synaptic  \                 |$g_{\mathrm{AMPA},\mathrm{GABA}}\left(t\right)=w\exp\left(-t/\tau\right)$ \                          |
|conductances	                  |$g_{\mathrm{NMDA}}\left(t\right)=w\left(\exp\left(-t/\tau_{\mathrm{slow}}\right)-\exp\left(-t/\tau_{\mathrm{fast}}\right)\right)$ |
+---------------------------------+-----------------------------------------------------------------------------------------------------+
|Channel dynamics                 |Hodgkin-Huxley formalism $\widetilde{I}\left(t\right)=\tilde{g}m^{k}h^{l}\left(V-E_{\mathrm{rev}}\right)$; \|
|                                 |$m,\,h$ follow a first-order activation scheme, $\frac{dx}{dt}=\phi\left[\alpha\left(V\right)\left(1-x\right)-\beta\left(V\right)x\right]=\left(x_{\infty}\left(V\right)-x\right)/\tau\left(V\right)$ with $x_{\infty}=\frac{\alpha}{\alpha+\beta}$, $\tau=\frac{1}{\phi\left(\alpha+\beta\right)}$ |
+---------------------------------+-----------------------------------------------------------------------------------------------------+
|$\mathrm{Ca}^{2+}$ concentration | $\frac{d\left[\mathrm{Ca}^{2+}\right]}{dt}=\alpha_{\mathrm{Ca}}\widetilde{I}_{\mathrm{Ca}} - \left[\mathrm{Ca}^{2+}\right]/\tau_{\mathrm{Ca}}$ |
+---------------------------------+-----------------------------------------------------------------------------------------------------+
|$\mathrm{Na}^{+}$ concentration  | $\frac{d\left[\mathrm{Na}^{+}\right]}{dt}=\alpha_{\mathrm{Na}}(\widetilde{I}_{\mathrm{Na}} + \widetilde{I}_{\mathrm{NaP}}) - R_{\mathrm{\mathrm{pump}}}\left(\left[\mathrm{Na}^{+}\right]^{3}/\left(\left[\mathrm{Na}^{+}\right]^{3}+15^{3}\right)-\left[\mathrm{Na}^{+}\right]_{\mathrm{eq}}^{3}/\left(\left[\mathrm{Na}^{+}\right]_{\mathrm{eq}}^{3}+15^{3}\right)\right)$ |
+---------------------------------+-----------------------------------------------------------------------------------------------------+

Table: **Neuron model and post-synaptic conductances.** Ionic concentrations are measured in $\mathrm{mM}$. {#tbl:neuron-general}


\pagebreak


+-----------------------------------------+-------------------------------------------------------------------------------------+
|**Fast sodium current, soma**          \ |activation variable $m$:   \                                                           | 
|$\widetilde{I}_{\mathrm{Na}}=\tilde{g}   |$\alpha=0.1\left(V+33\right)/\left(1-\exp\left(-\left(V+33\right)/10\right)\right)$ \|
|_{\mathrm{Na}}m^{3}h\left(V-E_           |$\beta=4\exp(-(V+53.7)/12)$ \                                                        |
|{\mathrm{Na}}\right)$                    |inactivation variable $h$: \                                                         |
|                                         |$\alpha=0.07\exp\left(-\left(V+50\right)/10\right)$ \                                |
|                                         |$\beta=1/\left(1+\exp\left(-\left(V+20\right)/10\right)\right)$ \                    |
|                                         |$\tau=\frac{1}{4\left(a+b\right)}$                                                   |
+-----------------------------------------+-------------------------------------------------------------------------------------+
|**Fast potassium current, soma**       \ |inactivation variable $k$: \                                                         | 
|$\widetilde{I}_{\mathrm{K}}=\tilde{g}    |$\alpha=0.01\cdot\left(V+34\right)/\left(1-\exp\left(-\left(V+34\right)/10\right)\right)$ \|
|_{\mathrm{K}}k^{4}\left(V-E_             |$\beta=0.125\cdot\exp\left(-\left(V+44\right)/25\right)$ \                           |
|{\mathrm{K}}\right)$                     |$\tau=\frac{1}{4\left(a+b\right)}$                                                   |
+-----------------------------------------+-------------------------------------------------------------------------------------+
|**Leakage current, soma**  \             |                                                                                     | 
|$\widetilde{I}_{\mathrm{L}}=\tilde{g}    |                                                                                     |
|_{\mathrm{L}}\left(V-E_                  |                                                                                     |
|{\mathrm{L}}\right)$                     |                                                                                     |
+-----------------------------------------+-------------------------------------------------------------------------------------+
|**Fast A-type $\mathrm{K}^{+}$           |activation variable $m$:   \                                                         | 
|current, soma  **\                       |$m_{\infty}=1/\left(1+\exp\left(-\left(V+50\right)/20\right)\right)$          \      |
|$\widetilde{I}_{\mathrm{A}}=\tilde{g}    |inactivation variable h: \                                                           |
|_{\mathrm{A}}m^{3}h\left(V-E_            |$h_{\infty}=1/\left(1+\exp\left(\left(V+80\right)/6\right)\right)$           \       |
|{\mathrm{K}}\right)$                     |$\tau=15$                                                                            |
+-----------------------------------------+-------------------------------------------------------------------------------------+
|**Non-inactivating $\mathrm{K}^{+}$      |activation variable $m$:   \                                                         |
|current, soma **\                        |$m_{\infty}=1/\left(1+\exp\left(-\left(V+34\right)/6.5\right)\right)$         \      |
|$\widetilde{I}_{\mathrm{KS}}=\tilde{g}   |$\tau=8/\left(\exp\left(-\left(V+55\right)/30\right)+\exp\left(\left(V+55\right)/30\right)\right)$ |
|_{\mathrm{KS}}m\left(V-E_                |                                                                                     |
|{\mathrm{K}}\right)$                     |                                                                                     |
+-----------------------------------------+-------------------------------------------------------------------------------------+
|**Persistent $\mathrm{Na}^{+}$           |activation variable $m$:   \                                                         |
|current, dendrite**\                     |$m_{\infty}=1/\left(1+\exp\left(-\left(V+55.7\right)/7.7\right)\right)$              |
|$\widetilde{I}_{\mathrm{Na}}=\tilde{g}   |                                                                                     |
|_{\mathrm{NaP}}m^{3}\left(V-E_           |                                                                                     |
|{\mathrm{Na}}\right)$                    |                                                                                     |
+-----------------------------------------+-------------------------------------------------------------------------------------+
|**Inwardly rectifying $\mathrm{K}^{+}$   |activation variable $m$:   \                                                         |
|current, dendrite **\                    |$m_{\infty}=1/\left(1+\exp\left(\left(V+75\right))/4\right)\right)$                  |
|$\widetilde{I}_{\mathrm{AR}}=\tilde{g}   |                                                                                     |
|_{\mathrm{AR}}m\left(V-E_                |                                                                                     |
|{\mathrm{K}}\right)$                     |                                                                                     |
+-----------------------------------------+-------------------------------------------------------------------------------------+
|**High-threshold $\mathrm{Ca}^{2+}$      |activation variable $m$:   \                                                         |
|current, dendrite **\                    |$m_{\infty}=1/\left(1+\exp\left(-\left(V+20\right)/9\right)\right)$                  |
|$\widetilde{I}_{\mathrm{Ca}}=\tilde{g}   |                                                                                     |
|_{\mathrm{Ca}}m^{2}\left(V-E_            |                                                                                     |
|{\mathrm{\mathrm{Ca}}}\right)$           |                                                                                     |
+-----------------------------------------+-------------------------------------------------------------------------------------+
|**$\mathrm{Na}^{+}$-dependent            |activation variable $m$:   \                                                         |
|$\mathrm{K^{+}}$                         |$m_{\infty}=0.37/\left(1+\left(38.7/\left[\mathrm{Na}^{+}\right]\right)^{3.5}\right)$|                                               |
|current, soma** \                        |                                                                                     |
|$\widetilde{I}_{\mathrm{KNa}}=\tilde{g}  |                                                                                     |
|_{\mathrm{KNa}}m\left(V-E_               |                                                                                     |
|{\mathrm{K}}\right)$                     |                                                                                     |
+-----------------------------------------+-------------------------------------------------------------------------------------+
|**$\mathrm{Ca}^{2+}$-dependent           |activation variable $m$:   \                                                         |
|$\mathrm{K^{+}}$                         |$m_{\infty}=\left[\mathrm{Ca}^{2+}\right]/\left(\left[\mathrm{Ca}^{2+}\right]+30\right)$|
|current, dendrite** \                    |                                                                                     |
|$\widetilde{I}_{\mathrm{KCa}}=\tilde{g}  |                                                                                     |
|_{\mathrm{KCa}}m\left(V-E_               |                                                                                     |
|{\mathrm{K}}\right)$                     |                                                                                     |
+-----------------------------------------+-------------------------------------------------------------------------------------+
 
Table: **Channel dynamics for excitatory neurons.** Membrane and reversal potentials are measured in $\mathrm{mV}$, conductances in $\mathrm{nS}$, 
ionic concentrations in $\mathrm{mM}$, and time constants in $\mathrm{ms}$. {#tbl:neuron-channels-ex}


\pagebreak


+-----------------------------------------+-------------------------------------------------------------------------------------+
|**Fast sodium current, soma**          \ |activation variable $m$:   \                                                         | 
|$\widetilde{I}_{\mathrm{Na}}=\tilde{g}   |$\alpha=0.5\left(V+35\right)/\left(1-\exp\left(-\left(V+35\right)/10\right)\right)$ \|
|_{\mathrm{Na}}m^{3}h\left(V-E_           |$\beta=20\exp\left(-\left(V+60\right)/18\right)$ \                                   |
|{\mathrm{Na}}\right)$                    |inactivation variable $h$: \                                                         |
|                                         |$\alpha=0.35\exp\left(-\left(V+58\right)/20\right)$ \                                |
|                                         |$\beta=5/\left(1+\exp\left(-\left(V+28\right)/10\right)\right)$ \                    |
|                                         |$\tau=\frac{1}{\left(a+b\right)}$                                                    |
+-----------------------------------------+-------------------------------------------------------------------------------------+
|**Fast potassium current, soma**       \ |inactivation variable $k$: \                                                         | 
|$\widetilde{I}_{\mathrm{K}}=\tilde{g}    |$\alpha=0.05\left(V+34\right)/\left(1-\exp\left(-\left(V+34\right)/10\right)\right)$\|
|_{\mathrm{K}}k^{4}\left(V-E_             |$\beta=0.625\exp\left(-\left(V+44\right)/80\right)$ \                                |
|{\mathrm{K}}\right)$                     |$\tau=\frac{1}{\left(a+b\right)}$                                                    |
+-----------------------------------------+-------------------------------------------------------------------------------------+
|**Leakage current, soma**  \             |                                                                                     | 
|$\widetilde{I}_{\mathrm{L}}=\tilde{g}    |                                                                                     |
|_{\mathrm{L}}\left(V-E_                  |                                                                                     |
|{\mathrm{L}}\right)$                     |                                                                                     |
+-----------------------------------------+-------------------------------------------------------------------------------------+

Table: **Channel dynamics for inhibitory neurons.** Membrane and reversal potentials are measured in $\mathrm{mV}$, conductances in $\mathrm{nS}$, 
ionic concentrations in $\mathrm{mM}$, and time constants in $\mathrm{ms}$. {#tbl:neuron-channels-in}



+-----------------------+--------+--------------------------------------------------------------------------------------------------------------------------------+
|$\mathrm{delay}$       |0.1     |synaptic delay ($\mathrm{ms}$)                                                                                                  |
+-----------------------+--------+--------------------------------------------------------------------------------------------------------------------------------+
|$W_{e\leftarrow        |2.72    |synaptic weights ($\mathrm{nS}$)                                                                                                |
|\mathrm{AMPA}}$        |        |                                                                                                                                |
+-----------------------+--------+--------------------------------------------------------------------------------------------------------------------------------+
|$W_{i\leftarrow        |0.13    |                                                                                                                                |
|\mathrm{AMPA}}$        |        |                                                                                                                                |
+-----------------------+--------+--------------------------------------------------------------------------------------------------------------------------------+
|$W_{e\leftarrow        |0.18    |                                                                                                                                |
|\mathrm{NMDA}}$        |        |                                                                                                                                |
+-----------------------+--------+--------------------------------------------------------------------------------------------------------------------------------+
|$W_{i\leftarrow        |0.03    |                                                                                                                                |
|\mathrm{NMDA}}$        |        |                                                                                                                                |
+-----------------------+--------+--------------------------------------------------------------------------------------------------------------------------------+
|$W_{e\leftarrow        |2.59    |                                                                                                                                |
|\mathrm{GABA}}$        |        |                                                                                                                                |
+-----------------------+--------+--------------------------------------------------------------------------------------------------------------------------------+
|$W_{i\leftarrow        |0.07    |                                                                                                                                |
|\mathrm{GABA}}$        |        |                                                                                                                                |
+-----------------------+--------+--------------------------------------------------------------------------------------------------------------------------------+
|$U$                    |0.01    |initial utilization of synaptic efficacy                                                                                        |
+-----------------------+--------+--------------------------------------------------------------------------------------------------------------------------------+
|$\tau_{\mathrm{rec}}$  |130     |recovery time constant ($\mathrm{ms}$) of available synaptic efficacy for NMDA synapses onto PY and FS neurons                  |
+-----------------------+--------+--------------------------------------------------------------------------------------------------------------------------------+
|$\tau_{\mathrm{facil}}$|0       |recovery time constant ($\mathrm{ms}$) of utilization of synaptic efficacy                                                      |
+-----------------------+--------+--------------------------------------------------------------------------------------------------------------------------------+
|$\tau_{\mathrm{AMPA}}$ |2       |time constant of AMPA channels ($\mathrm{ms}$)                                                                                  |
+-----------------------+--------+--------------------------------------------------------------------------------------------------------------------------------+
|$\tau_{\mathrm{GABA}}$ |10      |time constant of GABA channels ($\mathrm{ms}$)                                                                                  |
+-----------------------+--------+--------------------------------------------------------------------------------------------------------------------------------+
|$\tau_{\mathrm{NMDA}}  |100     |slow time constant of NMDA channels ($\mathrm{ms}$)                                                                             |
|^{\mathrm{slow}}$      |        |                                                                                                                                |
+-----------------------+--------+--------------------------------------------------------------------------------------------------------------------------------+
|$\tau_{\mathrm{NMDA}}  |2       |fast time constant of NMDA channels ($\mathrm{ms}$)                                                                             |
|^{\mathrm{fast}}$      |        |                                                                                                                                |
+-----------------------+--------+--------------------------------------------------------------------------------------------------------------------------------+

Table: **Parameter specification for our synaptic implementation.** Note that the synaptic conductances differ with respect to the original implementation. {#tbl:syn-values} 


\pagebreak


+-------------+----------------+-----------------------------------------------------------------+
|$E_{\mathrm  |120             |reversal potential of $\mathrm{Ca}^{2+}$ channels ($\mathrm{mV}$)|
|{Ca}}$       |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$E_{\mathrm  |55              |reversal potential of $\mathrm{Na}^{+}$ channels ($\mathrm{mV}$) |
|{Na}}$       |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$E_{\mathrm  |$-100(-90)$     |reversal potential of $\mathrm{K}^{+}$ channels ($\mathrm{mV}$)  |
|{K}}$        |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$E_{\mathrm  |$-60.95\pm0.3$  |leak reversal potential ($\mathrm{mV}$), mean and standard deviation of Gaussian distribution|
|{L}}$        |($-63.8\pm0.15$)|                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$E_{\mathrm  |0               |reversal potential AMPA and NMDA channels ($\mathrm{mV}$)        |
|{ex}}$       |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$E_{\mathrm  |$-70$           |reversal potential GABA channels ($\mathrm{mV}$)                 |
|{in}}$       |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$\tilde{C}_  |150(200)        |somatic membrane capacitance ($\mathrm{pF}$)                     |
|{\mathrm{m}  |                |                                                                 |
|}^{\mathrm   |                |                                                                 |
|{soma}}$     |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$\tilde{C}_  |350(-)          |dendritic membrane capacitance ($\mathrm{pF}$)                   |
|{\mathrm{m}  |                |                                                                 |
|}^{\mathrm   |                |                                                                 |
|{dendr}}$    |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$\tilde{g}_  |7500            |maximal conductance of fast $\mathrm{Na}^{+}$ channel ($\mathrm{nS}$)|
|{\mathrm{Na} |(7000)          |                                                                 |
|}^{\mathrm   |                |                                                                 |
|{soma}}$     |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$\tilde{g}_  |1575            |maximal conductance of fast $\mathrm{K}^{+}$ channel ($\mathrm{nS}$)|
|{\mathrm{K}  |(1800)          |                                                                 |
|}^{\mathrm   |                |                                                                 |
|{soma}}$     |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$\tilde{g}_  |$10\pm1$        |leak conductance ($\mathrm{nS}$)                                 |
|{\mathrm{L}  |($20.5\pm0.05$) |                                                                 |
|}^{\mathrm   |                |                                                                 |
|{soma}}$     |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$\tilde{g}_  |150             |maximal conductance of A-type fast $\mathrm{K}^{+}$ channel ($\mathrm{nS}$)|
|{\mathrm{KA} |(-)             |                                                                 |
|}^{\mathrm   |                |                                                                 |
|{soma}}$     |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$\tilde{g}_  |200             |maximal conductance of $\mathrm{Na}^{+}$-dependent $\mathrm{K}^{+}$ channel ($\mathrm{nS}$) |
|{\mathrm{KNa}|(-)             |                                                                 |
|}^{\mathrm   |                |                                                                 |
|{soma}}$     |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$\tilde{g}_  |86.4            |maximal conductance of non-inactivating $\mathrm{K}^{+}$ channel ($\mathrm{nS}$) |
|{\mathrm{KS} |(-)             |                                                                 |
|}^{\mathrm   |                |                                                                 |
|{soma}}$     |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$\tilde{g}_  |200             |maximal conductance of $\mathrm{Ca}^{2+}$-dependent $\mathrm{K}^{+}$ channel ($\mathrm{nS}$) |
|{\mathrm{KCa}|(-)             |                                                                 |
|}^{\mathrm   |                |                                                                 |
|{dendr}}$    |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$\tilde{g}_  |9               |maximal conductance of inwardly rectifying $\mathrm{K}^{+}$ channel ($\mathrm{nS}$) |
|{\mathrm{KAR}|(-)             |                                                                 |
|}^{\mathrm   |                |                                                                 |
|{dendr}}$    |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$\tilde{g}_  |24              |maximal conductance of persistent $\mathrm{Na}^{+}$ channel ($\mathrm{nS}$)  |
|{\mathrm{NaP}|(-)             |                                                                 |
|}^{\mathrm   |                |                                                                 |
|{dendr}}$    |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$\tilde{g}_  |150.5           |maximal conductance of high-threshold $\mathrm{Ca}^{2+}$ channel ($\mathrm{nS}$) |
|{\mathrm{Ca} |(-)             |                                                                 |
|}^{\mathrm   |                |                                                                 |
|{dendr}}$    |                |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+
|$g_{         |1750            |axial dendritic conductance between somatic and dendritic compartments  |
|\mathrm{ax}}$|(-)             |                                                                 |
+-------------+----------------+-----------------------------------------------------------------+

Table: **Parameter specification for our neuronal implementation.** The superscripts $\mathrm{soma}$ and $\mathrm{dendr}$ refer 
to the somatic and dendritic compartments. When parameters for excitatory and inhibitory neurons are different, parameters 
for the latter are given in brackets. "(-)" means parameter not used for the inhibitory neuron model. {#tbl:neuron-values} 

# Results

We here focus on the main activity regime of the model, namely spontaneously generated Up-Down oscillations, and do not reproduce 
the network response to pharmacological manipulations. First, we simulate the reimplemented model with the parameters given in the 
original paper. However, despite close reproduction of the model constituents, the simulated network activity is characterized by 
unreasonably high firing rates due to the dominance of NMDA conductances far exceeding the potency of the opposing GABAergic inhibition. 
To achieve a network state characterized by comparable excitation and inhibition, we modify the synaptic strengths $W$ of all synapse types 
(values are given in Table @tbl:syn-values). This results in periods of spontaneously generated activity (Up states) which propagate 
along the network in a wave-like fashion (Fig.$\,$@fig:result-spikes$\,\!$A) with a mean firing rate around 20 spikes/s for excitatory and 35 spikes/s
for inhibitory neurons at the peak of activity (Fig.$\,$@fig:result-spikes$\,\!$B). On the level of individual neurons, Up states are characterized 
by a depolarized membrane potential and an increase of the intracellular concentration of $\mathrm{Na}^{+}$ ions by $2$--$4\,\mathrm{mM}$ 
(Fig.$\,$@fig:result-spikes$\,\!$C,D for two neurons), similar to the results shown in Figure 2 of the original work. Synaptic and intrinsic 
neuronal conductances (Fig.$\,$@fig:result-intra$\,\!$B--I) are tightly coupled to the membrane depolarization (Fig.$\,$@fig:result-intra$\,\!$A) 
and show a dynamic range very similar to that reported in the original figure. When filtered with a rectangular kernel (see Methods), 
excitatory and inhibitory conductances are clearly coupled (Fig.$\,$@fig:result-intra$\,\!$E), like in the original model (see Figure 6). 

When computed as the reciprocal of the sum of all conductances (see Eq.$\,$@eq:G-model-def), membrane resistance is decreased during Up states 
relative to silent Down states (Fig.$\,$@fig:result-intra$\,\!$B,$\,$@fig:result-R), consistent with the results shown in Figure 5B of the original work. 
We also estimated membrane resistance by injecting hyperpolarizing DC pulses ($300\,\mathrm{pA}$ for $100\,\mathrm{ms}$) into neurons, constantly hyperpolarized by $250\,\mathrm{pA}$
current injection (see Fig.$\,$@fig:result-R$\,\!$A), similar to the way described in the original Figure 4. In agreement with the original paper, 
the resulting membrane resistance, computed as the voltage deflection at the end of the pulse divided by the pulse amplitude (Fig.$\,$@fig:result-R$\,\!$B 
dots) is quantitatively similar to the values obtained with the reciprocal of Eq.$\,$@eq:G-model-def (Fig.$\,$@fig:result-R solid line) at 
least for membrane potentials below $-80\,\mathrm{mV}$.


![**Spontaneous Up-Down oscillations generated in a network of 1280 neurons with a chain-like architecture.** **(A)** Spiking activity of excitatory 
(red) and inhibitory (blue) neurons propagates along the chain in a wave-like fashion. **(B)** Population firing rates show clear Up states, 
synchronized across excitatory and inhibitory neurons. Similar to Figure 3 of the original paper. **(C,D)** Intracellularly recorded membrane potentials 
(top) and concentration of intracellular $\mathrm{Na}^{+}$ ions (bottom) of PY neurons are similar to those reported in the original study (Figure 2). 
The neuron in (C) shows pacemaking activity with spiking during Down states, while the neuron in (D) is typical for the majority of neurons and is 
characterized by well-defined Up and Down states.](../code/fig4.png){#fig:result-spikes}


![**Membrane input resistance and various ionic conductances in the course of the slow oscillation on the example of a representative neuron.**
**(A)** Membrane voltage trace shows periods of high activity (Up states). **(B)** Total input resistance, measured as the reciprocal of the summed 
open membrane conductances. Excitatory and inhibitory synaptic conductances **(C, D)** are approximately proportional when binned with $20\,\mathrm{ms}$ bin 
width **(E)**. **(F-J)** The dynamics of various other conductances ($\tilde{g}_{\mathrm{NaP}}$; $\tilde{g}_{\mathrm{KS}}$; $\tilde{g}_{\mathrm{KAR}}$; 
$\tilde{g}_{\mathrm{KCa}}$; $\tilde{g}_{\mathrm{KNa}}$) closely resembles that reported in the original paper (Figure 5).](../code/fig5.png){#fig:result-intra}


![Accessing neuronal membrane resistance through the reciprocal of the sum of open channel conductances (gray plot in **B**) and through the injection 
of brief hyperpolarizing pulses (black dots in **B**) results in the quantitatively similar estimates. Hyperpolarizing pulses have an amplitude of $300\,\mathrm{pA}$
and duration of $100\,\mathrm{ms}$ **(A)**, while the neuron is hyperpolarized by a $250\,\mathrm{pA}$ current, similarly to the procedure described in the original 
Figure 4.](../code/fig6.png){#fig:result-R}


However, note that the authors define membrane conductance as 

$$ G_{\mathrm{m}}^{\mathrm{model}}=\underset{i}{\sum}G_{i}\left(V\right) $$ {#eq:G-model-def} 

\noindent where $G_{i}$ are the instantaneous conductances of all ionic channels, which are typically time- and voltage-dependent. This definition corresponds 
to the “instantaneous chord conductance” according to the classification given by [@Jack1983]. In contrast, the typical experimentally measured quantity 
corresponds to “steady-state slope conductance” from the same classification: 

$$ G_{\mathrm{m}}^{\mathrm{exp}}=\frac{\Delta I}{\Delta V} $$ {#eq:G-exp-def}

\noindent where $\Delta I$ is the extra current (injected into the soma) required to achieve a steady-state membrane potential shift $\Delta V=V-V_{0}$ from 
the initial level $V_{0}$ (prior to current injection). The steady state here refers to the situation where transient capacitive currents 
become negligible. To compare these two definitions we use the fact that in the steady-state case, the externally injected current $I$ is equal 
to minus the sum of all ionic currents, 

$$ I=\underset{i}{\sum}G_{i}\left(V\right)\cdot\left(V-E_{i}\right), $$ {#eq:I-as-sum}

\noindent where $E_{i}$ is the reversal potential of channel $i$. Combining Eqs$\,$@eq:G-exp-def and @eq:I-as-sum and using 

$$\begin{aligned}\Delta\left(G_{i}\left(V\right)\cdot\left(V-E_{i}\right)\right)&=\Delta\left(G_{i}\left(V\right)\right)\cdot\left(V-E_{i}\right)+
G_{i}\left(V\right)\cdot\Delta\left(V-E_{i}\right)\\\nonumber
&=\Delta\left(G_{i}\left(V\right)\right)\cdot\left(V-E_{i}\right)+G_{i}\left(V\right)\cdot\Delta V\end{aligned}$$
results in an overall membrane conductance 

$$ G_{\mathrm{m}}^{\mathrm{exp}}=\underset{i}{\sum}\frac{\Delta\left(G_{i}\left(V\right)\cdot\left(V-E_{i}\right)\right)}{\Delta V}=
\underset{i}{\sum}\frac{\Delta\left(G_{i}\left(V\right)\right)}{\Delta V}\left(V-E_{i}\right)+\underset{i}{\sum}G_{i}\left(V\right). $$ {#eq:G-m-def}

\noindent As one can see, the definition used in the original study (Eq.$\,$@eq:G-model-def) excludes the first term. To demonstrate the difference in 
these two definitions, we simulate isolated pyramidal neurons with DC inputs of different amplitudes injected for 20 seconds into the soma 
record the membrane potential (Fig.$\,$@fig:IV$\,\!$A) as well as all ionic channel conductances. A long pulse duration is chosen here to take into account 
the long time required for the model neuron to reach a steady state. Then we compute membrane conductance according to Eq.$\,$@eq:G-m-def and Eq.$\,$@eq:G-model-def
a few milliseconds before the DC input is switched off. The I-V curve shows a profound non-linearity (Fig.$\,$@fig:IV$\,\!$B) leading to $G_{\mathrm{m}}^{\mathrm{exp}}$
(the reciprocal of the slope of the I-V curve) reaching $0.012\,\mathrm{nS}$ (corresponding to $83\,\mathrm{G}\Omega$), while $G_{\mathrm{m}}^{\mathrm{model}}$ remains near $15\,\mathrm{nS}$
(corresponding to $0.07\,\mathrm{G}\Omega$) (Fig.$\,$@fig:IV$\,\!$C). Therefore, under the protocol mimicking a standard electrophysiological procedure, the instantaneous chord 
and steady-state slope conductances deviate by several orders of magnitude for the present neuron model. To measure the steady-state slope conductance 
(Eq.$\,$@eq:G-exp-def) during ongoing network activity in the model, we use the virtual hyperpolarization method (see Methods). Freezing the state of the neuron 
as done in this method is necessary because of the aforementioned long relaxation time of the model neurons, which necessitates pulse durations that exceed 
the length of an Up state and therefore prevents measuring steady-state resistances during freely evolving activity. The results of this method for the above 
case of prolonged DC stimulation closely match those obtained from the I-V curve (Fig.$\,$@fig:IV$\,\!$C). In the active network model, this method results in a membrane 
conductance around $2\,\mathrm{nS}$ (corresponding to $500\,\mathrm{M}\Omega$) during the Down state (Fig.$\,$@fig:R-mod$\,\!$B,C). During the Up states, however, membrane conductance tends 
to drop to zero and even becomes negative, indicating self-depolarizing dynamics in the subthreshold periods. Note that in Fig.$\,$@fig:R-mod$\,\!$B,C time periods of $20\,\mathrm{ms}$
before and after each spike event (shown by vertical gray lines) are not considered to avoid the contamination of subthreshold dynamics with spikes.
  

![**The method of accessing neuronal membrane conductance (or resistance) suggested by [@Compte03_2707] deviates from the approach based on the 
construction of an I-V curve.** **(A)** A $10\,\mathrm{s}$ DC injection (inside the region marked by vertical dashed lines) with various amplitudes into the 
somatic compartment of a PY neuron results in neuronal hyper-/depolarization. Relaxation of membrane potentials after offset of the DC input takes 
hundreds of ms. **(B)** The I-V curve shows a strong voltage dependence of the neuronal resistance (measured as the slope of the curve). **(C)**
Membrane conductance, measured with the method suggested by the authors (Eq.$\,$@eq:G-model-def; dashed curve), is around $15\,\mathrm{nS}$ (corresponding to a 
$70\,\mathrm{M}\Omega$ resistance). Membrane conductance, measured according to Eq.$\,$@eq:G-m-def (solid curve), reaches $0.012\,\mathrm{nS}$ (corresponding to a $83\,\mathrm{G}\Omega$
resistance). The virtual hyperpolarization method (circles, see Methods) closely reproduces the latter conductance measure. ](../code/fig7.png){#fig:IV}


![**Application of the virtual hyperpolarization method to the network model.** Neuronal depolarization during Up states **(A)** is associated with 
dominant strongly negative membrane conductance **(B,C)**, which is the result of self-depolarization due to intrinsic neuronal dynamics. Membrane 
conductance is calculated according to Eq.$\,$@eq:G-exp-def with $\mathrm{Na^{+}}$ and $\mathrm{Ca^{2+}}$ frozen **(B)** and steady-state **(C)**. 
Vertical gray lines show the time window of $20\,\mathrm{ms}$ just before and after individual spikes, which is masked in the determination of the conductance 
to avoid contamination by supra-threshold activity. ](../code/fig8.png){#fig:R-mod}

 
# Conclusion

After some minor modifications to the model, we are able to reproduce the majority of the original results concerning spontaneously generated Up-Down 
oscillations. In particular, the dynamics of membrane potentials, membrane resistances, channel conductances, intracellular concentrations of $\mathrm{Na}^{+}$
and $\mathrm{Ca}^{2+}$ ions, and spiking activity closely resemble those of the original model. In our implementation, the synaptic dynamics does not 
depend on the shapes of individual presynaptic action potentials, but just incorporates the average post-synaptic effect of an action potential. The 
fact that we can reproduce the emerging network phenomena suggests that this detail of the synapse model is not relevant on the network level. We provide 
a closer look at the method which the authors of the original study use to access neuronal conductance (or resistance). Their measure reflects an 
“instantaneous chord conductance”, which results in values around $20\,\mathrm{nS}$ throughout network activity. In experimental works, the typical measure approximates 
a “steady-state slope conductance” (see [@Jack1983] for classification), which results in a membrane conductance of around $2\,\mathrm{nS}$ during Down states and 
dominant negative conductances during Up states in the model. The latter suggests powerful non-synaptic self-excitation mechanisms, which dominate during 
subthreshold neuronal activation, deviating from typical results reported in the experimental literature. Therefore, the dynamics of membrane conductance 
in the model merits further investigation.

Overall, we confirm the correctness of the original implementation of the model. 


# Acknowledgements

Partly supported by by Helmholtz Portfolio Supercomputing and Modeling for the Human Brain (SMHB), EU Grant 269921 (BrainScaleS), and EU Grant 604102 
(Human Brain Project, HBP). All network simulations carried out with NEST (http://www.nest-simulator.org).


# References

