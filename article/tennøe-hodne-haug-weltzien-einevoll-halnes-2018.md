---
Title: "Fast-Activating Voltage- and Calcium-Dependent Potassium (BK) Conductance
        Promotes Bursting in Pituitary Cells: A Dynamic Clamp Study"
Author:
  - name: Simen Tennøe
    affiliation: 1, 2
  - name: Kjetil Hodne
    affiliation: 3
  - name: Trude M. Haug
    affiliation: 4
  - name: Finn-Arne Weltzien
    affiliation: 3
  - name: Gaute T. Einevoll
    affiliation: 1, 5, 6
  - name: Geir Halnes
    affiliation: 1, 5
Address:
  - code:    1
    address: Centre for Integrative Neuroplasticity, University of Oslo, Oslo, Norway
  - code:    2
    address: Department of Informatics, University of Oslo, Oslo, Norway
  - code:    3
    address: Department of Basic Sciences and Aquatic Medicine, Norwegian University of Life Sciences, Campus Adamstuen, Norway
  - code:    4
    address: Institute of Oral Biology, University of Oslo, Oslo, Norway
  - code:    5
    address: Faculty of Science and Technology, Norwegian University of Life Sciences, Ås, Norway
  - code:    6
    address: Department of Physics, University of Oslo, Oslo, Norway
Contact:
  - simenten@student.matnat.uio.no
Editor:
  - Nicolas P. Rougier
Reviewer:
  - Name Surname
  - Name Surname
Publication:
  received:  Nov, 8, 2018
  accepted:  Nov, 1, 2018
  published: Nov, 1, 2018
  volume:    "**4**"
  issue:     "**1**"
  date:      Nov 2018
  number: 1
Repository:
  article:   "http://github.com/rescience/rescience-submission/article"
  code:      "http://github.com/rescience/rescience-submission/code"
  data:
  notebook:
Reproduction:
  - "*Fast-Activating Voltage- and Calcium-Dependent Potassium (BK) Conductance
     Promotes Bursting in Pituitary Cells: A Dynamic Clamp Study*,
     J. Tabak, M. Tomaiuolo, A. Gonzalez-Iglesias,  L. Milescu and R. Bertram,
     Journal of Neuroscience 31.46 (2011), 10.1523/JNEUROSCI.3235-11.2011"
Bibliography:
  bibliography.bib

---

# Introduction

As part of the dynamic clamp study by Tabak et al. 2011 [@tabak2011],
a computational model was developed for the voltage dynamics of endocrine
pituitary cells in rats.
The model captured the spontaneous activity of these cells,
including the generation of Ca\textsuperscript{2+}-channel mediated spikes and pseudo-plateau bursts.
As an important achievement,
the model explained the paradoxical role that big conductance
K\textsuperscript{+} (BK) channels had in prolonging spike duration
and sometimes promoting burst firing in these cells [@vangoor2001],
contrary to what one would expect from a hyperpolarizing current.
The original model was implemented in XPP [@ermentrout2002].
The code for the model was made available online at
[https://www.math.fsu.edu/~bertram/software/pituitary/JNS_11b.ode](https://www.math.fsu.edu/~bertram/software/pituitary/JNS_11b.ode),
 while the code used in the analysis of the model outcome was not made available.

In the current paper, we have reimplemented the computational model by
Tabak et al. [@tabak2011] using the Python interface for the NEURON simulator [@hines1997],
a widely used simulator for multicompartmental neurons.
In addition, we have performed an uncertainty quantification and sensitivity
analysis of the model using the Uncertainpy Python package [@uncertainpy],
version 1.1.4 (Zenodo: [10.5281/zenodo.1473453](http://doi.org/10.5281/zenodo.1473453)).
The model implementation works with Python 2 and 3.
The results in this paper were created using Python 3.7.0 within a
Docker ([https://www.docker.com/](https://www.docker.com/)) environment.

The reimplemented model reproduced the characteristic firing patterns seen
in the original publication,
and we thus confirmed the original study.
The sensitivity analysis further presented a systematic overview of the model in
terms of how its characteristic response features depended on the various model
parameters.
Supporting the main conclusion from the original work,
the sensitivity analysis showed that the bursting propensity of the model was
highly sensitive to the BK conductance.
However, the analysis also revealed that the bursting propensity was sensitive
to additional parameters (conductances),
and thus that BK is not the sole determinant for whether the cell is bursty.


# Methods

When reimplementing the model by Tabak et al. [@tabak2011] we followed the descriptions in the original
publication,
using the original implementation for verification purposes.
We also had a brief communication with the original authors to obtain details
on the analysis part of the model.


## Model

The model by Tabak et al. [@tabak2011] was defined by the equation:

$$C \frac{dV}{dt} = - (I_{\mathrm{Ca}} + I_{\mathrm{K}} + I_{\mathrm{BK}} + I_{\mathrm{SK}} + I_{\mathrm{leak}} + I_{\mathrm{noise}} ),$${#eq:V}

\noindent
where $C$ is the membrane capacitance,
$V$ is the membrane potential,
and $I_X$ the current through a specific ion channel $X$.
The model included six different currents:

* $I_{\mathrm{Ca}}$ -- Voltage gated Ca\textsuperscript{2+} current.
* $I_{\mathrm{K}}$ -- Voltage gated K\textsuperscript{+} current.
* $I_{\mathrm{BK}}$ -- Big conductance K\textsuperscript{+} current.
* $I_{\mathrm{SK}}$ -- Small conductance K\textsuperscript{+} current.
* $I_{\mathrm{leak}}$ -- Leak current.
* $I_{\mathrm{noise}}$ -- Stochastic current representing channel noise.

\noindent
A current through an ion channel $X$ was given by the simplified relation:

$$I_{X} = G_{X}Y_X(V - E_{X}),$${#eq:I}

\noindent
where $G_{X}$ denotes the maximum ion channel conductance,
and $E_{X}$ denotes the reversal potential of the ion species conducted by channel $X$.
$Y_X$ denotes an ion channel specific gating function,
which was unity for $I_{\mathrm{leak}}$,
an instantaneous function of $V$ for $I_{\mathrm{Ca}}$,
an instantaneous function of the cytosolic Ca\textsuperscript{2+} concentration
for $I_{\mathrm{SK}}$,
and a dynamic function of $V$ and $t$ for the remaining ion channels
$I_{\mathrm{K}}$, $I_{\mathrm{BK}}$, and $I_{\mathrm{K}}$.

The original implementation used the total membrane capacitance (units pF)
and total membrane conductances (units nS),
while the NEURON simulator requires these entities to be specified per membrane
area with units $\mu \mathrm{F/cm}^2$ and S/cm\textsuperscript{2},
respectively.
NEURON also requires that the membrane area is defined.
To get the parameters on the form required by NEURON we defined an arbitrary
membrane area ($A$),
and divided the capacitance $C$ and ion channel conductances $G_X$ by $A$:

$$g_{X,\,\mathrm{NEURON}} = \frac{G_X}{A}, \qquad c_{\mathrm{NEURON}} = \frac{C}{A}.$${#eq:c_g}

\noindent
Combining Equation @eq:V, @eq:I and @eq:c_g shows that the model is independent
of the choice of $A$:

$$\frac{C}{A} \frac{dV}{dt} =  -\left(\frac{G_{X}}{A}Y_X(V - E_{X}) + \ldots \right).$${#eq:V_A}

The original model further included an equation for handling the intracellular
Ca\textsuperscript{2+} concentration,
which is relevant for the gating of SK channels:

$$\frac{d[Ca]}{dt} = - f_c (\alpha I_{\mathrm{Ca}} + k_c[Ca]), $${#eq:Ca}

\noindent
where $f_c$  denotes the fraction of free Ca\textsuperscript{2+} in the cytoplasm,
$k_c$ denotes the extrusion rate,
and the constant $\alpha$ converts an incoming current to a molar concentration.
$\alpha$ was converted to NEURON units by taking:

$$\alpha_{\mathrm{NEURON}} = A\alpha.$${#eq:alpha}

\noindent
Combining Equation @eq:c_g, @eq:Ca and @eq:alpha shows that this choice keeps
the model independent of the choice of $A$:

$$\frac{d[Ca]}{dt} = - f_c (A\alpha \frac{G_{\mathrm{Ca}}}{A}(V - E_{\mathrm{Ca}}) + k_c[Ca]).$${#eq:Ca_A}

We arbitrarily chose a cell body with a membrane area of $\pi \cdot 10^{-6}$ cm\textsuperscript{2},
i.e. with a diameter of 10 $\mu$m.
We used all equations from the original publication,
substituting $G_X$ with $g_{X,\,\mathrm{NEURON}}$,
$C$ with $c_{\mathrm{NEURON}}$, and $\alpha$ with $\alpha_{\mathrm{NEURON}}$.
The parameter values from the original publication and the converted parameter
values are summarized in Table @tbl:parameters.
Parameters not listed in this table were kept unchanged from the original publication.
To make the discussion and results easier to compare to the original publication,
we will refer to the original conductance values through the rest of this paper.

The noise was added by using a current clamp that injected a random current at
each time step in the simulation,
as described by the original publication.
Simulations with noise were run with a fixed time step of `dt = 0.01` ms,
which is the same time step used in the original publication.
When performing the sensitivity analysis,
the noise amplitude was set to zero ($A_{\mathrm{noise}} = 0$),
and the simulations were run using adaptive time steps.

We found one discrepancy between the parameters listed in the original
publication and the values found in the original source code.
The maximum conductance of K\textsuperscript{+} channels ($G_{\mathrm{K}}$)
was listed as 3.2 nS in the original publication,
while the value used in the original source code was 3 nS.
Both values were tested and $G_{\mathrm{K}} = 3$ nS gave results most similar
to the results in the original publication.
We therefore decided to use $G_{\mathrm{K}} = 3$ nS instead of the value listed
in the original publication.


Tabak                 Value       Unit                     NEURON             Value                 Unit
--------------------- ----------- -----------------------  ------------------ --------------------- -------------------------------------------------------------
                                                           `A`                $3.14 \cdot 10^{-6}$  cm\textsuperscript{2}
C                     10          pF                       `c`                1.6                   $\mu \mathrm{F/cm}^2$
$G_{\mathrm{Ca}}$     2           nS                       `g_Ca`             $6.37 \cdot 10^{-4}$  S/cm\textsuperscript{2}
$G_{\mathrm{K}}$      3           nS                       `g_K`              $9.55 \cdot 10^{-4}$  S/cm\textsuperscript{2}
$G_{\mathrm{BK}}$     2           nS                       `g_BK`             $0$                   S/cm\textsuperscript{2}
$G_{\mathrm{SK}}$     2           nS                       `g_SK`             $6.37 \cdot 10^{-4}$  S/cm\textsuperscript{2}
$G_{\mathrm{l}}$      0.2         nS                       `g_l`              $6.37 \cdot 10^{-5}$  S/cm\textsuperscript{2}
$\alpha$              0.0015      $\mu \mathrm{M/fC}^2$    `alpha`            $4.71 \cdot 10^{-3}$  $\mathrm{mM} \cdot \mathrm{cm}^2 \mathrm{/} \mu \mathrm{C}$
$A_{\mathrm{noise}}$  4           pA                       `noise_amplitude`  $0.004$               nA
--------------------- ----------- -----------------------  ------------------ --------------------- -------------------------------------------------------------

Table: The parameter values in Tabak et al. [@tabak2011] that were converted from currents and capacitance to currents and capacitance per membrane area due to requirements by the NEURON simulator. The original model parameter values are denoted Tabak while the parameter values in the reimplemented model are denoted NEURON, with names as used in the model implementation. {#tbl:parameters}


## Event detection

In the analysis,
we ran the model for 60000 ms
and discarded the first 10000 ms of the voltage trace to eliminate the
transient initial response.

The first step of the model analysis was to detect events (spikes or bursts)
in the model voltage trace.
To do this, the voltage was normalized so that its minimum value was set to 0
and its maximum was set to 1.
The start of an event was specified to be when the voltage crossed an onset
threshold (defined to be $0.55$),
and the end of an event to be when it next descended below another,
lower termination threshold (defined to be $0.45$).
An event includes the first point before it crossed the onset threshold and the
first point after it descended below the termination threshold.

The difference in onset and termination threshold was necessary to prevent
random fluctuations around the threshold (during upstroke or downstroke)
to be considered as independent events.
If the voltage trace started above the onset threshold,
we discarded the first part of the voltage trace until we got below the
termination threshold.
Similarly,
if an event did not fall below the termination threshold before the simulation
ended, that event was discarded.
Additionally,
we required that events have an amplitude of at least 10 mV.
This prevents the problem where the normalization step leads to detecting false
events with an amplitude less than 1 mV in cases where the model does not generate
any events and instead exhibits small (much less than 1 mV) fluctuations around
a steady state.

We used Uncertainpy to detect events, as the described threshold-detection
algorithm is available to us by using the `uncertainpy.Spikes` object with the
arguments `normalize=True`, `trim=False` and `min_amplitude = 10`.
Note that in Uncertainpy the `end_threshold` is given relative to the onset
threshold,
so to get a termination threshold = 0.45 we set `end_threshold = -0.1`.

An event was defined as a burst when its duration was longer than
a given threshold (60 ms).
The burstiness factor was defined as the fraction of the total number of events
that were considered as bursts.
All parameters used in the analysis are summarized in Table @tbl:parameters_analysis.

The description of the threshold-detection algorithm for detecting events
(bursts or spikes) was incomplete in the original publication.
We contacted the original authors,
who were helpful in describing the threshold-detection algorithm,
but who did not recall the exact numerical values of all threshold choices.
The onset threshold,
termination threshold and burst-duration threshold used
(Table @tbl:parameters_analysis) were therefore set to the values we found to
give the best agreement between our analysis outcome and
that in the original publication.

Parameter                    Value       Unit
---------------------------- ----------- -----------------------
Simulation time              $60000$     ms
Discard                      $10000$     ms
Event onset threshold        $0.55$
Event termination threshold  $0.45$
Burst threshold              $60$        ms
Minimum event amplitude      $10$        mV
---------------------------- ----------- -----------------------

Table: The parameters used in the analysis of the model. {#tbl:parameters_analysis}



## Uncertainty quantification and sensitivity analysis

We used Uncertainpy to further examine the model through an uncertainty
quantification and sensitivity analysis.
This enabled us to quantify how sensitive salient response properties of the
model is to changes in the various parameters.
In the sensitivity analysis,
the four conductances $G_{\mathrm{Ca}}$, $G_{\mathrm{K}}$, $G_{\mathrm{SK}}$,
and $G_{\mathrm{l}}$ were assigned uniform distributions within $\pm 50\%$
of their original values. $G_{\mathrm{BK}}$,
which had no default value in the original model,
was given a uniform distribution between 0 and 1 nS as this was the parameter
range explored in the original study.
We use polynomial chaos with the point collocation method
(the default of Uncertainpy) and a polynomial order of eight.
In the sensitivity analysis,
we wanted all the variance in the simulation
outcome to reflect parameter variations,
and the random noise was therefore turned off by setting $A_{\mathrm{noise}} = 0$.

We calculated the uncertainty and sensitivity of the five features of the model:

* Event rate, which is the event firing rate (named `spike_rate` in Uncertainpy).
* Average event peak, which is the average event peak voltage (named \newline  `average_AP_overshoot` in Uncertainpy).
* Average AHP (afterhyperpolarization) depth, which is the average minimum voltage between events.
* Burstiness factor, the fraction of events with a duration longer than 60 ms.
* Average duration, the average duration of an event.

Some of these features are not defined for all parameter combinations
(for example average AHP depth is not defined when there are no events).
The point collocation method still gives reliable results,
as long as the features are defined for a sufficiently large fraction of the
parameter combinations (in our case the lowest was ~91.5\%) [@eck2016].

Some of the outcomes from the sensitivity analysis were unexpected
and were explored further by varying selected parameters and documenting
how these variations affected the average event duration and burstiness factor
of the model.
The parameters varied in this additional analysis were $G_{\mathrm{BK}}$,
$G_{\mathrm{SK}}$, and $G_{\mathrm{K}}$,
all of which were varied within the range used in the uncertainty analysis.


# Results

We repeated all simulations and qualitatively reproduced Figure 1 and Figure 2
in Tabak et al. [@tabak2011].
The remaining results in the original publication were experimental results,
and therefore outside the scope of this reproduction.

The results shown in Figure @fig:figure1 correspond well to those in
Figure 1 of the original publication.
The original model and the reproduced version showed the same behavior when
increasing $G_\mathrm{BK}$.
As random noise was added to the simulations,
an exact replication could not be expected.

![Model predictions for the effect of various $G_\mathrm{BK}$ conductances on burstiness. \textbf{A}-\textbf{C} Left, membrane potential of the model. Right, distribution of event durations in the time interval from 1 to 5 s (of the 50 s simulated). The grey line indicates the threshold for what is considered a spike and what is considered a burst, and BF denotes the burstiness factor. \textbf{D} The burstiness factor increased with $G_\mathrm{BK}$. \textbf{E} The burstiness factor decreased with $\tau_\mathrm{BK}$.](figures/figure_1.eps){#fig:figure1}

The results shown in Figure @fig:figure2 correspond well to those in Figure 2
of the original publication.
The behavior we observe is similar to the behavior in the original publication.
When increasing $G_\mathrm{BK}$ the model went from having a low burstiness factor
(between 0 and 0.1) to having a high burstiness factor (between 0.9 and 1).
For any value of $G_\mathrm{BK}$,
the model evaluations tended to be either predominantly bursting
(i.e. most events were bursts) or predominantly spiking
(few events were bursts),
so that the number of model evaluations with an intermediate burstiness factor
between 0.1 and 0.9 (events changed between being bursts or not)
was always low (less than 20 evaluations).

A small deviation was found between the current analysis and the original work.
For the maximum value of $G_\mathrm{BK}$,
we got fewer model evaluations with low burstiness
than in the original work.
We do not know the precise cause of this difference,
though we can speculate that it may be due to smaller differences in the
performed analysis (e.g., onset and termination threshold definitions or implementations),
or due to underlying differences between the NEURON and XPP implementations
(e.g., NEURON uses backward Euler as the numerical integration scheme
while the XPP implementation uses forward Euler).


![Robustness of the burstiness of the model for three values of $G_{\mathrm{BK}}$ when changing $G_{\mathrm{Ca}}$, $G_{\mathrm{K}}$, $G_{\mathrm{SK}}$, and $G_{\mathrm{l}}$ uniformly within $\pm 50\%$ of their original values. \textbf{A} For $G_{\mathrm{BK}} \rightarrow 0$ nS, $67.5\%$ of the active models were spikers (burstiness factor $< 0.3$). \textbf{B} For $G_{\mathrm{BK}} \rightarrow 0.5$ nS, $33.8\%$ were spikers. \textbf{C} For $G_{\mathrm{BK}} \rightarrow 1$ nS, only $4.4\%$ were spikers.](figures/figure_2.eps){#fig:figure2}




## Uncertainty quantification and sensitivity analysis

The uncertainty quantification and sensitivity analysis of the model is shown
in Figure @fig:sensitivity.
The sensitivity was given as the total-order Sobol indices,
which quantify how much of the variance of the model each parameter
(accounting for all of its interactions with other parameters)
is responsible for [@homma1996].

![Uncertainty quantification and sensitivity analysis of a selected set of response features of the model. \textbf{A} Event rate denotes the event firing rate. \textbf{B} Average event peak denotes is the average event peak voltage. \textbf{C} Average AHP (afterhyperpolarization) depth denotes the average minimum voltage between two consecutive events. \textbf{D} Burstiness factor denotes the fraction of events with duration longer than 60 ms. \textbf{E} Average duration denotes the average duration of the events.](figures/sensitivity.eps){#fig:sensitivity}

The sensitivity analysis showed that the spike rate was sensitive to almost all
ion channel conductances,
but most so to $G_{\mathrm{K}}$ (Figure @fig:sensitivity\textbf{A}).
Such a role of the delayed rectifying K\textsuperscript{+} channel in
controlling the firing rate has been seen in other studies [@guan2013].

The event amplitude was mainly sensitive to $G_{\mathrm{Ca}}$
(Figure @fig:sensitivity\textbf{B}),
which is not surprising given that the events are generated by
$I_{\mathrm{Ca}}$.
However, it also had a relatively high sensitivity to $G_{\mathrm{BK}}$,
in line with what was found in the previous study [@tabak2011].

The average afterhyperpolarization depth was in turn
most sensitive to $G_{\mathrm{Ca}}$ (Figure @fig:sensitivity\textbf{C}).
This may seem counterintuitive,
as $I_{\mathrm{Ca}}$ is not a hyperpolarizing current.
However,
$I_{\mathrm{Ca}}$ is responsible for triggering all the three
hyperpolarizing currents ($I_\mathrm{K}$, $I_\mathrm{BK}$ and $I_\mathrm{SK}$)
that generate the afterhyperpolarization depth.
$I_\mathrm{K}$ and $I_\mathrm{BK}$ are activated by the voltage
deflection caused by $I_{\mathrm{Ca}}$,
while $I_\mathrm{SK}$ is activated by the Ca\textsuperscript{2+}
entering through $I_{\mathrm{Ca}}$.


The burstiness factor of the model was mainly sensitive to $G_{\mathrm{K}}$
and $G_{\mathrm{BK}}$ (Figure @fig:sensitivity\textbf{D}).
The sensitivity to $G_{\mathrm{BK}}$ confirms the findings in the original publication,
i.e. that BK channels promote bursting.
However, the large sensitivity to $G_{\mathrm{K}}$ is a novel insight for the current study
and indicates that also $G_{\mathrm{K}}$ was important for determining if the
model produced bursts or spikes.
This observation is tightly related to the explanation for how BK can act as a
burst promoter in the first place,
which is contrary to what one would expect from a hyperpolarizing current.
The explanation,
proposed by both Tabak et al. [@tabak2011] and the experimental
studies they were inspired by [@vangoor2001],
was that $G_{\mathrm{BK}}$ promoted bursting by reducing the peak amplitude of
events (as reflected in Figure @fig:sensitivity\textbf{B}),
thereby preventing full activation of the otherwise more strongly
hyperpolarizing delayed rectifier current ($I_{\mathrm{K}}$).
In this context, the sensitivity analysis simply shows that the indirect effect on
$I_{\mathrm{K}}$ obtained by varying $G_{\mathrm{BK}}$ was smaller than the
direct effect on $I_{\mathrm{K}}$ obtained by varying $G_{\mathrm{K}}$
(Figure @fig:sensitivity\textbf{D}).

Surprisingly,
the average event duration had a very low sensitivity to
$G_{\mathrm{BK}}$ (Figure @fig:sensitivity\textbf{E}),
and was instead most sensitive to $G_{\mathrm{SK}}$.
This was unexpected
since the burstiness was highly sensitive to $G_{\mathrm{BK}}$,
and a burst was defined as an event exceeding a certain duration.
An exploration of the counterintuitive relationship between
Figure @fig:sensitivity\textbf{D} and \textbf{E} is presented below.

## Parameter exploration

To explore the relationship between the results in
Figure @fig:sensitivity\textbf{D} and \textbf{E},
we examined the effects of varying $G_\mathrm{BK}$, $G_{\mathrm{K}}$,
and $G_{\mathrm{SK}}$ on the burstiness
and the average duration of events (Figure @fig:durations).
It should be noted that this figure only shows how the model responds when
changing two parameters at the time,
so the higher-order interactions included in the total-order Sobol
sensitivity indices are absent.

![The average duration of events while varying $G_\mathrm{BK}$ and either \textbf{A} $G_{\mathrm{K}}$ or \textbf{B} $G_{\mathrm{SK}}$. The areas in parameter space where the average duration of the events is longer than the burstiness factor threshold are in green, while the areas where the average duration is below this threshold are in yellow. Areas in blue produce no events and the average duration is then set to -1 for visualization purposes.](figures/durations.eps){#fig:durations}

Figure @fig:durations\textbf{A} shows the regions in the
$G_{\mathrm{BK}}$/$G_{\mathrm{K}}$ parameter plane where the model produced
regular spikes (yellow) and bursts (green).
For low ($<$ 2 nS) values of $G_{\mathrm{K}}$,
the cell was bursting regardless of the value of $G_{\mathrm{BK}}$.
Hence, for low values of $G_{\mathrm{K}}$,
the burstiness of the cell was insensitive to $G_{\mathrm{BK}}$.
In comparison, a sufficiently large change in $G_{\mathrm{K}}$
could switch the cell between a regular and bursty state for any (fixed) value of $G_{\mathrm{BK}}$.
These results thus fit well with the sensitivity analysis in Figure @fig:sensitivity\textbf{D},
which showed that the burstiness was more sensitive to $G_{\mathrm{K}}$ than to $G_{\mathrm{BK}}$.

We next fixed $G_{\mathrm{K}}$ at the default value 3 nS,
and explored how it could be the case that burstiness was sensitive to $G_{\mathrm{BK}}$
but not so much to $G_{\mathrm{SK}}$,
while event duration was sensitive to $G_{\mathrm{SK}}$ but not so much to
$G_{\mathrm{BK}}$ (Figure @fig:durations\textbf{B}).
As the figure shows,
for $G_{\mathrm{BK}}<$ 0.2 nS the cell was always regularly spiking,
while for $G_{\mathrm{BK}}>$ 0.8 nS, the
cell was always bursting, irregardless of the values of $G_{\mathrm{SK}}$.
In comparison,
changing $G_\mathrm{BK}$ (keeping $G_{\mathrm{SK}}$ fixed) could always
switch the cell between a regular and bursty state.
In equivalence with the analysis of Figure @fig:durations\textbf{A},
this explains why the burstiness was less sensitive to $G_{\mathrm{SK}}$
than $G_{\mathrm{BK}}$.
However,
although changes in $G_{\mathrm{BK}}$ more often led to changes in burstiness,
the effects on the event duration was modest.
That is,
for most (fixed) values of $G_{\mathrm{SK}}>$,
a change in $G_{\mathrm{BK}}$ could push the event duration from slightly below to
slightly above the burst-duration threshold,
but did not lead to larger changes in burstiness.
Oppositely,
reducing $G_\mathrm{SK}$ to the lower values in the explored range resulted
in burst durations of several thousands of milliseconds
(as long as $G_{\mathrm{BK}}>$ 0.2 nS).
Hence,
while $G_{\mathrm{BK}}$ was important for achieving a burst in the first place,
$G_{\mathrm{SK}}$ had a much larger impact on the duration of the burst.
This explains the difference in sensitivity between the average duration and
burstiness factor observed in Figure @fig:sensitivity\textbf{D} and \textbf{E}.


# Conclusion

<!-- Conclusion, at the very minimum, should indicate very clearly if you were able
to replicate original results. If it was not possible but you found the reason
why (error in the original results), you should explain it. -->

We were able to qualitatively reproduce all the computational results in
Tabak et al. [@tabak2011].
By performing an uncertainty quantification and sensitivity analysis we confirmed
the key conclusions in the original publication using a different simulator and
different analysis methods, which provided additional insight into how different
membrane mechanisms interact
to produce the characteristic response features of the model.
Overall,
the reproduction effort went smoothly, with a little help from the original authors
in describing the threshold-detection algorithm used in the analysis of the model.
The original model now exists as a model using the Python interface for NEURON,
which hopefully makes it accessible to a wider audience.

# References
