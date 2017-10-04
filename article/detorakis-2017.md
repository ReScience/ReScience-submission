---
Title: "A Generalized Linear Integrate-and-Fire Neural Model Produces Diverse Spiking Behaviors"
Author:
  - name: Georgios Detorakis
    affiliation: 1
Address:
  - code:    1
    address: Department of Cognitive Sciences, UC Irvine, Irvine CA, USA
Contact:
  - gdetorak@uci.edu (gdetor@protonmail.com)
Editor:
  - Name Surname
Reviewer:
  - Name Surname
  - Name Surname
Publication:
  received:  Sep,  1, 2015
  accepted:  Sep, 1, 2015
  published: Sep, 1, 2015
  volume:    "**1**"
  issue:     "**1**"
  date:      Sep 2015
  number: 1
Repository:
  article:   "http://github.com/rescience/rescience-submission/article"
  code:      "http://github.com/rescience/rescience-submission/code"
  data:      
  notebook:  
Reproduction:
  - "A Generalized Linear Integrate-and-Fire Neural Model Produces
     Diverse Spiking Behaviors, Stefan Mihalas and Ernst Niebur, Neural
     Computation 21, 704–718, 2009."
Bibliography:
  bibliography.bib

---


# Introduction

Integrate-and-fire neurons are being used extensively in the field of
neuroscience for modeling spiking behaviors [@dayan:2001]. In this work
we provide a reference implementation of [@mihalas:2009], where the
authors have introduced a generalization of the leaky integrate-and-fire
neuron model. The Mihalas-Niebur Neuron (MNN) model is a linear
integrate-and-fire neuron model capable of expressing a rich spiking
behavior based on a set of parameters.

An MNN model expresses tonic and phasic spiking, class $1$ and $2$,
spike frequency adaptation, accommodation, threshold variability,
rebound spike, integrator, input bistability, hyperpolarizing spiking
and bursting, tonic, phasic and rebound bursting, mixed mode,
afterpotentials, basal bistability, preferred frequency and spike
latency. Due to its simplicity, the MNN model has been used in
neuromorphic implementations such as [@folowosele:2011].

The model consists of linear differential equations, which describe the
membrane and threshold potentials and internal currents. All the results
provided in [@mihalas:2009] have been obtained by using only two
internal currents and thus we use the exact same number of internal
currents in this work. The subthreshold dynamics are defined by a set of
linear ordinary differential equations, while an instantaneous threshold
potential controls when the neuron fires an action potential (spike) in
a dynamic way. The ability of the MNN model to generate such a diverse
spiking behavior is due to the complex update rules. In this work the
MNN model has been implemented in Python (version 3.6.1) using Numpy
(version 1.13.1) and Matplotlib (version 2.0.2) packages.


# Methods

In order to implement the model described in [@mihalas:2009], we
discretized the dynamical system using the forward Euler integration
scheme. The time step is fixed to $0.1\, \mathrm{ms}$ for all the
simulations, and the total simulation time $t_f$ varies according to
figure 1 of the original paper. Our implementation differs from the one
in the original paper, since in [@mihalas:2009], authors numerically
solve equation $3.5$ (algebraic equation) under the constraint imposed
by inequality $3.4$ and thus they compute the spike times. On the other
hand, in this work we directly compute numerically the solution of the
dynamical system defined by equations $2.1$ and $2.2$ in [@mihalas:2009]
(see tables @tbl:2 and @tbl:3).

We provide all equations and parameters of the model in tables as it has
been suggested by [@nordlie:2009]. Table @tbl:1 provides the
summary of the model. Tables @tbl:2 and @tbl:3 give the
subthreshold dynamics (differential equations) describing the membrane
and the threshold potentials as well as the two internal currents and
the update rules. The parameters for all the simulations are given in
table @tbl:4, while the external current intensities and pulse
duration are provided in table @tbl:5. The parameters in this work
are exactly the same used in the original paper (table 1, pg. $711$). We
had to infer the time intervals and the total simulation times for the
pulses since they are not given explicitly in the original paper. Thus,
we extracted the time intervals from figure $1$ of [@mihalas:2009] by
visual inspection. The initial conditions are given in
table @tbl:6.

All simulations ran on a Dell OptiPlex $7040$, equipped with a sixth
generation i$7$ processor, $16\, \mathrm{GB}$ of physical memory and
running Arch Linux (x$86\_64$). The total execution time of all
simulations was $2.41$ seconds and the peak consumed memory was
$162\, {\mathrm{MB}}$[^1].



Model Summary                                                  
-------------- ------------------------------------- 
Populations    No population -- single neuron model    
Topology       --                                            
Connectivity   --                                               
Neuron Model   Linear Integrate-and-Fire Neuron        
Channel Models Linear, first order ODEs                
Synapse Model  --                                               
Plasticity     --                                               
Input          Constant current/rectangular pulses     
Measurements   Membrane potential, phase plane         
-------------- ------------------------------------- 
Table: **Summary of the model.** {#tbl:1}


Neuron Model
------------------------------------ ----------------------------------
Name                                 Mihalas-Niebur Neuron (MNN)
Type                                 Linear Leaky Integrate-and-Fire Neuron
Membrane Potential                   $$ \frac{\mathrm{d}V(t)}{\mathrm{d}t} = \frac{1}{C} \Big(I_e + I_1 + I_2 - G(V(t) - E_L)  \Big) $$
Instantaneous Threshold Potential    $$ \frac{\mathrm{d}\Theta(t)}{\mathrm{d}t} = a(V(t) - E_L) - b(\Theta(t) - \Theta_{\infty}) $$
Internal Currents                    $$ \frac{\mathrm{d}I_{1}(t)}{\mathrm{d}t} = -k_1I_1(t) $$
                                     $$ \frac{\mathrm{d}I_{2}(t)}{\mathrm{d}t} = -k_2I_2(t) $$
------------------------------------ ----------------------------------
Table: Description of the subthreshold dynamics of Mihalas--Niebur neuron
model. $V(t)$ and $\Theta(t)$ are the membrane and threshold potentials,
respectively. $E_L$ and $\Theta_{\infty}$ are the reversal potentials for the
membrane and threshold variables, respectively. $a, b, k_1, k_2$ and $G$ are
constant parameters. $I_e$ is the external current applied on the neuron
model. {#tbl:2}

\pagebreak

Variable     Rule                          
------------ ------------------------------
$V(t)$       $V_r$                        
$\Theta(t)$  $\max\{\Theta_r, \Theta(t) \}$  
$I_1(t)$     $R_1 \times I_1(t) + A_1$      
$I_2(t)$     $R_2 \times I_2(t) + A_2$     
------------ ------------------------------
Table: **Update rules.** $V_r$ and $\Theta_r$ are the reset values for the 
membrane and threshold potentials, respectively. $R_1, R_2, A_1$ and $A_2$ are
constants. {#tbl:3}


Figure           $$a (\mathrm{s^{-1}})$$ $$\frac{A_1}{C} (\mathrm{V/s})$$ $$\frac{A_2}{C} (\mathrm{V/s})$$ $$t_f \mathrm{s}$$ 
---------------- ----------------------- -------------------------------- -------------------------------- ------------------
1A                 $0$                    $0$                               $0$                             $0.2$
1B                 $0$                    $0$                               $0$                             $0.5$
1C                 $5$                    $0$                               $0$                             $0.2$
1D                 $5$                    $0$                               $0$                             $0.5$
1E                 $5$                    $0$                               $0$                             $1.0$
1F                 $5$                    $0$                               $0$                             $0.4$
1G                 $5$                    $0$                               $0$                             $1.0$
1H                 $5$                    $0$                               $0$                             $0.3$
1I                 $5$                    $0$                               $0$                             $0.4$
1J                 $5$                    $0$                               $0$                             $1.0$
1K                $30$                    $0$                               $0$                             $0.4$
1L                $30$                   $10$                             $-0.6$                            $0.4$
1M                 $5$                   $10$                             $-0.6$                            $0.5$
1N                 $5$                   $10$                             $-0.6$                            $0.5$
1O                 $5$                   $10$                             $-0.6$                            $1.0$
1P                 $5$                    $5$                             $-0.3$                            $0.5$
1Q                 $5$                    $5$                             $-0.3$                            $0.2$
1R                 $0$                    $8$                             $-0.1$                            $0.2$
1S                 $5$                   $-3$                              $0.5$                            $0.8$
1T               $-80$                    $0$                                $0$                           $0.05$
---------------- ----------------------- -------------------------------- -------------------------------- ------------------
Table: **Simulation Parameters**. Common parameters for all simulations: $b = 10\, \mathrm{s^{-1}}$, $G/C = 50\, \mathrm{s^{-1}}$, $k_1 = 200\, \mathrm{s^{-1}}$, $k_2 = 20\, \mathrm{s^{-1}}$, $\Theta_{\infty} = -0.05\, \mathrm{V}$, $R_1 = 0$, $R_2 = 1$, $E_l = -0.07\, \mathrm{V}$, $V_r = -0.07\, \mathrm{V}$, $\Theta_r = -0.06\, \mathrm{V}$. {#tbl:4}

\pagebreak


Figure Type                                         $$I_e/C (V/s)$$
------ -------------------------------------------- ------------------------
1A     ![image](figs/const.pdf){width="10.00000%"}  $1.5$
1B     ![image](figs/const.pdf){width="10.00000%"}  $1 + 10^{-6}$
1C     ![image](figs/const.pdf){width="10.00000%"}  $2$
1D     ![image](figs/const.pdf){width="10.00000%"}  $1.5$
1E     ![image](figs/pulse.pdf){width="10.00000%"}  $$1.5 (0.1 \mathrm{s}),\, 0 (0.5\mathrm{s}),\, 0.5 (0.1\mathrm{s}),\, 1 (0.1\mathrm{s}),\, 1.5 (0.1\mathrm{s}),\, 0 (0.1 \mathrm{s}) $$
1F     ![image](figs/pulse.pdf){width="10.00000%"}  $$1.5 (0.02 \mathrm{s}),\, 0 (0.18 \mathrm{s}),\, -1.5 (0.025 \mathrm{s}),\, 0 (0.025 \mathrm{s}),\, 1.5 (0.025 \mathrm{s}),\, 0 (0.125 \mathrm{s})$$
1G     ![image](figs/pulse.pdf){width="10.00000%"}  $$0 (0.05 \mathrm{s}),\, -3.5 (0.756 \mathrm{s}),\, 0 (0.194 \mathrm{s}) $$ 
1H     ![image](figs/const.pdf){width="10.00000%"}  $$ 2(1 + 10^{-6}) $$
1I     ![image](figs/pulse.pdf){width="10.00000%"}  $$ 1.5 (0.02 \mathrm{s}),\, 0 (0.01 \mathrm{s}),\, 1.5 (0.02 \mathrm{s}),\, 0 (0.25 \mathrm{s}),\, 1.5 (0.02 \mathrm{s}),\, 0 (0.02 \mathrm{s})$$ 
                                                    $$ 1.5 (0.02 \mathrm{s}),\, 0 (0.04 \mathrm{s}) $$
1J     ![image](figs/pulse.pdf){width="10.00000%"}  $$ 1.5 (0.1 \mathrm{s}),\, 1.7 (0.4 \mathrm{s}),\, 1.5 (0.1 \mathrm{s}),\, 1.7 (0.4 \mathrm{s}) $$
1K     ![image](figs/const.pdf){width="10.00000%"}  $-1$
1L     ![image](figs/const.pdf){width="10.00000%"}  $-1$
1M     ![image](figs/const.pdf){width="10.00000%"}  $2$
1N     ![image](figs/const.pdf){width="10.00000%"}  $1.5$
1O     ![image](figs/pulse.pdf){width="10.00000%"}  $$0 (0.1 \mathrm{s}),\, -3.5 (0.5 \mathrm{s}),\, 0 (0.4 \mathrm{s})$$
1P     ![image](figs/const.pdf){width="10.00000%"}  $2$
1Q     ![image](figs/pulse.pdf){width="10.00000%"}  $$2 (0.015 \mathrm{s}),\, 0 (0.185 \mathrm{s})$$
1R     ![image](figs/pulse.pdf){width="10.00000%"}  $$5 (0.01 \mathrm{s}),\, 0 (0.09 \mathrm{s}),\, 5 (0.01 \mathrm{s}),\, 0 (0.09 \mathrm{s})$$
1S     ![image](figs/pulse.pdf){width="10.00000%"}  $$5 (0.005 \mathrm{s}),\, 0 (0.005 \mathrm{s}),\, 4 (0.005 \mathrm{s}),\, 0 (0.385 \mathrm{s}),\, 5 (0.005 \mathrm{s}),\, 0 (0.045 \mathrm{s}) $$
                                                    $$4 (0.005\mathrm{s}),\, 0 (0.345 \mathrm{s}) $$
1T     ![image](figs/pulse.pdf){width="10.00000%"}  $$8 (0.002 {\mathrm{s}}),\, 0 (0.048 \mathrm{s})$$
------ -------------------------------------------- ------------------------
Table: **External current**. This table provides the external current
for each panel in Figure @fig:1. There are two types of external
currents, constants and pulses. In the case of pulses the duration of
each pulse is given in seconds along with its intensity. {#tbl:5}


Variable           Initial Value
------------------ -------------
$$V(t)$$           $$-0.07\, \mathrm{V}$$ / $$-0.03\, \mathrm{V}$$ (Figure 1H)
$$\Theta(t)$$      $$-0.05\, \mathrm{V}$$ / $$-0.03\, \mathrm{V}$$ (Figure 1H)
$$I_1(t)$$         $$0.01\, \mathrm{V}$$
$$I_2(t)$$         $$0.001\, \mathrm{V}$$
------------------ ----------------------
Table: **Initial Conditions** In all simulations the same initial conditions have
been used, except from the one illustrated in Figure @fig:1 H. {#tbl:6}


# Results

All three figures from the original article have been successfully
replicated. All the different spiking behaviors of the model are
illustrated in Figure @fig:1, where the black solid line indicates
the membrane potential ($V(t)$), the red dashed line illustrates the
instantaneous threshold potentials ($\Theta(t)$), and the gray line
shows the input to the neuron ($I_e/C$). The $x$-axis scales in all
panels are exactly the same as in the original paper (indicating the
total simulation time ($t_f$), while the $y$-axis scale differs from the
one in the original paper. In this work the $y$-axis scale is the same
for all the subplots ($[-95, -25]\, {\mathrm{mV}}$), except for
panels G and O ($[-145, -25]\, {\mathrm{mV}}$).

Figures @fig:2 and @fig:3 depict the phase space of the phasic
spiking ($V(t)$ and $\Theta(t)$) and phasic bursting ($V(t)$, $I_1(t)$,
and $I_2(t)$). In both figures the blue curves and the black dots
indicate the trajectory of the system and spike events, respectively. In
Figure @fig:2 the gray arrows show the evolution of the system
(vector field of the system). Figure @fig:3 has a different
orientation from the original one but both the original figure and the
replicated one illustrate the same trajectories and spike events of the system.
All the figures express the same qualitative behavior as the original figures
in [@mihalas:2009].

![**Neural responses of MNN.** Black solid lines indicate the
membrane potential ($V(t)$), the red dashed lines show the threshold
potentials ($\Theta(t)$), and the gray lines the external currents
applied on each case. [**A**]{} tonic spiking, [**B**]{} class $1$,
[**C**]{} spike frequency adaptation, [**D**]{} phasic spiking,
[**E**]{} accommodation, [**F**]{} threshold variability, [**G**]{}
rebound spike, [**H**]{} class $2$, [**I**]{} integrator, [**J**]{}
input bistability, [**K**]{} hyperpolarization induced spiking,
[**L**]{} hyperpolarization induced bursting, [**M**]{} tonic bursting,
[**N**]{} phasic bursting, [**O**]{} rebound burst, [**P**]{} mixed
mode, [**Q**]{} afterpotentials, [**R**]{} basal bistability, [**S**]{}
preferred frequency, [**T**]{} spike
latency.](figs/Figure01.pdf){#fig:1}

![**Phase space of phasic spiking.** Blue solid lines indicate the
trajectories of the model in the phase spiking behavior
(Figure @fig:1 D). The dashed line corresponds to $V(t) =
    \Theta(t)$, and the black dots represent spike events. The
parameters for this simulation are the same as in
Figure @fig:1 D.](figs/Figure02.pdf){#fig:2}

![**Phase space of phasic bursting.** Blue solid lines represent the
trajectories of the system and the black dots indicate spiking events.
The parameters for this simulation are the same as in
Figure @fig:1 N.](figs/Figure03.pdf){#fig:3}


# Conclusion

All figures in @mihalas:2009 have been successfully replicated with high
fidelity. Overall, the whole reproducing process was smooth and without
obscure points since most of the parameters are provided in the original
article. Only the time intervals for which the external current is
applied to the model and the initial conditions are not provided
explicitly. Therefore, we had to extract that information from figure
$1$ of the original article. To conclude, the article [@mihalas:2009]
has been successfully reproduced without any discrepancy.


[^1]: Python memory profiler used
    (<https://pypi.python.org/pypi/memory_profiler>).


# References

