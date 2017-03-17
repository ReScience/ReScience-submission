---
Title: "Modeling GABA Alterations in Schizophrenia: A Link Between Impaired Inhibition and Gamma and Beta Auditory Entrainment"
Author:
  - name: Metzner Christoph
    affiliation: 1
Address:
  - code:    1
    address: Biocomputation Research Group, Science and Technology Research Institute,University of Hatfield, Hatfield, Hertfordshire, UK

Contact:
  - c.metzner@herts.ac.uk
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
Repository:
  article:   "http://github.com/rescience/rescience-submission/article"
  code:      "http://github.com/rescience/rescience-submission/code"
  data:      
  notebook:  
Reproduction:
  - "Modeling GABA Alterations in Schizophrenia: A Link Between Impaired Inhibition and Gamma and Beta Auditory Entrainment, D. Vierling-Claassen, P. Siekmeier, S. Stufflebeam and N. Kopell,
	Journal of Neurophysiology(99):2656-2671, 2008, doi:10.1152/jn.00870.2007"  
Bibliography:
  bibliography.bib

---

# Introduction
We provide an implementation of [@Vierling2008], which models impaired auditory entrainment
in the gamma range for schizophrenic patients. Particularly, we only reimplement the simplified network model
and do not replicate the Genesis model which is also developed in the article. 
We focus on the main result: an increase in inhibitory decay time constants leads to
a reduction of power in the gamma range and an increase in power in the beta range, replicating
experimental findings for schizophrenic patients. Therefore, we reproduce Figures 4,5,6, and 7 of the original paper.
The original model is implemented using Matlab but the source code is not publicly available.
The model and analysis scripts are implemented using Python 2.7.9.


# Methods

The methods section should explain how you replicated the original results:

* did you use paper description
* did you contact authors ?
* did you use original sources ?
* did you modify some parts ?
* etc.

The model was implemented solely from the paper description, since the original code is not publicly
available. The model is a simple model consisting of two neural populations (excitatory and inhibitory cells).
Indivudal cells are modeled as theta neurons (for a detailed discussion of this neuron model see [@Boergers2003]).
Each population connects to itself and to the other with an all-to-all connectivity. Both populations
also have two sources of input, the oscillatory drive input and a background noise input. The drive input
periodically sends spikes to both populations with a given frequency, whereas the background noise input
sends noise spikes at times drawn from a Poisson distribution. 

Since the theta neuron model, and the connectivity as well as the input of the network model are fairly simple, the implementation was straightforward.
The model was implemented using Python 2.7.9. using numpy ?.?.?. Visualisation of results was also done in Python
using the matplotlib module (matplotlib ?.?.?).
Furthermore, since the model is computationally very inexpensive, we did not aim to provide the most efficient implementation but rather
an implementation that is clear, and easy to understand and use.

\pagebreak

Table: Model summary {#tbl:summary}

+----------------+-----------------------------------------------------------------+
|Populations     |One excitatory and one inhibitory population|
+----------------+-----------------------------------------------------------------+
|Topology        |None|
+----------------+-----------------------------------------------------------------+
|Connectivity    |All-to-all|
+----------------+-----------------------------------------------------------------+
|Neuron model    |Theta model|
+----------------+-----------------------------------------------------------------+
|Synapse model   |(Quasi-)Instantaneous rise, exponential decay|
+----------------+-----------------------------------------------------------------+
|External input  |Poisson noise and periodic drive to both populations|
+----------------+-----------------------------------------------------------------+
|Recordings      |Theta variables (both populations); 'MEG' signal (summed EPSCs at exc. neurons)|
+----------------+-----------------------------------------------------------------+



Table: Model description {#tbl:descritpion}

-------------- -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Neuron model    $\frac{d \theta}{dt}=1-\cos \theta + (b+S+N(t))\cdot(1+ \cos \theta)$
                $\theta$ voltage variable; $b$ applied current; $S$ total synaptic input;$N(t)$ time-varying noise;
Synaptic input  $S_k = \sum_{j=1}^n \alpha_j \cdot g_{jk} \cdot s_{jk}$
                $\alpha_j=\pm 1$ controlling excitation and inhibition; $g_{jk}$ synaptic strength from cell $j$ to cell $k$; $s_{jk}$ synaptic gating variable from cell $j$ to cell $k$
Synapse model   $\frac{ds_{jk}}{dt}= - \frac{s _{jk}}{\tau _j} + e ^{- \eta \cdot (1+ \cos \theta _j)} \cdot \frac{1-s _{jk}}{\tau _R}$
                $\tau_j$ synaptic decay time; $\tau_R$ synaptic rise time
Drive model     Single excitatory pacemaker cell firing at the click
                train frequency and providing input to all cells
Noise model     $N=H(t-t_n) \cdot \frac{A \cdot g _{gmax} \cdot (e^{-(t-t _n)/ \tau _1} - e^{-(t-t _n)/ \tau _2} )}{\tau _1 - \tau _2}$
                EPSC for noise spike at time $t_n$
-------------- --------------------------------------------------------------------------------------------------------------------------------------------------------------------------



Table: Model parameters {#tbl:parameters}


Parameter               Definition                             Value 
----------------------- -------------------------------------- ------------------------------------------
$n_E$                   Exc. population size                   20
$n_I$                   Inh.population size                    10 
$\tau_R$                Synaptic rise time                     0.1 
$\tau_{exc}$            Excitatory decay time                  2.0 
$\tau_{inh}$            Inhibitory decay time (control)        8.0 
$\tau_{inh}$            Inhibitory decay time (schizophrenia)  28.0 
$g_{ee}$                E-E synaptic strength                  0.015 
$g_{ei}$                E-I synaptic strength                  0.025 
$g_{ie}$                I-E synaptic strength                  0.015 
$g_{ii}$                I-I synaptic strength                  0.02
$g_{de}$                Synaptic strength of drive to E cells  0.3 
$g_{di}$                Synaptic strength of drive to I cells  0.08
----------------------- -------------------------------------  ------------------------------------------








# Results

Reproduce figures 4-7 of the paper

![Replication of Figure 4: Raw simulated MEG signals (averaged over 20 trials) for the control and the schizophrenic network at the three different driving frequencies.](Replication-Figure4.eps){#fig:Vierling4}

![Replication of Figure 5: Power spectra of the averaged MEG signals from @fig:Vierling4](Replication-Figure5.eps){#fig:Vierling5}

![Replication of Figure 6: Single trial from the control network. 40 Hz drive. Network entrains to 40 Hz, as can be seen in the frequency diagram, raster plot and MEG trace.](Replication-Figure6.eps){#fig:Vierling6}

![Replication of Figure 7: Single trial from the schizophrenia network. 40 Hz drive. Network entrains to 40 Hz but also shows a strong 20 Hz component, as can be seen in the frequency diagram, raster plot and MEG trace. Especially the inhibitory neurons only entrain to a 20 Hz rhythm (see raster plot).](Replication-Figure7.eps){#fig:Vierling7}

# Conclusion

Conclusion, at the very minimum, should indicate very clearly if you were able
to replicate original results. If it was not possible but you found the reason
why (error in the original results), you should exlain it.


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
A reference to citation @markdown.

![Figure caption](rescience-logo.pdf){#fig:logo}

$$ A = \sqrt{\frac{B}{C}} $$ {#eq:1}


# References
