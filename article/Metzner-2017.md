---
Title: "Modeling GABA Alterations in Schizophrenia: A Link Between Impaired Inhibition and Gamma and Beta Auditory Entrainment"
Author:
  - name: Metzner Christoph
    affiliation: 1
Address:
  - code:    1
    address: Centre of Computer Science and Informatics Research,University of Hatfield, Hatfield, Hertfordshire, UK

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
experimental findings for schizophrenic patients. Therefore, we reproduce Figures 4,5,6,7,10 and 11 of the original paper.
The original model is implemented using Matlab but the source code is not publicly available.
The model and analysis scripts are implemented using Python 2.7.9.


# Methods
The model was implemented solely from the paper description, since the original code is not publicly
available. The model is a simple model consisting of two neural populations (excitatory and inhibitory cells).
Indivudal cells are modeled as theta neurons (for a detailed discussion of this neuron model see [@Boergers2003]).
Each population connects to itself and to the other with an all-to-all connectivity. Both populations
also have two sources of input, the oscillatory drive input and a background noise input. The drive input
periodically sends spikes to both populations with a given frequency, whereas the background noise input
sends noise spikes at times drawn from a Poisson distribution. Table @tbl:summary summarizes the network model
whereas Table @tbl:description describes the equations underlying the different parts of the model. Tables
@tbl:modparameters and @tbl:simparameters list the parameters of the model and the simulations, their definitions and values, respectively.

The model was implemented using Python 2.7.9. using numpy ?.?.?. Visualisation of results was also done in Python
using the matplotlib module (matplotlib ?.?.?).
Furthermore, since the model is computationally very inexpensive, we did not aim to provide the most efficient implementation but rather
an implementation that is clear, and easy to understand and use.



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

+----------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|Neuron model    | $\frac{d \theta}{dt}=1-\cos \theta + (b+S+N(t))\cdot(1+ \cos \theta);$|
|                | $\theta$ voltage variable; $b$ applied current; $S$ total synaptic input;$N(t)$ time-varying noise;|
+----------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|Synaptic input  | $S_k = \sum_{j=1}^n \alpha_j \cdot g_{jk} \cdot s_{jk}$;|
|                | $\alpha _j=\pm 1$ controlling excitation and inhibition; $g_{jk}$ synaptic strength from cell $j$ to cell $k$; $s_{jk}$ synaptic gating variable from cell $j$ to cell $k$|
+----------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|Synapse model   | $\frac{ds_{jk}}{dt}= - \frac{s _{jk}}{\tau _j} + e ^{- \eta \cdot (1+ \cos \theta _j)} \cdot \frac{1-s _{jk}}{\tau _R}$;|
|                | $\tau _j$  synaptic decay time; $\tau_R$ synaptic rise time|
+----------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|Drive model     | Single excitatory pacemaker cell firing at the click|
|                | train frequency and providing input to all cells|
+----------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|Noise model     | $N=H(t-t_n) \cdot \frac{A \cdot g _{gmax} \cdot (e^{-(t-t _n)/ \tau _{exc}} - e^{-(t-t _n)/ \tau _R} )}{\tau _{exc} - \tau _R}$;|
|                | EPSC for noise spike at time $t_n$|
+----------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+



Table: Model parameters {#tbl:modparameters}


Parameter               Definition                                Value 
----------------------- ----------------------------------------- ------------------------------------------
$n_E$                   Exc. population size                      20
$n_I$                   Inh.population size                       10 
$\tau_R$                Synaptic rise time                        0.1 
$\tau_{exc}$            Excitatory decay time                     2.0 
$\tau_{inh}$            Inhibitory decay time (control)           8.0 
$\tau_{inh}$            Inhibitory decay time (schizophrenia)     28.0 
$g_{ee}$                E-E synaptic strength                     0.015 
$g_{ei}$                E-I synaptic strength                     0.025 
$g_{ie}$                I-E synaptic strength                     0.015 
$g_{ii}$                I-I synaptic strength                     0.02
$g_{de}$                Synaptic strength of drive to E cells     0.3 
$g_{di}$                Synaptic strength of drive to I cells     0.08
$b$                     Applied current (regardless of cell type) -0.1
$Ag_{max}$              Scaling factor for noise EPSCs            0.5
----------------------- ----------------------------------------- ------------------------------------------

Table: Simulation parameters {#tbl:simparameters}


Parameter              Value 
---------------------- ----------
Time step ($dt$)       0.061
Number of data points  8192 ($=2^{13}$)
Total time             500ms
---------------------- ----------


# Results

As explained in the introduction, we only replicated the simple model from [@Vierling2008], and **not** the GENESIS model.
We aimed to reproduce Figures 4 (raw, simulated MEG signal) and 5 (power spectra for MEG signals from Figure 4)  from [@Vierling2008], which show the main results of the model:

1. Control network:

..1.1. 40 Hz drive: Strong entrainment to the drive, no power in frequency bands apart from 40 Hz
 
..1.2. 30 Hz drive: Strong entrainment to the drive, no power in frequency bands apart from 30 Hz

..1.3. 20 Hz drive: Entrainment to the drive, however, more power in the harmonic, 40 Hz band than in the 20 Hz band

2. Schizophrenia network:

..2.1. 40 Hz drive: Weaker entrainment to the drive (compared to control network), emergence of a subharmonic component (at 40 Hz)

..2.2. 30 Hz drive: Strong entrainment to the drive, no power in frequency bands apart from 30 Hz

..2.3. 20 Hz drive: Stronger entrainment to the drive (compared to control network), less power in the harmonic band


Figures @fig:Vierling4 and @fig:Vierling5 show the output of the replicated model for the same simulations as for Figures 4 and 5 from the original article. The main charactersitics
described above can be clearly seen. However, in our model these main features are a little bit less pronounced than in the original model. Since the network
model receives Poissonian noise (which is quite strong), this difference may simply stem from a difference in noise.
Furthermore, we have to mention the differences in amplitude for the simulated MEG signals (and accordingly the power spectra thereof) between the original model
and our replication, which we believe stems from a scaling in the original model. However, since the original source code is not available, we cannot verify that.

![Replication of Figure 4: Raw simulated MEG signals (averaged over 20 trials) for the control and the schizophrenic network at the three different driving frequencies.](Replication-Figure4.eps){#fig:Vierling4}

![Replication of Figure 5: Power spectra of the averaged MEG signals from @fig:Vierling4](Replication-Figure5.eps){#fig:Vierling5}

After having looked at the model output averaged over 20 trials with different noise, we also show single trial data which exemplify the main model features
, as was done in the original article.
Figures @fig:Vierling6 and fig:Vierling7 show the model output in response to 40 Hz drive input for the control and the schizophrenia network, respectively. 
The strong entrainment in the control case, the reduction of entrainment and the emergence of a subharmonic 20 Hz component are again clearly visible. However, as before,
in our model implementation the emergent 20 Hz component is less pronounced than in the original implementation (best seen in the excitatory population activitydisplayed in the raster plot of @fig:Vierling7).

![Replication of Figure 6: Single trial from the control network. 40 Hz drive. Network entrains to 40 Hz, as can be seen in the frequency diagram, raster plot and MEG trace.](Replication-Figure6.eps){#fig:Vierling6}

![Replication of Figure 7: Single trial from the schizophrenia network. 40 Hz drive. Network entrains to 40 Hz but also shows a strong 20 Hz component, as can be seen in the frequency diagram, raster plot and MEG trace. Especially the inhibitory neurons only entrain to a 20 Hz rhythm (see raster plot).](Replication-Figure7.eps){#fig:Vierling7}

Figures @fig:Vierling10 and @fig:Vierling11 show the model output in response to 20 Hz drive for the control and the schizophrenia network, respectively. Again, main features of the original model
are faithfully reproduced.

![Replication of Figure 10: Single trial from the control network. 20 Hz drive. Network entrains to 20 Hz but also shows a 40 Hz component, as can be seen in the frequency diagram, raster plot and MEG trace.](Replication-Figure10.eps){#fig:Vierling10}

![Replication of Figure 11: Single trial from the schizophrenia network. 20 Hz drive. Network entrains to 20 Hz without 40 Hz component, as can be seen in the frequency diagram, raster plot and MEG trace. note that the 40 Hz power in the power spectrum is a harmonic.](Replication-Figure11.eps){#fig:Vierling11}



# Conclusion

Overall, we believe we have faithfully reproduced the main fetaures of the simple model from [@Vierling2008]. However, we note that overall features a re a little bit less pronounced in
our reimplementation compared to the original model. This may simply stem from a difference in noise between the two sets of simulations.




# References