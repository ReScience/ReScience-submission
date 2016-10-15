---
Title: "Ionic Current Model of a Hypoglossal Motoneuron"
Author:
  - name: Aaron R. Shifman
    affiliation: 1
Address:
  - code:    1
    address: Department of Biology, University of Ottawa, Ottawa, Ontario, Canada
Contact:
  - ashifman@uottawa.ca
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
  date:      Sep 2016
Repository:
  article:   "http://github.com/rescience/rescience-submission/article"
  code:      "http://github.com/rescience/rescience-submission/code"
  data:      
  notebook:  
Reproduction:
  - "*Ionic Current Model of a Hypoglossal Motoneuron*, Purvis LK., Butera RJ., (2005) Journal of Neurophysiology, 93, 723--733"
Bibliography:
  article.bib

---

# Introduction

This work serves as a reference implementation of the model proposed in [@purvis], which is a biophysical model of the rat hypoglossal motor neuron. The model was originally solved (integrated) in XPP and analysed in MATLAB. This version has been developed in python both for the fact that python is relatively ubiquitous, and availability of scientific packages (NumPy, SciPy, \dots). A version of this model exists in cellML format [here](https://models.cellml.org/workspace/purvis_butera_2005), however the cellML model is not able to replicate the original results (action potentials do not occur at 1nA bias current), nor could it create action potentials with similar appearance. Due to the nature of the model debugging is nearly impossible (multiple variables with the same name i.e. "theta_1").

# Methods
This work was prepared using the author's description of the model [@purvis], and personal communication with the authors. Some slight inconsistencies were identified in the original model which had negligible effect on the results of the paper. These have been corrected in this implementation.

The mathematical description of the model was taken directly from page 2, and current dynamics were taken from the appendix. The constants (parameters) were taken from table 2 and the initial conditions were taken from table 1. It should be noted the the given initial conditions were not the equilibrium point (2 significant display figures were insufficient). In order to address this, the model was allowed to evolve (in the absence of bias current) for 20 seconds, and the rest state (and hence initial conditions) was taken to be the solution after 20 seconds.

While the initial implementation used a Dormad-Prince adaptive solver, this implementation used the lsoda algorithm (provided in the scipy.integrate package), with absolute and relative tolerances of $10^{-10}$, and a maximum time step of $\Delta T=0.05$ ms (convergence of all results was tested with $\Delta T=0.01$ ms and no qualitative differences were found). The maximum time step is necessary in order to account for the discontinuity imposed by the step function current (not mentioned in the original paper).  
This implementation used python 3.5.2 (64 bit) on Linux with packages: SciPy v0.17.1, NumPy v1.11.1, matplotlib v1.5.1.

# Results
All results from the original paper were faithfully reproduced, however for space concerns only certain figures are included in the manuscript, although the included software will reproduce all original figures. 

The first figure in the original paper is the voltage response to a 1 ms, 1 nA current pulse. This is reproduced in figure \ref{fig:fig1}. As seen, the results are visually identical to the original figure, implying that the ODE solver, as well as all parameters, and functions were implemented correctly.

![Model response (action potential) to a 1 nA pulse active for 1 ms. Various properties are outlined in the response. Figure equivalent to figure 1 in the original text \label{fig:fig1} ](../code/figures/figure_1.pdf)

The subsequent analysis performed in the original text is an analysis of the timescales and amplitudes of the individual ionic currents (figure 2). This implementation faithfully replicates this analysis in figure \ref{fig:fig2}. However it is important to note that for the I-SK current the peak amplitude in this implementation is ~0.8 nA, whereas the peak amplitude in the original paper is ~1 nA. Since the shape of the time series is accurately reproduced, and the action potential appears to be nearly identical (figure \ref{fig:fig1}) it is most likely that the original figure had an error in the axis label.

![Ionic currents underlying an action potential created by a 1 nA pulse active for 1 ms. Figure equivalent to figure 2 in the original text \label{fig:fig2}](../code/figures/figure_2.pdf)

The model response to elongated current application is illustrated in figure \ref{fig:fig3}, this quite accurately reflects the results presented in the corresponding plot in the original paper (figure 5).

![Model response to sustained 0.22 nA pulse. Figure equivalent to figure 5 in the original text \label{fig:fig3}](../code/figures/figure_5.pdf)

The last figure in the text, figure 9, could not implemented solely based on the description given in the text. In the text the mAHP length is defined as being "measured from the repolarizing phase of the action potential to the point in the AHP that reaches the resting potential". Through clarification with the authors, mAHP length is defined as the length of time between the negative and positive crossing of the resting potential value.

This is implemented in figure \ref{fig:fig4}, this very closely reflects the original figure (digitized values provided in the code), however there are some values that deviate slightly. The most likely source of this discrepancy is some convergence error in the original figure. If not enough time was given for membrane potential to stabilize to its new value (given the new parameters), the definition of mAHP length would be slightly skewed. Given that this implementation explicitly gave sufficient time for the membrane potential to settle, and consequently creates a more parsimonious result (eliminates large discontinuity in trend), it is likely that the original implementation required a longer convergence.

![mAHP length as a funtion of H current amplitude. Figure equivalent to figure 9 in the original text \label{fig:fig4}](../code/figures/figure_9.pdf)

# Conclusion

With confidence, this implementation has indicated that the the original results were both true, and correctly implemented. This implementation did identify some subtle discrepancies in the original manuscript, which were addressed in this implementation.

# Acknowledgement
Funding for this work was provided in the form of an OGS grant to ARS from the government of Ontario

# References
