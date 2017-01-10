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

This work serves as a reference implementation of the model proposed in [@purvis], which is a biophysically detailed model of the ionic currents involved in the rat hypoglossal motoneuron (HM). HMs innervate the tongue and are involved in breathing, and issues with HMs have been implicated in sleep apnea[@purvis]. A biophysical model is derived from known properties, and allows researchers to describe in detail the role and interactions of all components of the phenomenon. The model was originally solved (integrated) in XPP and analysed in MATLAB. This version has been developed in Python both for the fact that Python is relatively ubiquitous, and the availability of scientific packages (NumPy, SciPy, \dots). A version of this model exists in cellML format[@cellml], however the cellML model is not able to replicate the original results (action potentials do not occur at 1nA bias current), nor could it create action potentials with similar appearance. Due to the nature of the model, debugging is nearly impossible (multiple variables with the same name i.e. "theta_1").

# Methods
This work was prepared using the authors' description of the model [@purvis], and personal communication with the authors. Some slight inconsistencies were identified in the original paper which had a negligible effect on the results of the paper. These have been addressed in this implementation.

The mathematical description of the model was taken directly from page 2 of the original manuscript, and current dynamics were taken from the appendix. The constants (parameters) were taken from table 2 and the initial conditions were taken from table 1. It should be noted the the given initial conditions, while near the equilibrium point (rest), were not at equilibrium. In order to address this, the model was allowed to evolve (in the absence of bias current) until numerical equilibrium was reached, and this was taken to be the initial condition.

While the initial implementation used a Dormad-Prince adaptive solver, this implementation used the LSODA algorithm (provided in the scipy.integrate package), with absolute and relative tolerances of $10^{-10}$, and a maximum time step of $\Delta T=0.05$ ms (convergence of all results was tested with $\Delta T=0.01$ ms and no qualitative differences were found). The maximum time step is necessary in order to account for the relatively short step function current (not mentioned in the original paper).
This implementation used python 3.5.2 (64 bit) on Linux with packages: SciPy v0.17.1, NumPy v1.11.1, matplotlib v1.5.1.

# Results
All results from the original paper have been implemented, however some slight inconsistencies are noted and will be discussed later.

The first figure in the original paper is the voltage response to a 1 ms, 1 nA current pulse. This is reproduced in figure \ref{fig:fig1}. As seen, the results are very similar to the original figure, however the amplitude of the ADP seems a bit below that of the original figure. This discrepancy will be addressed in the discussion.

![Model response (action potential) to a 1 nA pulse active for 1 ms. Various properties are outlined in the response. \label{fig:fig1} ](../code/figures/figure_1.pdf)

The subsequent analysis performed in the original text is an analysis of the timescales and amplitudes of the individual ionic currents. This implementation is recreated in figure \ref{fig:fig2}. This implementation is very close to the original figure, however it is important to note that the I-SK current in this implementation has a peak of ~0.8 nA, whereas the peak amplitude in the original paper is ~1 nA, there are also slight differences in the calcium currents. 

![Ionic currents underlying an action potential created by a 1 nA pulse active for 1 ms. \label{fig:fig2}](../code/figures/figure_2.pdf)

Figure 3 in the text is an analysis of the dependency of the ADP on calcium and pre-pulse hyperpolarization. This has been replicated in figure \ref{fig:fig3}. Again there are slight discrepancies which will be addressed later.

![ADP dependancy on calcuim and voltage. Maximal T-type channel conductance is varied, and resting potential is clamped at a hyperpolarizing potential. \label{fig:fig3}](../code/figures/figure_3.pdf)

Validating the model against known experiments, in figure 4 in the text the authors simulate the application of apamin by eliminating the SK-current. This result is very well reproduced in figure \ref{fig:fig4}. The authors also demonstrate spike frequency adaptation which is known to occur with attenuated SK-current.

![Comparison of model to known experimental results. Removal of SK-current eliminates the ADP and attenuated SK-currents cause spike frequency adaptation. \label{fig:fig4}](../code/figures/figure_4.pdf)

Next, the model response to elongated current application is illustrated in figure \ref{fig:fig5}, this quite accurately reflects the results presented in the corresponding plot in the original paper (figure 5). The action potential plateau is accurately recreated and the length of the plateau is given by the length of the stimulus (neuron returns to rest after stimulus is removed). Discrepancies in the length of the plateau are purely due to discrepancies in the stimulus duration (the duration used in figure \ref{fig:fig5} is a guess that appears close).

![Model response to sustained 0.22 nA pulse which lasts for 2.3 seconds. \label{fig:fig5}](../code/figures/figure_5.pdf)

The next two figures are an analysis of frequency adaptation of the neurons under varying T-type calcium currents. In the text figure 7 is described as having a T-type current 10x that of figure 6. While true, it would be more accurate to say that figure 6 has a T-type calcium current 0.1x that of figure 7, as 0.1$\mu S$ is the "standard" conductivity and figure 6 has a conductivity of 0.01$\mu S$. These results are quite well replicated in figures \ref{fig:fig6} and \ref{fig:fig7}. Of note the "double peak" as the model stabilizes onto its sustained spiking behavior is accurately replicated.

![Model response to sustained current step with maximal T-type conductance of 0.01$\mu S$ \label{fig:fig6} ](../code/figures/figure_6.pdf)

![Model response to sustained current step with maximal T-type conductance of 0.1$\mu S$ \label{fig:fig7} ](../code/figures/figure_7.pdf)

Following this, spike adaptation is quantified as a function of both T-type and N-type calcium currents and varying the definition of total calcium current. These results are very well replicated in figure \ref{fig:fig8}.

![Spike frequency adaptation to various calcium currents conductances. \label{fig:fig8} ](../code/figures/figure_8.pdf)

The last figure in the text, figure 9, could not implemented solely based on the description given in the text. In the text the mAHP length is defined as being "measured from the repolarizing phase of the action potential to the point in the AHP that reaches the resting potential". Through clarification with the authors, mAHP length is defined as the length of time between the negative and positive crossing of the resting potential value.

This is implemented in figure \ref{fig:fig9}, this very closely reflects the original figure (digitized values provided in the code), however there are some values that deviate slightly, again this result will be discussed later.

![mAHP length as a funtion of H-current amplitude. \label{fig:fig9}](../code/figures/figure_9.pdf)

# Discussion

By and far, this implementation replicates the original paper. Slight discrepancies in the details however are apparent. For example in figure \ref{fig:fig1} the ADP amplitude is a few mV different in this implementation compared to the original. Also for example in figure \ref{fig:fig2} I-SK is ~0.8nA in comparison to ~1nA in the original. These discrepancies continue through several other figures. In order to track down these discrepancies the original authors provided the XPP files they used to create the figures. The XPP simulation had the same equations and parameters as this implementation and closely matches this implementation (code for figure 2 provides the option to compare the I-SK that XPP calculates with this implementation). Given the fact that the paper had the same parameters as this implementation and the XPP simulation matches this implementation, its likely that some simulations had a slightly different parameter set. For example figures \ref{fig:fig2} and \ref{fig:fig6} have different parameter sets, and there are discrepancies in figure \ref{fig:fig2} but not in figure \ref{fig:fig6}. This may imply that the discrepancies may be in the differing parameters.

This implementation has identified subtle discrepancies in the original paper. However this replication proves that all major results: replication of frequency adaptations, parametric variation, etc... are valid despite these discrepancies.

# Conclusion

With confidence, this implementation has indicated that the the original results were both true, and correctly implemented. This implementation did identify some subtle discrepancies in the original manuscript, which were addressed in this implementation.

# Acknowledgement
I would like to Dr. Butera for providing the XPP files for the original paper. Funding for this work was provided in the form of an OGS grant to ARS from the government of Ontario.

# References
