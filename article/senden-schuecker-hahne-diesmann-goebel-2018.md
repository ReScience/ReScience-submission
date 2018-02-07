---
Title: "A neural model of the saccade generator in the reticular formation"
Author:
  - name: Mario Senden
    affiliation: 1, 2
  - name: Jannis Schuecker,
    affiliation: 3
  - name: Jan Hahne,
    affiliation: 4
  - name: Markus Diesmann,
    affiliation: 3, 5, 6
  - name: Rainer Goebel,
    affiliation: 1, 2, 7
Address:
  - code:    1
    address: Department of Cognitive Neuroscience, Faculty of Psychology and Neuroscience,
			Maastricht University, 6201BC Maastricht, The Netherlands
  - code:    2
    address: Maastricht Brain Imaging Centre, Faculty of Psychology and Neuroscience, 
			Maastricht University, P.O. Box 616, 6200 MD Maastricht, The Netherlands
  - code:    3
    address: nstitute of Neuroscience and Medicine (INM-6) and Institute for Advanced Simulation (IAS-6) and JARA Institute Brain Structure-Function Relationships (INM-10), Jülich Research Centre, Jülich, Germany
  - code:    4
    address: School of Mathematics and Natural Sciences, Bergische Universit\"at Wuppertal,
			Wuppertal, Germany
  - code:    5
    address: Department of Psychiatry, Psychotherapy and Psychosomatics, Medical Faculty, 
             RWTH Aachen University, 52062 Aachen, Germany
  - code:    6
    address: Department of Physics, Faculty 1, RWTH Aachen University, 52062 Aachen, Germany
  - code:    7
    address: Department of Neuroimaging and Neuromodeling, Netherlands Institute for Neuroscience, 
			an Institute of the Royal Netherlands Academy of Arts and Sciences (KNAW), 1105BA Amsterdam, The Netherlands
Contact:
  - mario.senden@maastrichtuniversity.nl
Editor:
  - Name Surname
Reviewer:
  - Name Surname
  - Name Surname
Publication:
  received:  Month,  Day, 2018
  accepted:  Month, Day, 2018
  published: Month, Day, 2018
  volume:    "**1**"
  issue:     "**1**"
  date:      Month 2018
Repository:
  article:   "http://github.com/rescience/rescience-submission/article"
  code:      "http://github.com/rescience/rescience-submission/code"
  data:
  notebook:
Reproduction:
  - "*A neural model of the saccade generator in the reticular formation*, G. Gancarz,
	S. Grossberg, Neural Networks, 1159-1174, 1998"
Bibliography:
  bibliography.bib

---

# Introduction

We provide an implementation of the saccade generator (SG); a rate neuron model of neural circuitry in the reticular formation proposed by Gancarz & Grossberg [@Gancarz1998]. The SG is an integral part of the eye movement system [@Grossberg2012] and as such vital for developing large-scale architectures of visuo-motor integration. For this reason it is of interest to implement the model in a publicly available, widely used, and actively developed neural simulation framework such as NEST [@Gewaltig2007]. We show that the model translates well to the NEST framework as our implementation faithfully reproduces all simulation results reported in the original publication. Our code uses the Python interface [@Eppler2008] for legibility with both model and analysis scripts being implemented using Python 2.7.12.

# Methods

The SG model described by Gancarz & Grossberg [@Gancarz1998] consists of a horizontal and a vertical component each with two long-lead burst neurons (LLBNs), excitatory burst neurons (EBNs), inhibitory burst neurons (IBNs), and tonic neurons (TNs). Within each component, the two directions (left-right, down-up) interact antagonistically. Additionally, both components share a single omnipause neuron (OPN) which tonically inhibits each EBN to suppress unwanted saccades. In what follows we briefly describe the dynamics of neurons within the horizontal component (equivalent descriptions apply to the vertical component).

### Long-lead burst neurons

Each long-lead burst neuron ($L$) receives excitatory external input $I$ and inhibitory feedback from the ipsilateral IBNs ($B$):

$$
\begin{array}{ll}
\tau\frac{dL_l}{dt} &= -1.3L_l+I_l-2B_l \\\\
\tau\frac{dL_r}{dt} &= -1.3L_r+I_r-2B_r \\
\end{array}
$$

### Excitatory burst neurons

Each excitatory burst neuron ($E$) receives excitatory input from the ipsilateral LLBN, a constant arousal signal equal to $1$, and inhibitory input from the omnipause neuron ($P$):

$$
\begin{array}{ll}
\tau\frac{dE_l}{dt} &= -3.5E_l+(2-E_l)(5E_l+1)-(E_l+1)(10L_r+20g(P))) \\\\
\tau\frac{dE_r}{dt} &= -3.5E_r+(2-E_r)(5E_r+1)-(E_r+1)(10L_l+20g(P))) \\
\end{array}
$$
with a nonlinear gain function $g(x)$ given by

$$
g(x) = \frac{x^4}{0.1^4+x^4} \textrm{.}
$$

### Inhibitory burst neurons

Each inhibitory burst neuron ($B$) receives excitatory input from the ipsilateral EBN:

$$
\begin{array}{ll}
\tau\frac{dB_l}{dt} &= -2.4B_l+3E_l \\\\
\tau\frac{dB_r}{dt} &= -2.4B_r+3E_r \\
\end{array}
$$

### Omnipause neuron

The omnipause neuron ($P$) receives inhibitory input from all LLBNs:

$$
\tau\frac{dP}{dt} = -0.2P+(1-P)(1.2+J)-3.5(P+0.4)(g(L_l)+g(L_r)+g(L_u)+g(L_d))
$$
where $J$ represents external electrical stimulation.

### Tonic neurons

Each tonic neuron ($T$) receives excitatory input from the ipsilateral EBN and inhibitory input from the contralateral EBN:

$$
\begin{array}{ll}
\tau\frac{dT_l}{dt} &= 0.1(E_l-E_r) \\\\
\tau\frac{dT_r}{dt} &= 0.1(E_r-E_l) \textrm{.} \\
\end{array}
$$
Horizontal eye position ($\theta$) depends on activity of the right TN; i.e. $\theta = 260(T_r-0.5)$.

In implementing this model, we largely followed the descriptions provided in the original publication with a number of well-motivated exceptions. Specifically, the original model description has two features which cannot be straightforwardly translated to NEST. First, a nonlinear gain function is applied to a subset of inputs to EBNs and the OPN while a linear gain function is applied to their remaining inputs. Since NEST only applies a single gain function per neuron to each of its inputs, we opted for using a linear gain function for EBNs and the OPN. Furthermore, we passed those inputs requiring an additional nonlinear gain function through an auxiliary unit instantaneously applying the desired nonlinearity before passing the result on to EBNs and the OPN. Second, constant input to a neuron was not hard-coded but rather provided by an appropriately weighted bias node. Neither of these changes lead to discrepancies with original results.

In all simulations we used the Exponential Euler method for numerical integration of rate neurons [@Hahne2016] at a time step of $0.05\,\mathrm{ms}$ and a time constant of $50\,\mathrm{ms}$. It should be noted that rates of all neurons were initialized to zero and were thus not at resting equilibrium. In order to address this, the model was allowed to evolve for $100\,\mathrm{ms}$, to relax towards equilibrium before applying any input. Furthermore, we always simulated the full model; i.e. both its horizontal and vertical components even if input was applied only to one of the two.

# Results

All results from the original publication have been implemented. Our results accord very well with those reported by Gancarz & Grossberg [@Gancarz1998]. Only a single discrepancy was observed which we examined in depth.

The first simulation showcases the evolution of activity for each neuron type in the horizontal SG to a constant input applied to the left LLBN. The original publication does not report exact activation values observed for each neuron rendering a quantitative analysis of the accuracy of our replication impossible. However, qualitatively activation profiles shown in figure 3 in the original publication and those shown in figure @fig:fig_1 show good correspondence. In both implementations, the left LLBN showed a prelude of activity with left EBN bursts beginning after the onset of LLBN activity. The right EBN produced a small burst at the end of a saccade. Furthermore, increases in activity of the left TN were mirrored by decreases in the right. Finally, OPN activity dropped to zero during production of a saccade.

![**Activity profiles in the left (A) and right (B) SG.** 
All activities are in response to constant input ($\mathrm{I=1}$) applied to the left LLBN for $265\,\mathrm{ms}$.](figures/fig1.eps){#fig:fig_1 height="45%" width="75%"}

The second simulation shows the relation between input strength and burst amplitude for LLBNs and EBNs. With increasing input strength, both amplitude and duration of activation increased in long-lead and excitatory burst neurons. As before, these results accord well with those shown in figure 5 of the original publication. 

![**Activity profiles in LLBN (A) and EBN (B).** 
Increased input strength resulted in larger LLBN and EBN burst size with inputs equal to 1 (blue), 1.75 (green), and 2.5 (red) each applied to the left LLBN for $85\,\mathrm{ms}$.](figures/fig2.eps){#fig:fig_2 height="36%" width="75%"}

The third simulation generates saccades in response to different input strengths applied to the horizontal and vertical SG (figure @fig:fig_3). As in the original publication, saccades were generally straight with a slight tendency to curve. We extracted quantitative estimates of saccade amplitude for the original publication from their figure 6 and calculated the root-mean-squared error (RMSE) between these estimates and saccade amplitudes produced by our model. The RMSE was 0.16\textdegree\, 0.17\textdegree\, 0.18\textdegree\, 0.22\textdegree\, 0.14\textdegree\ for the blue, green, red, turqoiuse, and purple curves, respectively. Given that saccade amplitudes are between 10\textdegree\ and 15\textdegree\, these RMSE values indicate a close match between the two implementations.

![**Oblique saccades.** 
Inputs to the right and upward LLBNs were $\mathrm{{I}_{r}=0.67}$ & $\mathrm{{I}_{u}=0.08}$ (blue); $\mathrm{{I}_{r}=0.7}$ & $\mathrm{{I}_{u}=0.22}$ (green); $\mathrm{{I}_{r}=0.74}$ & $\mathrm{{I}_{u}=0.4}$ (red); $\mathrm{{I}_{r}=0.75}$ & $\mathrm{{I}_{u}=0.6}$ (turqoiuse); and $\mathrm{{I}_{r}=0.7}$ & $\mathrm{{I}_{u}=0.9}$ (purple) and were applied for $75\,\mathrm{ms}$. As expected, saccades were fairly straight with a slight tendency to curve.](figures/fig3.eps){#fig:fig_3 height="48.75%" width="48.75%"}

The fourth simulation implements a staircase of three saccades in response to continuous input (figure @fig:fig_4). Our results agree with those reported in the original publication (their figure 7) with saccades being of equal length and the smaller component (horizontal) being stretched to produce straight oblique saccades. 

![**Saccadic staircase.** 
Three saccades in a staircase continued in the same direction as the initial saccade. Inputs to the right and upward LLBNs were equal to $\mathrm{{I}_{r}=0.2}$ and $\mathrm{{I}_{u}=0.33}$ for $250\,\mathrm{ms}$. Eye position was sampled every $2\,\mathrm{ms}$.](figures/fig4.eps){#fig:fig_4 height="45%" width="33%"}

In the fifth simulation the average activity of the left EBN was obtained for a series of saccades with different directions. Figure @fig:fig_5 shows a polar plot of average activity corresponding to each saccade reflecting the neuron's tuning curve. The tuning curve we observed for the left EBN exhibits a cardioid-like shape as was the case in the original publication (figure 8).


				   0   45  72  90  108 135 162 180 198 225 252 270 288 315
------------------ --- --- --- --- --- --- --- --- --- --- --- --- --- ---
$\mathrm{{I}_{l}}$ .00 .00 .00 .00 .20 .45 .63 .70 .63 .45 .20 .00 .00 .00
$\mathrm{{I}_{r}}$ .70 .45 .20 .00 .00 .00 .00 .00 .00 .00 .00 .00 .20 .45
$\mathrm{{I}_{d}}$ .00 .00 .00 .00 .00 .00 .00 .00 .20 .45 .63 .70 .63 .45
$\mathrm{{I}_{u}}$ .00 .45 .63 .70 .63 .45 .20 .00 .00 .00 .00 .00 .00 .00
------------------ --- --- --- --- --- --- --- --- --- --- --- --- --- ---

Table: Direction specific inputs to SG to produce EBN tuning curve. {#tbl:table_1}

![**EBN tuning curve.** 
Tuning curve of the left EBN exhibiting a cardioid-like shape. Inputs to the SG producing the desired saccades can be found in table @tbl:table_1. Each of these inputs was applied for $50\,\mathrm{ms}$.](figures/fig5.eps){#fig:fig_5
height="37.5%" width="37.5%"}

\pagebreak
The sixth simulation reported in the original publication was designed to replicate results of Stanford \textit{et al.} [@Stanford1996]. These authors stimulated the superior colliculus (SC) at various frequencies and measured the resulting saccade amplitude, duration, and velocity; showing that ampltidude saturates before velocity. Our implementation of the SG model was capable of replicating these results. However, reproducing simulation results reported by Gancarz & Grossberg [@Gancarz1998] with our implementation was complicated by the fact that the stimulation protocol given by the authors lead to the production of two rather than a single saccade for larger stimulation intensities. Furthermore, the authors did not report their criteria for identifying saccade on- and offsets. Stanford \textit{et al.} [@Stanford1996] used velocity criteria to determine onset ($v\textgreater30\,\mathrm{deg/s}$) and offset ($v\textless30\,\mathrm{deg/s}$) of a saccade. While this provided us with explicit criteria, velocity did not always drop below $30\,\mathrm{deg/s}$ after the first saccade before rising again with the second. To determine the offset of a saccade in those cases, we found the local minimum between the end of the first and the beginning of the second saccade. With these criteria in place, we stimulated the SC. Results of our simulation are shown in figure @fig:fig_6. We extracted quantitative results of the original publication from their figure 9 and calculated the RMSE between these estimates and our results for each curve. The RMSE was 1.24\textdegree\, $0.46\,\mathrm{ms}$, and $4.45\,\mathrm{deg/s}$ for saccade amplitude, duration, and peak velocity, respectively. Our results accord thus very well with those reported in the original publication. This includes the observation that saccade amplitude and especially saccade duration were largest for a stimulation intensity of 1.2. It is conceivable that this intensity marks the point after which the SG generator starts producing two rather than a single saccade with the eye movement at this intensity merging two saccades and thus being stretched out.

![**Effect of stimulation frequency on saccade (A) amplitude, (B) duration, and (C) velocity.** 
Saccade amplitude (A) saturated before saccade velocity (C). A range of unitless values varying between 1 and 2.4 at increments of 0.2 were applied to the SC (reflecting stimulation at different frequencies). The connection weight from SC to LLBN was $\mathrm{W=2}$ and stimulation duration was $125\,\mathrm{ms}$.](figures/fig6.eps){#fig:fig_6 height="37.5%" width="75%"}

The seventh simulation shows that saccade velocity and duration can be traded while keeping amplitude constant. To produce a high-velocity saccade, the SC was stimulated at a high frequency. Similarly, to produce a low-velocity saccade, the SC was stimulated at a low frequency. Figure @fig:fig_7 shows the results of our simulation. In line with the original publication, the saccade amplitude reflected by TN activity was identical after high- and low-frequency saccades thus confirming the reported effect. However, we observe differences in the specific details of our implementation with respect to the original, starting with the shape of the input curves. While the rise of the input curves in the original publication and our implementation coincide well, the decay in our curves started later and proceeded with a shorter time constant ($\approx25\,\mathrm{ms}$, i.e. half the simulation time constant). This resulted in more total input to the model in our implementation and hence lead to the generation of a second saccade not observed in the original publication. When repeating our simulation using an ad-hoc solution to create matching input curves (by halving the time constant after external stimulation), we observed results identical to those reported in the original publication. We proceeded to investigate the origin of the discrepancy with regard to the shapes of the input curves. First, we tested whether discrepancies were specific to our NEST implementation by solving the expressions describing the input analytically. Specifically, according to equations A13-A15 of [@Gancarz1998] the input $I$ is described by
$$
\begin{array}{lll}
	\tau \frac{dA}{dt} &=-A+F(t)\\
	f(x) &=
	\left\{
	\begin{array}{lll}
		0 \quad &\textrm{if} \quad x<0 \\
		x \quad &\textrm{if} \quad 0<x<1 \\
		1 \quad &\textrm{if} \quad x \geq 1 \\
	\end{array}
	\right. \\
	I(t) &=Wf(A)
\end{array}
$$
\
with the stimulus
$$
F(t)=
\left\{
\begin{array}{lll}
	0 \quad &\textrm{if} \quad t<{t}_{on} \\
	F \quad &\textrm{if} \quad {t}_{on}<t<{t}_{off} \textrm{.} \\
	0 \quad &\textrm{if} \quad t \geq {t}_{off} \\
\end{array}
\right.
$$
\
Given this, the analytic expression for $I(t)$ is
$$
I(t)=
\left\{
\begin{array}{lll}
	0 \quad &\textrm{if} \quad t<{t}_{on} \\
	W \cdot f(F \cdot (1-{e}^{-(t-{t}_{on})/ \tau})) \quad &\textrm{if} \quad {t}_{on}<t<{t}_{off} \textrm{.}\\
	W \cdot f(F \cdot (1-{e}^{-({t}_{off}-{t}_{on})/ \tau}) \cdot {e}^{-(t-{t}_{off})/ \tau}) \quad &\textrm{if} \quad t \geq {t}_{off} \\
\end{array}
\right.
$$ {#eq:1}
Evaluating the analytic expression given parameter values reported in the original manuscript exactly reproduced our numerical results and thus produced curves equally deviatng from those shown in the original publication. Next, we contacted one of the original authors to rule out the possibility that we used erroneous parameter settings but no mistake was pointed out to us. The origin of the discrepancy thus remains unclear. However, it is likely it stems from the description of the input rather than from problems with the SG model itself.

![**Trading saccade velocity for duration.** 
Duration and velocity of a saccade can be traded while keeping amplitude constant. To produce a high velicity saccade an input of $\mathrm{F=3}$ was applied to the SC for $68\,\mathrm{ms}$ ($82\,\mathrm{ms}$ in the original publication). To produce a low velicity saccade an input of $\mathrm{F=1.3}$ for $117\,\mathrm{ms}$. Given the shape of the input curve produced by our implementation of the original equations (A), a second saccade can be observed for both high and low velocity saccades. If the input curves reported in the original manuscript were recreated in an ad-hoc fashion (B), the responses of all model neurons matched those in [@Gancarz1998]](figures/fig7.eps){#fig:fig_7
height="75%" width="75%"}

The eighth simulation reported by Gancarz & Grossberg [@Gancarz1998] shows how strong sustained input to the SG produces smooth eye movements. Our results, shown in figure @fig:fig_8, reproduced these findings as they strongly resembled those shown in figure 11 of the original publication. Specifically, an initial burst exhibited by the EBN was followed by sustained lower activity. This was due to inhibitory feedback being insufficient to silence the LLBN when input remains continuously strong and resulted, in turn, in the observed smooth eye movement.

![**Smooth staircase eye movements.** 
Activity profiles of SG neurons accompanying smooth eyes movement as a result of strong ($\mathrm{I=3}$) sustained ($300\,\mathrm{ms}$) input.](figures/fig8.eps){#fig:fig_8 height="46.2%" width="60%"}

![**Interrupted saccade resulting from OPN stimulation.** 
OPN stimulation interrupted the saccade, which remained accurate nonetheless. An  input of ($\mathrm{I=0.7}$) was applied to the LLBN for $100\,\mathrm{ms}$. At $45\,\mathrm{ms}$ after onset of the input, the OPN was stimulated ($\mathrm{J=1.8}$) for $5\,\mathrm{ms}$. The dashed line shows TN activity for an uninterrupted saccade.](figures/fig9.eps){#fig:fig_9 height="70%" width="46.2%"}

\clearpage

The final simulation showcases the evolution of activity exhibited by SG neurons when the OPN is briefly electrically stimulated while a constant input is applied to the LLBN. External stimulation temporarily restored activation in the OPN and hence also inhibition of the EBN, leading to an interruption of the saccade. As is shown in figure @fig:fig_9, the saccade remained accurate despite this disruption. This is in agreement with results shown in figure 12 of the original publication.

# Conclusion

The reproduced results show very good qualitative correspondence with those reported by Gancarz and Grossberg [@Gancarz1998] and accorded well quantitatively whenever such information was available. The only discrepancy between our and the original implementation was that the input to the saccade generator took more time to decay in our implementation. While the exact cause of this discrepancy is unclear, we could rule out that it was specific to the NEST implementation. Furthermore, it does not affect reproduction of the main finding of simulation seven; namely that saccade velocity and duration can be traded while keeping its amplitude constant.

In conclusion, our reproduction confirms the results of the original publication and shows that an implementation in the NEST framework is feasible. This allows for the straightforward integration of the saccade generator into models of the visuo-motor system forming the interface between sensory processing and motor control. 

# Acknowledgments 
All network simulations carried out with NEST (http://www.nest-simulator.org). This research was funded by the European Union’s Horizon 2020 Research and Innovation Programme under Grant Agreement No. 7202070 (HBP SGA1).

# References


