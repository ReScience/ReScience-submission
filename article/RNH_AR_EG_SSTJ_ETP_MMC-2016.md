---
Title: "Optimization of a free water elimination two-compartment model for diffusion tensor imaging."
Author:
  - name: Rafael Neto Henriques
    affiliation: 1
  - name: Ariel Rokem
    affiliation: 2
  - name: Eleftherios Garyfallidis
    affiliation: 3
  - name: Samuel St-Jean
    affiliation: 4
  - name: Eric Thomas Peterson
    affiliation: 5
  - name: Marta Morgado Correia
    affiliation: 1
Address:
  - code:    1
    address: MRC Cognition and Brain Sciences Unit, Cambridge, Cambridgeshire, UK
  - code:    2
    address: The University of Washington eScience Institute, Seattle, WA, USA
  - code: 3
    address: Indiana University School of Informatics and Computing, Indiana, IA, USA
  - code: 4
    address: University Medical Center Utrecht, Utrecht, NL
  - code: 5
    address: Biosciences, SRI International, Menlo Park, CA, USA
Contact:
  - rafaelnh21@gmail.com
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
  - "Hoy, A.R., Koay, C.G., Kecskemeti, S.R., Alexander, A.L. (2014).
Optimization of a free water elimination two-compartment model for diffusion
tensor imaging. NeuroImage 103, 323-333. doi: 10.1016/j.neuroimage.2014.09.053
"
Bibliography:
  bibliography.bib

---

# Introduction

Diffusion-weighted Magnetic Resonance Imaging (DW-MRI)  is a
biomedical imaging technique that allows for the non-invasive acquisition of *in vivo* data
from which tissue microstructure can be inferred. Diffusion tensor imaging (DTI), one of the most
commonly used DW-MRI techniques in the brain, models diffusion anisotropy of tissues using a
second-order tensor known as the diffusion tensor (DT) [@Basser1994-hg], [@Basser1994-zd].

DTI-based measures such as fractional anisotropy (FA) and mean diffusivity (MD)
are usually used to assess properties of brain microstructure. For example,
FA is thought to be an indicator of different microstructural properties
(such as packing density of axons, and density of myelin in nerve fibers [@Beaulieu2002-tl]),
but also indicates white matter coherence (i.e. the alignment of axons within a measurement voxel).
However, because a measurement voxel can contain partial volumes of different 
types of tissue, these measures are not always specific to one particular type of tissue.
In particular, diffusion anisotropy in voxels near cerebral ventricle and parenchyma can be
underestimated by partial volume effects of cerebrospinal fluid (CSF).

To remove the influence of the freely diffusing CSF and quantify changes that are specifically related to brain tissue, the DTI model
can be extended to separately take into account the contributions of tissue and CSF by representing
the tissue compartment with an anisotropic diffusion tensor and the CSF compartment as an isotropic
free water diffusion coefficient. Recently, two procedures were
proposed by Hoy and colleagues to fit this two compartments model to
diffusion-weighted data at different b-values (i.e. different diffusion gradient-weightings) [@Hoy2014-lk,].
Although these procedures have been shown to provide diffusion based measures that are stable
to different degrees of free water contamination, the authors noted that their
original algorithms were "implemented by a group member with no formal programming
training and without optimization for speed" [@Hoy2014-lk].

In this work, we provide the first open-source reference implementation of the
free water contamination DTI model. All implementations are made in Python based on the descriptions provided
in Hoy et al.'s original article. For speed optimization, all necessary standard
DT processing steps use previously optimized functions freely available with the software
package Diffusion Imaging in Python ([Dipy](http://nipy.org/dipy/),  [@Garyfallidis2014-zo])
and the optimization algorithms provided by the open-source software for mathematics,
science, and engineering [Scipy](http://scipy.org/).

# Methods

The free water elimination DTI (fwDTI) model describes the measured diffusion-weighted signal $s_i$
with a simple bi-exponential expansion of DTI:

$$ s_i = s_0 \left [ f \exp\left ( -b D_{iso} \right ) + 
      (1-f) \exp\left ( -b g_i D g_i \right )\right ] $$ {#eq:1} where $s_0$ is
the signal when no diffusion sensitization gradient is applied, $f$ is the volume fraction
of free water contamination, $D_{iso}$ is the free water isotropic diffusion coefficient which is set to $3.0 \times 10^{-3}  mm^{2}/s$,
$D$ is the diffusion tensor of the tissue, $b$ is a parameter that depends on the diffusion gradient shape, and $g_i$ is
the diffusion gradient direction.

##Implementation of fitting algorithms 
  
Since no source implementation was previously provided by Hoy and colleagues,
our implementation relies on the equations provided in the original article.



**Weighted-Linear Least Square (WLS).** Two errors were found in the equations describing the first proposed algorithm. Firstly, the free-water adjusted
diffusion-weighted signal formula (original article's Methods subsection
"FWE-DTI") should be written as: $$
y_{ik} = \ln\left\{ \frac{s_i - s_0 f_k\exp(-bD_{iso})}{(1-f_k)} \right \}$$ {#eq:2} instead of: $$
y_{ik} = \ln\left\{ \frac{s_i - s_0 \exp(-bD_{iso})}{(1-f_k)} \right \}$$ {#eq:3} where $f_k$ is a grid search free water volume fraction sample.

Secondly, according to the general linear least squares solution [@Jones2010-pg],
the parameters matrix is estimated using the weighted linear least squares solution
of the free-water elimination model:

$$\gamma = (W^TS^2W)^{-1}W^{T}S^{2}y$$ {#eq:4} instead of:
$$\gamma = (W^TS^2W)(SW)^{T}Sy$$ {#eq:5} where $\gamma$ contains the fwDTI model parameters
$\gamma=[D_{xx},D_{xy},D_{yy},D_{xz},D_{yz},D_{zz},\ln(s_0)]$,
$y$ is a matrix containing the elements of $y_{ik}$ computed from equation 2,
$S$ is a diagonal matrix with diagonal set to the $s_i$ samples, and $W$ is a matrix
computed from the $m$ diffusion-weighted directions $g_i$ and b-values:

$$
W =
\begin{bmatrix}
-b_1 g^2_{1x} & -2b_1 g_{1x}g_{1y} & -b_1 g^2_{1y} & -2b_1 g_{1x}g_{1z} & -2b_1 g_{1y}g_{1z} & -b_1 g^2_{1z} & 1      \\
\vdots        & \vdots             & \vdots        & \vdots             & \vdots             & \vdots        & \vdots \\
 b_m g^2_{mx} & -2b_m g_{mx}g_{my} & -b_m g^2_{my} & -2b_m g_{mx}g_{mz} & -2b_m g_{my}g_{mz} & -b_m g^2_{mz} & 1
\end{bmatrix}
$$ {#eq:6}

To ensure that the WLS method converges to the local minimum, $f$ grid search sampling is performed 
over larger interval ranges relative to original article. Particularly, for the second and third iterations used
to refine the parameters precision, $f$ is resampled over intervals of 0.2 and 0.02 instead of interval
sizes of 0.1 and 0.01 proposed by Hoy and colleagues. On the other hand, the sample step size was maintained to
0.01 and 0.001 respectively. 

Moreover, since the WLS objective function is sensitive to the squared error
of the model weights ($\omega_i=s_i$):

$$F_{WLS} = \frac{1}{2} \sum_{i=1}^{m} \left
 [ \omega_i \left ( y_{i} -\sum_{j=1}^{7}W_{ij}\gamma_{j}\right ) \right ]^{2}$$ {#eq:7} when evaluating which ($f$, $D_tissue$) pair is associated with smaller residuals the NLS
objective function is used instead:

$$F_{NLS} = \frac{1}{2} \sum_{i=1}^{m} \left
 [s_{i} - S_{0} f\exp(-b_iD_{iso})
- (1-f)\exp(-\sum_{j=1}^{7}W_{ij}\gamma_{j})\right ]^{2}$$ {#eq:8}

Similarly to the original article [@Hoy2014-lk], this procedure is only used
to obtain the intial guess for the free water elimination parameters, which were then
used to initialize a fwDTI model non-linear convergence solver (see below).


**Non-Linear Least Square Solution (NLS).** To improve computation speed, instead of using the modified Newton's algorithm proposed in the original article,
the non-linear convergence was done using Scipy's wrapped modified Levenberg-Marquardt algorithm
(function `scipy.optimize.leastsq` of [Scipy](http://scipy.org/)).

To constrain the model parameters to plausible ranges, some variable transformations
can be applied to the non-linear objective function. These were implemented as optional features
that can be controlled through user-provided arguments. To restrict the range of the volume fraction
to values between 0 and 1, the variable $f$ in equation 8 can be replaced by $\sin(f_t - \pi/2)/2+1/2$
and non-linear convergence is performed as a function of $f_t$. To ensure that the diffusion tensor is
positive definite, diffusion parameters can be converted to the Cholesky decomposition elements as described in [@Koay2006-zo].

In addition to the `scipy.optimize.leastsq` function, a more recently implemented version of Scipy's
optimization function `scipy.optimize.least_square` (available as of Scipy's version 0.17) was also tested.
The latter directly solves the non-linear problem with predefined
constraints in a similar fashion to what is done in the original article, however
our experiments showed that this procedure does not overcome the performance of
`scipy.optimize.leastsq` in terms of accuracy, and requires more computing time
(see supplementary_notebook_1.ipynb for more details).

To speed up the performance of the non-linear optimization procedure, the Jacobian of the free water
elimination DTI model was analytically derived and incorporated in the non-linear procedure (for details
of the Jacobian derivation see supplementary_notebook_2.ipynb). Due to increased mathematical complexity, 
our analytical Jacobian derivation is not compatible with the Cholesky decomposition. This
variable transformation is therefore not used by default.

**Removing problematic estimates.** For cases where the ground truth free water volume fraction is 1 (i.e. voxels
containing only free water), the tissue's diffusion tensor component can erroneously fit
the free water diffusion signal rather than placing the free water signal in the free water compartment,
and therefore incorrectly estimate the water volume fraction close to 0 rather than 1.
To remove these problematic cases, for all voxels with standard DTI mean diffusivity values larger than $2.7 \times 10^{-3} mm^{2}/s$, the free
water volume fraction is set to one while all other diffusion tensor
parameters are set to zero. This mean diffusivity threshold was arbitray adjusted to 90%
of the theoretical free water diffusion value, however this can be adjusted by
changing the optional input 'mdreg' in both WLS and NLS free water elimination
procedures.

**Implementation dependencies.** In addition to the dependency on Scipy, both
free water elimination fitting procedures require modules from Dipy [@Garyfallidis2014-zo],
since these contain all necessary standard diffusion tensor fitting functions.
Although the core algorithms for the free water elimination model are implemented here separately from Dipy,
a version of these will be incorporated as a sub-module of Dipy's model
reconstruction module ([https://github.com/nipy/dipy/pull/835](https://github.com/nipy/dipy/pull/835)). In addition, the
implemented procedures also requires the python pakage [NumPy](http://www.numpy.org/),
which is also a dependency of both Scipy and Dipy.

## Simulations
In their original study, Hoy and colleagues simulated a measurement along 32 diffusion
directions with diffusion weighting b-values of 500 and 1500 $s/mm^{2}$ and with six b-value=0 $s/mm^{2}$ images.
These simulations correspond to the results reported in Figure 5 of the original article.
We conducted Monte Carlo simulations using the multi-tensor simulation
module available in Dipy and using identical simulated acquisition parameters.
As in the original article, fitting procedures are tested for voxels with five different FA values
and with constant diffusion trace of $2.4 \times 10^{-3} mm^{2}/s$.
The eigenvalues used for the five FA levels are reported in Table @tbl:table.

Table: Eigenvalues values used for the simulations {#tbl:table}

FA            0                      0.11                   0.22                   0.3                    0.71
------------ ---------------------- ---------------------- ---------------------- ---------------------- ----------------------
$\lambda_1$  $8.00 \times 10^{-4}$  $9.00 \times 10^{-4}$  $1.00 \times 10^{-3}$  $1.08 \times 10^{-3}$  $1.60 \times 10^{-3}$
$\lambda_2$  $8.00 \times 10^{-4}$  $7.63 \times 10^{-4}$  $7.25 \times 10^{-4}$  $6.95 \times 10^{-4}$  $5.00 \times 10^{-4}$
$\lambda_3$  $8.00 \times 10^{-4}$  $7.38 \times 10^{-4}$  $6.75 \times 10^{-4}$  $6.25 \times 10^{-4}$  $3.00 \times 10^{-4}$

For each FA value, eleven different degrees of free water contamination were
evaluated (f values equally spaced from 0 to 1). To assess the robustness of the
procedure, Rician noise with signal-to-noise ratio (SNR) of 40 relative to the b-value=0 $s/mm^{2}$ images was
used. For each FA and f-value pair, simulations were performed for 120
different diffusion tensor orientations. Simulations for each diffusion tensor
orientation were repeated 100 times making a total of 12000 simulated
iterations for each FA and f-value pair.

## *In vivo* data

Similarly to the original article, the procedures are also tested using *in vivo* human brain data [@valabregue2015],
that can be automatically downloaded by Dipy's functions (see run_invivo_data.py code script).
The original dataset consisted of 74 volumes of images acquired for a
b-value of 0 $s/mm^{2}$ and 578 volumes diffusion weighted images acquired along 16 diffusion gradient directions
for b-values of 200 and 400 $s/mm^{2}$ and along 182 diffusion gradient directions for b-values
of 1000, 2000 and 3000 $s/mm^{2}$. In this study, only the data for b-values up to 2000 $s/mm^{2}$
are used to decrease the impact of non-Gaussian diffusion effects which are not
taken into account by the free water elimination model. We also processed the data with the standard DTI tensor model
(as implemented in Dipy) in order to compare the results with the free water elimination model.

# Results

The results from the Monte Carlo simulations are shown in Figure @fig:simulations. As reported
in the original article, FA values estimated using the free water elimination model match the tissue's ground truth
values for free water volume fractions $f$ ranging around 0 to 0.7 (top panel of
Figure @fig:simulations). However, FA values seem to be overestimated for higher volume fractions. This bias is more
prominent for lower FA values in which overestimations are visible from lower free water volume
fractions. The lower panels of Figure @fig:simulations suggest that the free water elimination model produces
accurate free water volume fraction for the full range of volume fraction ground truth values. All the features observed
here are consistent with Figure 5 of the original article.

![Fractional Anisotropy (FA) and free water volume fraction ($f$) estimates obtained with Monte Carlo simulations
using the free water elimination fitting procedures. The top panel shows the FA median and interquartile range
for the five different FA ground truth levels and plotted as a function of the ground truth water volume fraction.
The bottom panels show the estimated volume fraction $f$ median and interquartile range as a function of its ground truth values
(right and left panels correspond to the higher and lower FA values, respectively). This figure reproduces
Fig. 7 of the original article.](fwdti_simulations.png){#fig:simulations}

*In vivo* tensor statistics obtained from the free water elimination and standard DTI models
are shown in Figure @fig:invivo. Complete processing of all these measure took less than 1 hour
on an average Desktop and Laptop PC (~2GHz processor speed), while the reported processing time
by Hoy et al. was around 20 hours. The free water elimination model seems to produce higher values
of FA in general and lower values of MD relative to the metrics obtained from the standard DTI model.
These differences in FA and MD estimates are expected due to the suppression
of the isotropic diffusion of free water. As similarly reported in the original article,
high amplitudes of FA are observed in the periventricular gray matter which might be related to
inflated values in voxels with high $f$ values. These can be mitigated by excluding voxels with
high free water volume fraction estimates (see supplementary_notebook_3.ipynb), similarly to
what is suggested by Hoy and colleagues [@Hoy2014-lk].

![*In vivo* diffusion measures obtained from the free water DTI and standard
   DTI. The values of FA for the free water DTI model, the standard DTI model and
   their difference are shown in the top panels (A-C),
   while respective MD values are shown in the bottom panels (D-F). In addition
   the free water volume fraction estimated from the free water DTI model is shown in
   panel G.](In_vivo_free_water_DTI_and_standard_DTI_measures.png){#fig:invivo}


# Conclusion

Despite the changes done to reduce the algorithm's execution time, the
implemented procedures to solve the free water elimination DTI model have comparable performance
in terms of accuracy to the original methods described by Hoy and colleagues [@Hoy2014-lk].
Based on similar Monte Carlo simulations with the same SNR used in the original article,
our results confirmed that the free water elimination DTI model is able to remove confounding effects
of fast diffusion for typical FA values of brain white matter. Similarly to
what was reported by Hoy and colleagues, the proposed procedures seem to generate
biased values of FA for free water volume fractions near 1. Nevertheless,
our results confirm that these problematic cases correspond to regions that are not typically
of interest in neuroimaging analysis (voxels associated with cerebral ventricles)
and might be removed by excluding voxels with measured volume fractions above a reasonable
threshold such as 0.7.

# Author Contributions

Conceptualization: RNH, AR, MMC.
Data Curation: RNH, AR, EG, SSTJ.
Formal Analysis: RNH. 
Funding Acquisition: RNH, AR.
Investigation: RNH.
Methodology: RNH, AR, EG.
Project Administration: RNH, MMC, AR, EG.
Resources: RNH, MMC, AR.
Software: RNH, AR, ETP, EF, SSTJ.
Supervision: MMC, AR.
Validation: AR, SSTJ, EG.
Visualization: RNH.
Writing - Original Draft Preparation: RNH.
Writing - Review & Editing: RNH, AR, MMC, ETP, SSTJ.


# Acknowledgments

Rafael Neto Henriques was funded by Fundação para a Ciência e Tecnologia FCT/MCE (PIDDAC) under grant SFRH/BD/80114/2012.

Ariel Rokem was funded through a grant from the Gordon \& Betty Moore Foundation and the Alfred P. Sloan Foundation to the University of Washington eScience Institute.

Thanks to Romain Valabregue, CENIR, Paris for providing the data used here.


# References
