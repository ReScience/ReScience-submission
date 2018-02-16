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

Since the first describing of spike timing-dependent plasticity (STDP) [@Bi1998],
different description of STDP are published to reproduce different experimental findings.
The early implementations, so called pair-based STDP learning rules,
failed on reproducing some experimental observations,
like from triplet or quadruplets experiments [@Pfister2006].

Here, we introduce a reimplementation of the @Clopath2010 STDP model,
what is enable to reproduce experimental findings of triplet studies.
The proposed model try to be more biological plausible than previous models with
a theoretical approach, to be a voltage based STDP model.
This means, that the occur of long term depression (LTD) or long term potentiation (LTP) depends on the
postsynaptic membrane voltage.
Clopath and colleagues could show that their learning rule
can develop stable weights, as it is necessary for learning receptive fields of V1 simple cells.
To achieve stable weights, they implemented a homeostatic mechanism to
adjust the amount of generated LTD, over the relation between the average postsynaptic
membrane potential and a reference value.
Furthermore, their model lead to two different connection structures,
depending on the spiking behavior of the neurons [@Clopath2010].
If neurons fire high at the same time, they build strong bidirectional connections (rate code).
If neurons fire in a specific temporary order, they connection structure follow that order (temporal code).
They used a adaptive exponential integrate-and-fire (AdEx) neuron model.

Their model is reimplemented in Python (v2.7) and with help of the neuronal
simulator ANNarchy [@Vitay2015] (v4.6). For the analysis and the figures
we used numpy (v1.11.0) and matplotlib (v1.5.1). Not only the voltage based
STDP learing rule is reimplemented, even the AdEx neuron model.
In the supplementary material of the original publication is the matlab source code for a simple example
published. Besides them, it exists a matlab implementation to demonstrate the
stable learning of weights on modeldb (http://modeldb.yale.edu/144566).
The matlab implementation was used as a reference for the here presented one,
and used for evaluation of the correctness.

# Methods

The methods section should explain how you replicated the original results:

* did you use paper description
* did you contact authors ?
* did you use original sources ?
* did you modify some parts ?
* etc.

If relevant in your domain, you should also provide a new standardized
description of the work.


# Results

Results should be compared with original results and you have to explain why
you think they are the same or why they may differ (qualitative result vs
quantitative result). Note that it is not necessary to redo all the original
analysis of the results.


# Conclusion

Conclusion, at the very minimum, should indicate very clearly if you were able
to replicate original results. If it was not possible but you found the reason
why (error in the original results), you should explain it.


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

![Figure caption for Free!](rescience-logo.pdf){#fig:logo}

$$ A = \sqrt{\frac{B}{C}} $$ {#eq:1}


# References
