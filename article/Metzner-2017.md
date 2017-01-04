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
experimental findings for schizophrenic patients. Therefore, we reproduce Figs 4,5,6, and 7 of the original paper.
The original model is implemented using Matlab but the source code is not publicly available.
The model and analysis scripts are implemented using Python 2.7.9.


# Methods

The methods section should explain how you replicated the original results:

* did you use paper description
* did you contact authors ?
* did you use original sources ?
* did you modify some parts ?
* etc.


Table: Model summary {#tbl:table}

+----------------+-----------------------------------------------------------------+
|Populations     |One exc. and one inh.population|
+----------------+-----------------------------------------------------------------+
|Topology        |None|
+----------------+-----------------------------------------------------------------+
|Connectivity    |All-to-all|
+----------------+-----------------------------------------------------------------+
|Neuron model    |Theta model|
+----------------+-----------------------------------------------------------------+
|Synapse model   | Instantaneous rise (see [Boergers2003]), exponential decay|
+----------------+-----------------------------------------------------------------+
|Neuron model    |Theta model|
+----------------+-----------------------------------------------------------------+
|Neuron model    |Theta model|
+----------------+-----------------------------------------------------------------+


A reference to table @tbl:table.

If relevevant in your domain, you should also provide a new standardized
description of the work.


# Results

Results should be compared with original results and you have to explain why
you think they are the same or why they may differ (qualitative result vs
quantitative result). Note that it is not necessary to redo all the original
analysis of the results.


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
