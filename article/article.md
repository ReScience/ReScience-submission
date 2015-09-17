---
Title: "Chaos in a long-term experiment with a plankton community"
Author:
  - name: Owen Petchey
    affiliation: 1
  - name: Marco Plebani
    affiliation: 1
  - name: Frank Pennekamp
    affiliation: 1
Address:
  - code:    1
    address: Institute of Evolutionary Biology and Environmental Studies, University of Zurich, Zurich, Switzerland
Contact:
  - owen.petchey@ieu.uzh.ch
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
  - "Benincà, E., Huisman, J., Heerkloss, R., Jöhnk, K.D., Branco, P., Van Nes, E.H., Scheffer, M. & Ellner, S.P. (2008) Chaos in a long-term experiment with a plankton community. Nature, 451, 822–825 DOI: 10.1038/nature06512"
Bibliography:
  article.bib

---

# Introduction

The original paper describes analyses of fluctuations in the abundance of organisms in a plankton community derived from the Baltic Sea, housed in a laboratory environment. The length of the time series (samples every few days for 2,300 days) allowed for analyses revealing that the observed dynamics exhibited characteristics consistent with chaos produced by species interactions. The article concludes that stability is not required for persistence of complex food webs, and that long-term prediction of abundances may be fundamentally impossible. The demonstration of chaotic dynamics and limited forecast horizons (sensu @Petchey2015) are important in the field of ecology, since the ability to predict dynamics is an open question with considerable applied importance @Petchey2015 @Mouquet2015.


# Methods

This reproduction started with the raw data (source given below) and used information from the original paper, the supplementary information [the Supplement to the Nature paper](http://www.nature.com/nature/journal/v451/n7180/extref/nature06512-s1.pdf), and communications with Elisa Benincà and Stephen Ellner. The latter provided code used to produce results in the original paper, and its use in this reproduction is indicated below.

## Scope of the reproduction

An attempt was made to reproduce the majority of the results in the original article. Instances where we did not attempt to reproduce a result are detailed below.

## The data

The data are available as an Excel file supplement to [an Ecology Letters publication](http://onlinelibrary.wiley.com/doi/10.1111/j.1461-0248.2009.01391.x/abstract) @Beninca2009. The Excel file contains several datasheets. Two are particularly important, as they are the source of the raw data (one contains original species abundances, the one with the nutrient concentrations). We saved these two datasheets as comma separated value text files. In the code associated with this reproduction, these data files are read from the associated github repository.

Another datasheet in the Ecology Letters supplement contains transformed variables, though this was not used in this reproduction, as we reproduced the described transformations and applied them to the raw data.

We also received some data direct from Steve Ellner, see below for details. 

## Reproduction environment

The R language and environment for statistical computing and graphics was used to make the reproduction. Additional R packages required are specified in the code associated with this reproduction.

The code for this reproduction resides in an R markdown document, as well as a source file containing some required functions. Some of the code takes several minutes to run, so an intermediate data file is provided with results from this code.

# Results

## Population dynamics

The raw data show populations dynamics at least very similar to those in figure 1b-g of the original publication (figure @fig:dynamics).

![Observed population dynamics.](figures/unnamed-chunk-21-1) {#fig:dynamics}

## Data transformation

Did our transformation give data matching the data we received from Ellner, and that in the Ecology Letters supplement?

## Spectral analyses



## Correlations among species abundances



## Lyapunov exponents by direct method



## Lyapunov exponents by indirect method



## Predictability decay

Only done for non-linear models (the original article also does this for linear models).




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

![Figure caption](rescience-logo.pdf) {#fig:logo}

$$ A = \sqrt{\frac{B}{C}} $$ {#eq:1}


# Acknowledgements

This reproduction was made as part of the Reproducible Research in Ecology, Evolution, Behaviour, and Environmental Studies (RREEBES) Course, lead by Owen Petchey at the University of Zurich. More information about the course [here](https://github.com/opetchey/RREEBES/blob/master/README.md) on github.


# References
