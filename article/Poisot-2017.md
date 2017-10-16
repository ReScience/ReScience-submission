---
Title: "On the coexistence of specialists and generalists"
Author:
  - name: Timothée Poisot
    affiliation: 1
Address:
  - code:    1
    address: Département de Sciences Biologiques, Université de Montréal, Montréal, QC, Canada
Contact:
  - timothee.poisot@umontreal.ca
Editor:
  - Name Surname
Reviewer:
  - Name Surname
  - Name Surname
Publication:
  received:  Sep, 1, 2015
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
  - "On the coexistence of specialists and generalists, David Sloan Wilson \& Jin Yoshimura, The American Naturalist 144:4, 692-707, 1994."
Bibliography:
  bibliography.bib
---

# Introduction

The coexistence of specialists and generalist within ecological communities is a
long-standing question. @wils94csg have suggested that this coexistence can be
understood when examined in the light of (i) differential fitness loss
associated to specialism, (ii) active habitat selection, (iii) negative density
dependence due to competition, and (iv) stochastic changes in habitat quality,
that allow combinations of species to persist even though coexistence would not
be possible in a purely deterministic world. Here I propose an implementation of
this model in *Julia* [@beza17jfa], and show that it is able to reproduce most
figures from the original manuscript.

The @wils94csg model describes three species across two patches of habitat.
Species 1 is a specialist of habitat 1, species 2 is a specialist of habitat 2,
and species 3 is a generalist. This results in the maximum density that these
species can reach in both habitats:

\begin{equation}
\mathbf{K} = \begin{bmatrix} K_1 & aK_1 \\ a_K2 & K_2 \\ bK_1 & b_K2
\end{bmatrix} \,.
\end{equation}

In this matrix, $K_1$ is the quality of habitat 1, $K_2$ is the quality of
habitat 2, $a$ is the fitness cost of the specialist in its non-optimal
environment, and $b$ is the fitness cost of generalism. Note that $1 > b > a >
0$.

# Methods

The methods section should explain how you replicated the original results:

* did you use paper description
* did you contact authors ?
* did you use original sources ?
* did you modify some parts ?
* etc.

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
