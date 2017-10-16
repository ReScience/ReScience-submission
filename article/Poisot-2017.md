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
  - "On the coexistence of specialists and generalists, David Sloan Wilson and Jin Yoshimura, The American Naturalist 144:4, 692-707, 1994."
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
\mathbf{K} = \begin{bmatrix}
  K_1 & aK_1 \\
  a_K2 & K_2 \\
  bK_1 & b_K2
\end{bmatrix} \,.
\end{equation}

In this matrix, $K_1$ is the quality of habitat 1, $K_2$ is the quality of
habitat 2, $a$ is the fitness cost of the specialist in its non-optimal
environment, and $b$ is the fitness cost of generalism. Note that $1 > b > a >
0$.

Species distribute themselves across habitats in a way that minimizes the
negative effect of other species on their fitness. This is modelled by each
species having a value $p_i$, which is the proportion of its species choosing
habitat 1. Values of $\mathbf{p}$ are found by measuring the negative density
effect of each species in each habitat:

\begin{equation}
D_{l1} = \frac{\sum_{i\in l,m,m}p_iN_i}{K_{l1}}
\end{equation}

and

\begin{equation}
D_{l2} = \frac{\sum_{i\in l,m,m}(1-p_i)N_i}{K_{l1}} \,.
\end{equation}

To find a value of $p_l$, we fix $p_m$ and $p_n$, and iterate over $p_l \in
[0;1]$, to find the value of $p_l$ minimizing $|D_{l1}-D_{l2}|$. This procedure
is repeated about 20 times.

Once the values of $\mathbf{p}$ are found, we can measure the density of
individuals in both habitats. Before we do so, there is a proportion $g$ of
individuals that select habitats at random. Given a total population size of
$N_i$, there are $N_i(g/2)$ individuals will go to either habitat, and
$N_i(1-g)p_i$ will pick habitat 1. With this information, we can write the
matrix describing habitat selection:

\begin{equation}
\mathbf{N} = \begin{bmatrix}
  N_1(\frac{g}{2}+(1-g)p_1) & N_1(\frac{g}{2}+(1-g)(1-p_1))\\
  N_2(\frac{g}{2}+(1-g)p_2) & N_2(\frac{g}{2}+(1-g)(1-p_2))\\
  N_3(\frac{g}{2}+(1-g)p_3) & N_3(\frac{g}{2}+(1-g)(1-p_3))
\end{bmatrix} \,.
\end{equation}


# Methods

For some non-stochastic situations, it is possible to calculate expected values
by hand. The original manuscript does provide some of these values, and they
were used to test this implementation. Original sources were not available, and
no attempts were made to contact the authors.

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
