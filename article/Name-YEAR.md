---
Title: "This is the title"
Author:
  - name: Name Surname
    affiliation: 1
  - name: Name Surname,
    affiliation: 2, 3
Address:
  - code:    1
    address: Affiliation Dept/Program/Center, Institution Name, City, State, Country
  - code:    2
    address: Affiliation Dept/Program/Center, Institution Name, City, State, Country
  - code:    3
    address: Affiliation Dept/Program/Center, Institution Name, City, State, Country
Contact:
  - corresponding-author@mail.com
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
  number: 1
Repository:
  article:   "http://github.com/rescience/rescience-submission/article"
  code:      "http://github.com/rescience/rescience-submission/code"
  data:
  notebook:
Reproduction:
  - "Original article (title, authors, journal, doi)"
Bibliography:
  bibliography.bib

---

# Introduction

**The introduction should introduce the original paper and put it in context
(e.g. is it an important paper in the domain ?). You must also specify if there
was an implementation available somewhere and provide a link to it if relevant
(and in such a case, you have to specify if the proposed replication is based
on this original implementation). You should also introduce your implementation
by listing language, tools, libraries, etc. and motivate choices if relevant.**

Parasitism is a special case of predation. In both interactions, a species
(parasitoid or predator) feeds on the other species (host or prey), acting as a
regulating factor (@Anderson78). However, the population dynamics of both system
are very different. @Thompson24 was the first to propose a model to describe
this host-parasitoid system. In his model, parasites are limited by the number
of eggs they lay. Depending on the relative increase rate of hosts and
parasites, either both population increase indefinitely or decrease to
extinction. Later, @Nicholson35 proposed other models for which the rate of
increase of parasite is limited by their capacity to find hosts. These were the
basis for many other models where parasites act as regulating factors
(@Hassell78; @Rockwood15).

In 1983, Dempster proposed that natural enemies may not be an important
regulating factor in insect dynamics. In fact, he failed to detect
density-dependence due to natural enemies in most of the studies on Lepidoptera
he reviewed. His proposition really contrasted with what was thought at the
time. In response to this article, @Hassell85 analyzed a insect dynamic model in
which the only regulating factor was natural enemies. He showed that the
difficulties to detect the density-dependent effect of natural enemies was due
to time delays and stochasticity. This paper is still considered a classic in
fields of insect and parasitoid-host population dynamics. It introduced an
important argument on the role of natural enemies on insect populations, a
controversial topic that aroused ecologists to debate for almost a decade
(@Turchin90).

We used information from @Hassell85 to replicate the model. We were able to
replicate the results central to the article. In addition, we made new analysis
to detect the density-dependent effect of natural enemies in the stochastic
model. To our best knowledge, the original implementation was not available. The
code for the simulations and the figures were written in *Julia*.

# Methods

The formulas used in this paper to show the difficulty of detecting natural
enemies as regulating factors are the same that were used in the original paper
of @Hassell85 . First of all, the equation of the host population dynamics is
given as

$$ N(*t+1*) = F * N(*t*) * f(*Nt,Pt*) * D $$ {#eq:1}

where N(*t*) and N(*t+1*) represent the host population at generation *t* and at
the next generation, *F* is the rate of increase of the population and *D* is
the density independent probability of survival of the hosts (mortality). The
specialist parasitoids population dynamics are represented by

$$ P(*t+1*) = c * N(*t*) * [1-f(*Nt, Pt*)] $$ {#eq:2}

where P(*t*) and P(*t+1*) are the number of parasitoids at generation t* and at
the next one, while *c* is the number of female parasitoids emerging from each
host parasitized. In both Eq. 1 and Eq. 2, f(*Nt,Pt*) represents the probability
of escaping mortality from natural enemies (parasitoids) and is given by Eq. 3.

$$ f(*Nt,Pt*) = [1 + (a * P(*t*)) / (m * (1 + a * T*h* *N(*t*)))] ^-m $$ {#eq:3}

where *a* is the per capita searching efficiency of the parasitoids, *m* is the
extent of clumping of the parasitoids attacks and T*h* is the handling time as a
proportion of the total time. This paper also explores the relationship between
the hosts and a generalist parasitoid population. This population dynamic
follows the equation

$$ P(*t*) = h * (1 - exp^(-N(*t*)/b)) $$ {#eq:4}

where *h* is the saturation number of parasitoids and *b* is the rate of
approaching the saturation number.

In the previous equations, the parameters were assumed to remain constant
between successive generations. However, it might not be the case for every one
of them because they can be depending on the density of parasitoids and hosts.
The searching efficiency of the parasitoids is one of them. To calculate the
searching efficiency *A* at a generation *t*, the following equation is used :

$$ A = (1/P(*t*)) * ln(N(*t*)/S) $$ {#eq:5}

where S is the number of hosts that survived parasitism. Finally, to assess the
mortality linked to natural enemies (*k*), Eq. 6 will be used.

$$ k = Log(N(*t*)/S) $$ {#eq:6}

The objective in this paper is to reproduce every result from the original
publication. Every figure from Hassell's paper will be reproduced, except for
the Figure 2 and Figure 7, which only represent the curves from certain
equations when a parameter is changed. Therefore, they are not necessary in
order to show how difficult it is to detect the regulating effect of natural
enemies on a host population. For every figure reproduced, we will use the exact
same values that were used in the original paper for the different parameters.

The software used to run the models and to generate the figures is Julia version
0.6.2 (@Bezanson17) All the coding used to replicate the original paper will be
available with the article.


# Results

Results should be compared with original results and you have to explain why
you think they are the same or why they may differ (qualitative result vs
quantitative result). Note that it is not necessary to redo all the original
analysis of the results.

### Deterministic model

The reproduced population dynamics for the host population and for the
specialist parasitoids (Figures 3a and 3b) are very similar to the ones in
Hassell's paper. When the level of clumping of parasitoid attacks is high, both
populations decrease during the first 10 generations before they stabilise until
the end of the simulation, with the host population twice as big as the
specialist parasitoid population. When the extent of clumping is lower, both
populations show decreasing oscillations during the 50 simulated generations.
The difference between this article and the original publication resides in the
@fig:figure3 and Figure 3d, where we standardized the axes of the graphs. At first
sight, in Hassell's paper, the k-values in Figure 3c seem almost as high as in
Figure 3d. However, when using the same scale for the two graphs, we can see
that the k-values in the case of high clumping show much less variation than
when the lower clumping is used.

The simulations performed to study the relationship between the population
dynamics of the hosts and the generalist parasitoids are also very similar to
the ones done by Hassell. This show a density dependant relationship between the
two populations, where the natural enemy regulates the host population until
they both reach a stable equilibrium. This is also the conclusion when we
observe the relation between the two populations in Figures 3a and 3c.

### Stochastic model

# Discussion

We replicated all the results from @Hassell85. However (as expected) we did not
find the exact same results for the stochastic model. This can be explain by the
stochasticity in the model and, maybe, the way Hassell calculated the mortality
per generation (k-value). The latter was not explicitly explained in the
original paper, so we had to calculate it from our own interpretation. We added
two figures (Fig.7 and Fig.8) to take into account the variability in the
stochastic models and we had results that matches those from @Hassell85.

In all, we came up with the same conclusion that the original paper : the
density-dependent effect from natural enemies is obscured by time delay and/or
stochasticity.

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

![(a) and (b) Deterministic simulations of the host and specialist parasite population dynamics (Eq.1, Eq.2, Eq.3) using two different level of clumping in the parasitoid attacks: (a) m = 0.2; (b) m = 0.8. The other parameters used are the same in both (a) and (b): F = 4, D = 0.5, c = 1, a = 0.5 and T*h* = 0. (c) and (d) The relationship between the mortality caused each generation by parasitism (k-values) and the log10 host density for the fifty first generations, linked to (a) and (b) respectively.](figures/figure3.pdf){#fig:figure3}


![(a) Deterministic simulations of the host and generalist parasite population dynamics (Eq.1, Eq.3, Eq.4) with the following parameters: F = 4, D = 0.5, h = 10, b = 25, a = 0.5 and T*h* = 0 and m = 0.5. (b) The relationship between the mortality caused each generation by parasitism (k-values) and the log10 host density for the fifty first generations of the population dynamics in (a).](figures/figure4.pdf){#fig:figure4}


![(a) - (c) Deterministic simulations of the host and specialist parasite population dynamics (same as in Figure 3a) except for one parameter that is treated as a normally distributed stochastic variable: (a) D = 0.5 ± 0.5, (b) c = 0.5 ± 0.5 and (c) a = 0.5 ± 0.5. The other parameters are the same as in Figure 3a. (d) - (f) The relationship between the mortality caused each generation by parasitism (k-values) and the log10 host density for the fifty first generations, linked to (a), (b) and (c) respectively. The regression statistics for each relationship go as follows: (d) *y* = 0.184 + 0.211*x*; *r^2* = 0.171. (e) *y* = 0.329 + 0.054*x*; *r^2* = 0.254. (f) *y* = 0.190 + 0.158*x*; *r^2* = 0.120.](figures/figure5.pdf){#fig:figure5}

![(a) - (c) Deterministic simulations of the host and generalist parasite population dynamics (same as in Figure 4a) except for one parameter that is treated as a normally distributed stochastic variable: (a) D = 0.5 ± 0.5, (b) h = 10 ± 5 and (c) a = 0.5 ± 0.5. The other parameters are the same as in Figure 3a. (d) - (f) The relationship between the mortality caused each generation by parasitism (k-values) and the log10 host density for the fifty first generations, linked to (a), (b) and (c) respectively. The regression statistics for each relationship go as follows: (d) *y* = 0.146 + 0.192*x*; *r^2* = 0.243. (e) *y* = 0.306 + 0.015*x*; *r^2* = 0.214. (f) *y* = 0.149 + 0.160*x*; *r^2* = 0.107.](figures/figure6.pdf){#fig:figure6}

![(a) - (c) Distribution of the values of correlation coefficient (r^2) obtained for the 5000 simulations done in Figure 5(d) to Figure 5(f) respectively. The dotted line represents the value of r^2 obtained by Hassell in his original paper.](figures/figure7.pdf){#fig:figure7}

![(a) - (c) Distribution of the values of correlation coefficient (r^2) obtained for the 5000 simulations done in Figure 6(d) to Figure 6(f), respectively. The dotted line represents the value of r^2 obtained by Hassell in his original paper.](figures/figure8.pdf){#fig:figure8}


$$ A = \sqrt{\frac{B}{C}} $$ {#eq:1}

# References
