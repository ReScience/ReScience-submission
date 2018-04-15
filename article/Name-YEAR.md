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

$$N_{t+1} = F \times N_t \times f(N_t,P_t) \times D$$ {#eq:1}

where N(*t*) and N(*t+1*) represent the host population at generation *t* and at
the next generation, *F* is the rate of increase of the population and *D* is
the density independent probability of survival of the hosts (mortality). The
specialist parasitoids population dynamics are represented by

$$P_{t+1} = c \times N_t \times [1-f(N_t, P_t)]$$ {#eq:2}

where $P_t$ and $P_{t+1}$ are the number of parasitoids at generation t* and at
the next one, while *c* is the number of female parasitoids emerging from each
host parasitized. In both @eq:1 and @eq:2, $f(N_t,P_t)$ represents the
probability of escaping mortality from natural enemies (parasitoids) and is
given by @eq:3.

$$f(N_t,P_t) = [1 + (a \times P_t) / (m \times (1 + a \times T \times h \times N_t))]^{-m} $$ {#eq:3}

where $a$ is the per capita searching efficiency of the parasitoids, $m$ is the
extent of clumping of the parasitoids attacks and $T\times h$ is the handling
time as a proportion of the total time. This paper also explores the
relationship between the hosts and a generalist parasitoid population. This
population dynamic follows the equation

$$P_t = h \times \left(1 - \text{exp}\left(-\frac{N_t}{b}\right)\right)$$ {#eq:4}

where $h$ is the saturation number of parasitoids and $b$ is the rate of
approaching the saturation number.

In order to determine if the natural enemies can be declared as
density-dependent factors, the host population mortality $k$ will be plotted
against population density for each simulated generation. The correlation
coefficient $r$ of the resulting scatter plot will indicate the strength of the
density-dependence of natural enemies. The higher $r$ is, the strongest the
relation between hosts and parasites is. The host mortality is given by

$$k = \text{log}_{10}\frac{N_t}{S}$$ {#eq:5}

where $S$ is the number of hosts that survived parasitism. This number is given
by the host population density multiplied by the probability of escaping
mortality from natural enemies (@eq:6).

$$S = N_t \times f(N_t,P_t)$$ {#eq:6}

The objective in this paper is to reproduce every result from the original
publication. Every figure from the original paper will be reproduced, except for
the Figure 2 and Figure 7, which only represent the functions for some equations
when a parameter is changed. Therefore, they are not necessary in order to show
how difficult it is to detect the regulating effect of natural enemies on a host
population. For every figure reproduced, we will use the exact same values that
were used in the original paper for the different parameters.

The software used to run the models and to generate the figures is Julia version
0.6.2 (@Bezanson17) All the code used to replicate the original paper is
available alongside the article.

# Results

Results should be compared with original results and you have to explain why
you think they are the same or why they may differ (qualitative result vs
quantitative result). Note that it is not necessary to redo all the original
analysis of the results.

### Deterministic model

The reproduced population dynamics for the host population and for the
specialist parasitoids (@fig:figure3 (a) and (b)) are very similar to the ones in
Hassell's paper. When the level of clumping of parasitoid attacks is high, both
populations decrease during the first 10 generations before they stabilise until
the end of the simulation, with the host population twice as big as the
specialist parasitoid population. When the extent of clumping is lower, both
populations show decreasing oscillations during the 50 simulated generations.
The difference between this article and the original publication resides in the
@fig:figure3 (c) and (d), where we standardized the axes of the graphs. At first
sight, in Hassell's paper, the k-values in Figure 3c seem almost as high as in
Figure 3d. However, when using the same scale for the two graphs, we can see
that the k-values in the case of high clumping show much less variation than
when the lower clumping is used.

The simulations performed to study the relationship between the population
dynamics of the hosts and the generalist parasitoids are also very similar to
the ones done by Hassell (@fig:figure4). This shows a density dependant relationship between the
two populations, where the natural enemy regulates the host population until
they both reach a stable equilibrium. This is also the conclusion when we
observe the relation between the two populations in @fig:figure3 (a) and (c).

### Stochastic model

When the models include stochasticity, the resulting population dynamics between hosts and parasites are much more different than when the models are deterministic. Overall, the results obtained in the replications are pretty similar than the ones from Hassell, except for @fig:figure6 (d), where the correlation between the mortality and the host density is not nearly as strong as in the original paper. Whether we look at the model for the specialist parasitoids or the one for the generalist parasitoids, both show irregular oscillations in the host-parasitoid population dynamics (@fig:figure5 and @fig:figure6, (a), (b) and (c)). The addition of stochastic parameters prevents the stabilization of host populations and makes it more difficult to identify parasitoids as a density-dependent control factor. The relationship between the k-values and the host density are similar in Hassell's publication and in ours (except in @fig:figure6 (d)). The simulations done in this paper tend to have a determination coefficient (R^2) a little higher than the one found by Hassell, but they are generally really close.

Because the inclusion of stochastic parameters in the population dynamics causes variability in the outputs, the results from two successive simulations can be very different. In order to account for this variability and to show how it can affect the population dynamics of the hosts and parasites, @fig:figure7 and @fig:figure8 show the extent of the variation of the correlation coefficient ($r$) obtained in 5000 different simulations. A dotted line was added to represent the value of the correlation coefficient ($r$) that came out of the deterministic models. The values of r vary greatly for every stochastic parameter, but the one thing that stays constant is that in every case, the mean value of $r$ is lower in the stochastic models than in the deterministic models. This is in agreement with Hassell's results, and shows that stochasticity makes it harder to see the density dependence effect of the parasites, whether they are generalist or specialist.

# Discussion

We replicated all the results from @Hassell85. However (as expected) we did not
find the exact same results for the stochastic model. This can be explain by the
stochasticity in the model and, maybe, the way Hassell calculated the mortality
per generation (k-value). The latter was not explicitly explained in the
original paper, so we had to calculate it from our own interpretation. We added
two figures (@fig:figure7 and @fig:figure8) to take into account the variability in the
stochastic models and we had results that matches those from @Hassell85.

In all, we came up with the same conclusion that the original paper : the
density-dependent effect from natural enemies is obscured by time delay and/or
stochasticity.

# Conclusion

Conclusion, at the very minimum, should indicate very clearly if you were able
to replicate original results. If it was not possible but you found the reason
why (error in the original results), you should exlain it.

To sum up, we determine with the observation that, the distribution of correlation coefficient when stochasity present, is even lower than when there is no stochasticity, which show a weak density dependent relationship as we mentioned in the discussion. Besides this, We were able to replicate the results central to the article, and we draw the same conclusion as @Hassell85 said, with the present of stochasticity and time delay, the natural enemy as density dependent factor is considered to be undetectable from the analyse of existing life table.

![(a) and (b) Deterministic simulations of the host and specialist parasite population dynamics (Eq.1, Eq.2, Eq.3) using two different level of clumping in the parasitoid attacks: (a) m = 0.2; (b) m = 0.8. The other parameters used are the same in both (a) and (b): F = 4, D = 0.5, c = 1, a = 0.5 and T*h* = 0. (c) and (d) The relationship between the mortality caused each generation by parasitism (k-values) and the log10 host density for the fifty first generations, linked to (a) and (b) respectively.](figures/figure3.pdf){#fig:figure3}


![(a) Deterministic simulations of the host and generalist parasite population dynamics (Eq.1, Eq.3, Eq.4) with the following parameters: F = 4, D = 0.5, h = 10, b = 25, a = 0.5 and T*h* = 0 and m = 0.5. (b) The relationship between the mortality caused each generation by parasitism (k-values) and the log10 host density for the fifty first generations of the population dynamics in (a).](figures/figure4.pdf){#fig:figure4}


![(a) - (c) Deterministic simulations of the host and specialist parasite population dynamics (same as in Figure 3a) except for one parameter that is treated as a normally distributed stochastic variable: (a) D = 0.5 ± 0.5, (b) c = 0.5 ± 0.5 and (c) a = 0.5 ± 0.5. The other parameters are the same as in Figure 3a. (d) - (f) The relationship between the mortality caused each generation by parasitism (k-values) and the log10 host density for the fifty first generations, linked to (a), (b) and (c) respectively. The regression statistics for each relationship go as follows: (d) *y* = 0.184 + 0.211*x*; *r^2* = 0.171. (e) *y* = 0.329 + 0.054*x*; *r^2* = 0.254. (f) *y* = 0.190 + 0.158*x*; *r^2* = 0.120.](figures/figure5.pdf){#fig:figure5}

![(a) - (c) Deterministic simulations of the host and generalist parasite population dynamics (same as in Figure 4a) except for one parameter that is treated as a normally distributed stochastic variable: (a) D = 0.5 ± 0.5, (b) h = 10 ± 5 and (c) a = 0.5 ± 0.5. The other parameters are the same as in Figure 3a. (d) - (f) The relationship between the mortality caused each generation by parasitism (k-values) and the log10 host density for the fifty first generations, linked to (a), (b) and (c) respectively. The regression statistics for each relationship go as follows: (d) *y* = 0.146 + 0.192*x*; *r^2* = 0.243. (e) *y* = 0.306 + 0.015*x*; *r^2* = 0.214. (f) *y* = 0.149 + 0.160*x*; *r^2* = 0.107.](figures/figure6.pdf){#fig:figure6}

![(a) - (c) Distribution of the values of correlation coefficient ($r$) obtained for the 5000 simulations done in Figure 5(d) to Figure 5(f) respectively. The dotted line represents the value of $r$ obtained in the deterministic model.](figures/figure7.pdf){#fig:figure7}

![(a) - (c) Distribution of the values of correlation coefficient ($r$) obtained for the 5000 simulations done in Figure 6(d) to Figure 6(f), respectively. The dotted line represents the value of $r$ obtained in the deterministic model.](figures/figure8.pdf){#fig:figure8}

# References
