---
Title: "Insect natural enemies as regulating factors"
Author:
  - name: Dominique Caron
    affiliation: 1
  - name: Vincent Lessard
    affiliation: 1
  - name: Qile Wu
    affiliation: 1
  - name: Timothée Poisot
    affiliation: 1
Address:
  - code:    1
    address: Département de sciences biologiques, Université de Montréal, Montréal, Québec, Canada
  - code:    2
    address: Québec Centre for Biodiversity Sciences, Montréal, Québec, Canada
Contact:
  - timothee.poisot@umontreal.ca
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
  article:   "https://github.com/BIO6032/2018_replication_hassell_1985/article"
  code:      "https://github.com/BIO6032/2018_replication_hassell_1985/code"
  data:
  notebook:
Reproduction:
  - "Hassell, M. P. Insect Natural Enemies as Regulating Factors. Journal of Animal Ecology, vol. 54, no. 1, 1985, pp. 323–334."
Bibliography:
  bibliography.bib

---

# Introduction

Parasitism is a special case of predation. In both interactions, a species
(parasitoid or predator) feeds on the other species (host or prey), acting as a
regulating factor (@Anderson78). However, the population dynamics of both
systems are very different. @Thompson24 was the first to propose a model to
describe this host-parasitoid system. In his model, parasites are limited by the
number of eggs they lay. Depending on the relative increase rate of hosts and
parasites, either both population increase indefinitely or decrease to
extinction. Later, @Nicholson35 proposed other models for which the rate of
increase of the parasites is limited by their capacity to find hosts. These were
the basis for many other models where parasites act as regulating factors
(@Hassell78; @Rockwood15).

In 1983, @Dempster83 proposed that natural enemies may not be an important
regulating factor in insect dynamics. In fact, he failed to detect
density-dependence due to natural enemies in most of the studies on Lepidoptera
he reviewed. His proposition really contrasted with what was thought at the
time. In response to this article, @Hassell85 analyzed an insect dynamic model
in which the only regulating factor was natural enemies. He showed that the
difficulties to detect the density-dependent effect of natural enemies was due
to time delays and stochasticity. This paper is still considered a classic in
fields of insect and host-parasitoid population dynamics. It introduced an
important argument on the role of natural enemies on insect populations, a
controversial topic that aroused ecologists to debate for almost a decade
(@Turchin90).

We used information from @Hassell85 to replicate the model. We were able to
replicate the results central to the article. In addition, we made new analyses
for the stochastic models that bring additional support to @Hassell85 arguments.
To our best knowledge, the original implementation was not available. The code
for the simulations and the figures were written in *Julia* v0.6.2.

# Methods

The mathematical formulation used in this paper to show the difficulties of
detecting natural enemies as regulating factors are the same that were used in
the original paper by @Hassell85 . First of all, the host population dynamics
are given as

$$N_{t+1} = F \times N_t \times f(N_t,P_t) \times D$$ {#eq:1}

where $N_t$ and $N_{t+1}$ represent the host population at generation $t$
and at the next generation, $F$ is the rate of increase of the population and
$D$ is the density independent probability of survival of the hosts (mortality).
The specialist parasitoids population dynamics are represented by

$$P_{t+1} = c \times N_t \times [1-f(N_t, P_t)]$$ {#eq:2}

where $P_t$ and $P_{t+1}$ are the number of parasitoids at generation $t$ and at
the next one, while $c$ is the number of female parasitoids emerging from each
host parasitized. In both @eq:1 and @eq:2, $f(N_t,P_t)$ represents the
probability of escaping mortality from natural enemies (parasitoids) and is
given by @eq:3.

$$f(N_t,P_t) = [1 + (a \times P_t) / (m \times (1 + a \times T_h \times N_t))]^{-m} $$ {#eq:3}

where $a$ is the per capita searching efficiency of the parasitoids, $m$ is the
extent of clumping of the parasitoids attacks and $T_h$ is the handling
time as a proportion of the total time. This paper also explores the
relationship between the hosts and a generalist parasitoid population.
Generalist parasite dynamics follows the equation

$$P_t = h \times \left(1 - \text{exp}\left(-\frac{N_t}{b}\right)\right)$$ {#eq:4}

where $h$ is the saturation number of parasitoids and $b$ is the rate of
approaching this saturation number.

To determine if the natural enemies can be declared as density-dependent
factors, the host population mortality due to parasitism ($k-value$) is plotted
against population density for each simulated generation. The correlation
coefficient $r$ of the resulting scatter plot indicates the strength of the
density-dependence of natural enemies. The higher $r$ is, the strongest the
relation between hosts and parasites is. The host mortality is given by

$$k_\text{value} = \text{log}_{10}\frac{N_t}{S}$$ {#eq:5}

where $S$ is the number of hosts that survived parasitism. This number is given
by the host population density multiplied by the probability of escaping
mortality from natural enemies (@eq:6).

$$S = N_t \times f(N_t,P_t)$$ {#eq:6}

The objective in this paper is to reproduce the main results of the original
publication. Therefore, we did not reproduce figures 1, 2 and 7. Figure 1 was a
schematic representation of an insect life cycle. Figure 2 and Figure 7
represent relationship between parameters and population size/proportion of
parasited host. Therefore, they are not necessary in order to show how difficult
it is to detect the regulating effect of natural enemies on a host population.
For every figure reproduced, we used the exact same values that were used in
the original paper for the different parameters.

The software used to code and run the models and to generate the figures is
*Julia* version 0.6.2 (@Bezanson17). All the code used to replicate the original
paper is available alongside the article.

# Results

### Deterministic model

The reproduced population dynamics for the host population and for the
specialist parasitoids (@fig:figure3 (a) and (b)) are very similar to the ones
in Hassell's paper. When the level of clumping of parasitoid attacks is high,
both populations decrease during the first 10 generations before they stabilize,
with the host population twice as big as the specialist parasitoid population.
When the extent of clumping is lower, both populations show decreasing
oscillations during the 50 simulated generations. The difference between this
article and the original publication resides in the @fig:figure3 (c) and (d),
where we standardized the axes of the graphs. At first sight, in Hassell's
paper, the $k-values$ in @fig:figure3 (c) seem almost as high as in @fig:figure3 (d). However,
when using the same scale for the two graphs, we can see that the $k-values$ in
the case of high level of clumping show a lot less variation than with a lower
level of clumping.

The simulations performed to study the relationship between the population
dynamics of the hosts and the generalist parasitoids are also very similar to
the ones done by Hassell (@fig:figure4). This shows a density dependant
relationship between the two populations, where the natural enemy regulates the
host population until they both reach a stable equilibrium. This is also the
conclusion when we observe the relation between the two populations in
@fig:figure3 (a) and (c).

### Stochastic model

When the models include stochasticity, the resulting population dynamics between
hosts and parasites are very different compared to the deterministic models.
Overall, the results obtained in the replications are very similar to the
original paper. Whether we look at the model for the specialist parasitoids or the one for the
generalist parasitoids, both show irregular oscillations in the host-parasitoid
population dynamics (@fig:figure5 and @fig:figure6, (a), (b) and (c)). The
addition of stochastic parameters prevents the stabilization of host populations
and makes it more difficult to identify parasitoids as a density-dependent
control factor, except in @fig:figure6 (a). The relationships between the $k-values$ and the host density
are similar in Hassell's publication and in ours.
However, the regression for these relationships in our replication tend to have
a determination coefficient ($R^2$) higher than the one found by Hassell, but
they are generally really close. Discrepancies in the coefficient of
determination can be explained by different routines for pseudo random number
generation.

Moreover, the oscillations we obtained in @fig:figure6 (b) are not the same range compared to the original paper. More precisely, the oscillations we obtained with $h$ stochastic are smaller compared to the original paper (@fig:figure6 (b)).

Because the inclusion of stochastic parameters in the population dynamics causes
variability in the outputs, the results from two successive simulations can be
very different. In order to account for this variability and to show how it can
affect the population dynamics of the hosts and parasites, we added @fig:figure7
and @fig:figure8. These figures show the extent of the variation of the
correlation coefficient ($r$) obtained in 5000 different simulations (as opposed
to a single simulation in the original article). A dotted line was added to
represent the value of the correlation coefficient ($r$) that came out of the
deterministic models. The values of $r$ vary greatly for every stochastic
parameter. In every case, the mean value of $r$ is lower in the stochastic
models than in the deterministic models. This is in agreement with Hassell's
results, and shows that stochasticity makes it harder to detect the density
dependence effect of the parasites, whether they are generalist or specialist.

# Discussion

Overall, we were able to replicate most results from @Hassell85. We found the
exact same results for the deterministic model. We standardized the limits for
the axes, which was not the case in the original paper. This allows a more
convenient comparison of the different results.

As expected, we did not find the exact same dynamic for the stochastic model.
The figures we added (@fig:figure7 and @fig:figure8) showed how adding
stochasticity into the model can cause great variability in the output. For
example, in the specialist parasitoid model with a stochastic density-dependent
host survival ($D$), the correlation we found (@fig:figure7 (a)) was sometimes
very weak ($r \approx 0.2$) and some other times almost as strong as the
deterministic model ($r \approx 0.7$). Also, the correlation between the
mortality from parasitism ($k$-value) and host density ($N$) found in the
stochastic model was almost always weaker than in the deterministic model
(@fig:figure7 and @fig:figure8). Therefore, the results we added strongly
support the main argument from the original paper: adding stochasticity almost
always obscures the density-dependent effect of natural enemies.

The discrepancies we noted in the dynamics of the stochastic model with the
generalist parasite and a stochastic saturation number of parasitoid ($h$; in
@fig:figure6 (b)) are difficult to explain. It seems unlikely that it is caused
by either our implementation of the model, or by errors during runtime, since we
successfully replicated results from all the other numerical experiments.
Without the original implementation, we can only speculate on the difference
between the original implementation and ours. This combination of parameters is
the only one with a stochastic parameter normally distributed with a mean not
equal to its associated standard deviation. We tried with an $h$ normally
distributed with a mean and a standard deviation of 5 (@fig:figure9), and the
results were a lot more similar to the one from @Hassell85 than what we
originally had (@fig:figure6 (b)). Again, this is only an hypothesis on the kind
of error that could explain the differences between the original paper and ours.
Other possible sources of error include the pseudo-random numbers generator
used, errors in the original implementation, or errors in the parameter values
reported in the original publication; sadly, these are virtually impossible to
rule out.

The mathematical model from the original paper was well detailed, which allowed
us to create our own implementation. The equation for the number of survivors
from parasitism ($S$) was the only one we needed to deduce from our own
interpretation. This variable is used in the computation of the mortality
($k_\text{value}$) which is a well documented index. Therefore, this has not
limited us in the replication of the article, and the fact that the
deterministic simulations match these of the original paper suggests that we
used the same formulation for $S$.

# Conclusion

We were able to replicate the original results. Even if we did not find the
exact same dynamics for the stochastic models, we draw the same conclusions :
the density-dependent effect from natural enemies is obscured by time delays
and/or stochasticity. This makes it very difficult to detect natural enemies as
regulating factors from life table data. In addition, we added density plots for
the correlation coefficient from 5000 iterations of each stochastic model. This
allowed us to determine that the differences between our results and Hassell's
were explained by the stochasticity of the models. Also, these new results add a
strong support to the arguments of the original paper. To conclude, the
reproduction of the reference article @Hassell85 was successful and we hope it
adds to the legacy left by this significant paper in the history of population
dynamics.

![(a) and (b) Deterministic simulations of the host and specialist parasite population dynamics (Eq.1, Eq.2, Eq.3) using two different level of clumping in the parasitoid attacks: (a) m = 0.2; (b) m = 0.8. The other parameters used are the same in both (a) and (b): F = 4, D = 0.5, c = 1, a = 0.5 and T*h* = 0. (c) and (d) The relationship between the mortality caused each generation by parasitism (k-values) and the log10 host density for the fifty first generations, linked to (a) and (b) respectively.](figures/figure3.pdf){#fig:figure3}


![(a) Deterministic simulations of the host and generalist parasite population dynamics (Eq.1, Eq.3, Eq.4) with the following parameters: F = 4, D = 0.5, h = 10, b = 25, a = 0.5 and T*h* = 0 and m = 0.5. (b) The relationship between the mortality caused each generation by parasitism (k-values) and the log10 host density for the fifty first generations of the population dynamics in (a).](figures/figure4.pdf){#fig:figure4}


![(a) - (c) Deterministic simulations of the host and specialist parasite population dynamics (same as in Figure 3a) except for one parameter that is treated as a normally distributed stochastic variable: (a) D = 0.5 ± 0.5, (b) c = 0.5 ± 0.5 and (c) a = 0.5 ± 0.5. The other parameters are the same as in Figure 3a. (d) - (f) The relationship between the mortality caused each generation by parasitism (k-values) and the log10 host density for the fifty first generations, linked to (a), (b) and (c) respectively. The regression statistics for each relationship go as follows: (d) *y* = 0.049 + 0.148*x*; *r²* = 0.564. (e) *y* = 0.082 + 0.113*x*; *r²* = 0.034. (f) *y* = 0.083 + 0.158*x*; *r²* = 0.082.](figures/figure5.pdf){#fig:figure5}


![(a) - (c) Deterministic simulations of the host and generalist parasite population dynamics (same as in Figure 4a) except for one parameter that is treated as a normally distributed stochastic variable: (a) D = 0.5 ± 0.5, (b) h = 10 ± 5 and (c) a = 0.5 ± 0.5. The other parameters are the same as in Figure 3a. (d) - (f) The relationship between the mortality caused each generation by parasitism (k-values) and the log10 host density for the fifty first generations, linked to (a), (b) and (c) respectively. The regression statistics for each relationship go as follows: (d) *y* = 0.111 + 0.174*x*; *r²* = 0.843. (e) *y* = -0.010 + 0.307*x*; *r²* = 0.187. (f) *y* = -0.106 + 0.444*x*; *r²* = 0.262.](figures/figure6.pdf){#fig:figure6}


![(a) - (c) Distribution of the values of correlation coefficient ($r$) obtained for the 5000 simulations done in Figure 5(d) to Figure 5(f) respectively. The dotted line represents the value of $r$ obtained in the deterministic model.](figures/figure7.pdf){#fig:figure7}

![(a) - (c) Distribution of the values of correlation coefficient ($r$) obtained for the 5000 simulations done in Figure 6(d) to Figure 6(f), respectively. The dotted line represents the value of $r$ obtained in the deterministic model.](figures/figure8.pdf){#fig:figure8}


![(a) Deterministic simulations of the host and generalist parasite population dynamics (same as in Figure 6b) except that we used *h* = 5 ± 5. (b) The relationship between the mortality caused each generation by parasitism (k-values) and the log10 host density for the fifty first generations, linked to (a). The regression statistics go as follows: *y* = 0.030 + 0.192*x*; *r²* = 0.095.](figures/figure9.pdf){#fig:figure9}

# References
