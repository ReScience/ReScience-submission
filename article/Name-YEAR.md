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

Parasitism is a special case of predation. In both interactions, a species (parasitoid or predator) feeds on the other species (host or prey), acting as a regulating factor (Anderson and May, 1978). However, the population dynamics of both system are very different. Thompson (1924) was the first to propose a model to describe this host-parasitoid system. In his model, parasites are limited by the number of eggs they lay. Depending on the relative increase rate of hosts and parasites, either both population increase indefinitely or decrease to extinction. Later, Nicholson and Bailey (1935) proposed other models for which the rate of increase of parasite is limited by their capacity to find hosts. These were the basis for many other models where parasites act as regulating factors (Hassell, 1978; Rockwook, 2015). 

In 1983, Dempster proposed that natural enemies may not be an important regulating factor in insect dynamics. In fact, when reviewing data from 24 studies on Lepidoptera, he failed, in most cases, to detect density-dependence due to natural enemies. That view challenged what was thought at that time. In response to this article, Hassell (1985) analyzed a insect dynamic model in which the only regulating factor was natural enemies. He showed that the difficulties to the density-dependent effect of natural enemies was due to time delays and stochasticity. This paper is still considered a classic in fields of insect and parasitoid-host population dynamics. **AJOUTER DES EXEMPLES DE CONTRIBUTION DE L'ARTICLE**

We used information from Hassel (1985) to replicate the model. We recreated the figures central to the article and increased the number of simulations for the stochastic model. To our best knowledge, the original implementation was not available. The code for the simulations and the figures were written in *Julia*.


# Methods

The methods section should explain how you replicated the original results:

* did you use paper description
* did you contact authors ?
* did you use original sources ?
* did you modify some parts ?
* etc.

If relevevant in your domain, you should also provide a new standardized
description of the work.

The formulas used in this paper to show the difficulty of detecting natural enemies as regulating factors are the same that were used in the original paper (Hassell, 1985). First of all, the equation of the host population dynamics is given as

$$ N(*t+1*) = F * N(*t*) * f(*Nt,Pt*) * D $$ {#eq:1}

where N(*t*) and N(*t+1*) represent the host population at generation *t* and at the next generation, *F* is the rate of increase of the population and *D* is the density independent probability of survival of the hosts (mortality). The specialist parasitoids population dynamics are represented by 

$$ P(*t+1*) = c * N(*t*) * [1-f(*Nt, Pt*)] $$ {#eq:2}

where P(*t*) and P(*t+1*) are the number of parasitoids at generation t* and at the next one, while *c* is the number of female parasitoids emerging from each host parasitized. In both Eq. 1 and Eq. 2, f(*Nt,Pt*) represents the probability of escaping mortality from natural enemies (parasitoids) and is given by Eq. 3.

$$ f(*Nt,Pt*) = [1 + (a * P(*t*)) / (m * (1 + a * T*h* *N(*t*)))] ^-m $$ {#eq:3}

where *a* is the per capita searching efficiency of the parasitoids, *m* is the extent of clumping of the parasitoids attacks and T*h* is the handling time as a proportion of the total time. This paper also explores the relationship between the hosts and a generalist parasitoid population. This population dynamic follows the equation

$$ P(*t*) = h * (1 - exp^(-N(*t*)/b)) $$ {#eq:4}

where *h* is the saturation number of parasitoids and *b* is the rate of approaching the saturation number.

In the previous equations, the parameters were assumed to remain constant between successive generations. However, it might not be the case for every one of them because they can be depending on the density of parasitoids and hosts. The searching efficiency of the parasitoids is one of them. To calculate the searching efficiency *A* at a generation *t*, the following equation is used :

$$ A = (1/P(*t*)) * ln(N(*t*)/S) $$ {#eq:5}

where S is the number of hosts that survived parasitism. Finally, to assess the mortality linked to natural enemies (*k*), Eq. 6 will be used.

$$ k = Log(N(*t*)/S) $$ {#eq:6}



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
